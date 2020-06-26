namespace ProteomIQon

open Deedle
open FSharpAux.IO
open FSharpAux
open FSharp.Stats
open BioFSharp

module TableSort =

    let createSchema (fieldType: string) (fields: string[]) =
        fields
        |> Array.map (fun s -> sprintf "%s=%s"s fieldType)
        |> String.concat ","

    let dropCols (columns: string[]) (frame: Frame<'a,string>): Frame<'a,string> =
        let rec loop i (frameAcc: Frame<'a,string>) =
            if i < columns.Length then
                let newFrame =
                    frameAcc
                    |> Frame.dropCol columns.[i]
                loop (i + 1) newFrame
            else
                frameAcc
        loop 0 frame

    let getAggregatedPeptidesVals (peptidesPresent: Set<string>) (peptidesMapped:Series<string,string[]>) (data: Frame<string,string>) (columnName:string): Series<string,float> =
        peptidesMapped
        |> Series.mapValues (fun peptides ->
            peptides
            |> Array.map (fun peptide ->
                if peptidesPresent.Contains(peptide) then
                    Some (data.Item(columnName).Item(peptide))
                else
                    None
            )
        )
        |> Series.mapValues (fun x ->
            x
            |> Array.choose id
            |> Array.median)

    let peptideEvidenceClassToFloat (ec: string) =
        match ec with
        | "Unknown" -> 0.
        | "C1a" -> 1.
        | "C1b" -> 2.
        | "C2a" -> 3.
        | "C2b" -> 4.
        | "C3a" -> 5.
        | "C3b" -> 6.
        | _ -> failwith "invalid peptide evidence class"

    let filterFrame (filterOn: Domain.FilterOnField[]) (frame: Frame<'a,string>) =
        let rec filterLoop (i: int) (frameAcc: Frame<'a,string>) =
            if i = filterOn.Length then
                printfn "2"
                frameAcc
            else
                let frame' =
                    let currentField = filterOn.[i]
                    printfn "1"
                    frameAcc
                    |> Frame.filterRowValues (fun s ->
                        if currentField.FieldName = "EvidenceClass" then
                            match s.TryGetAs<string>(currentField.FieldName) with
                            | OptionalValue.Missing -> failwith "EvidenceClass can not be missing"
                            | OptionalValue.Present v ->
                                match (currentField.UpperBound, currentField.LowerBound) with
                                | None, None -> true
                                | Some uB, None -> (peptideEvidenceClassToFloat v) < uB
                                | None, Some lB -> (peptideEvidenceClassToFloat v) > lB
                                | Some uB, Some lB ->
                                    (peptideEvidenceClassToFloat v) < uB && (peptideEvidenceClassToFloat v) > lB
                        else
                            match s.TryGetAs<float>(currentField.FieldName) with
                            | OptionalValue.Missing -> true
                            | OptionalValue.Present v ->
                                match (currentField.UpperBound, currentField.LowerBound) with
                                | None, None -> true
                                | Some uB, None -> v < uB
                                | None, Some lB -> v > lB
                                | Some uB, Some lB -> v < uB && v > lB
                    )
                printfn "3"
                filterLoop (i+1) frame'
        filterLoop 0 frame

    let sortTables (quantFiles: string[]) (protFiles: string[]) outDirectory (param: Domain.TableSortParams) =

        // creates a schema that sets the column type of every column containing values to float
        let quantColsWithValues =
            Array.append param.QuantColumnsOfInterest (param.QuantFieldsToFilterOn |> Array.map (fun x -> x.FieldName))
            |> Array.distinct
        let quantColumnTypes =
            createSchema "float" quantColsWithValues
        let quantColumnsOfInterest =
            if param.Labeled then 
                param.QuantColumnsOfInterest
                // appends "Ratio" column since it's not included in the original columns with values
                |> Array.append [|"Ratio"|]
            else
                param.QuantColumnsOfInterest
        let tables =
            Array.map2 ( fun quantFile protFile ->
                // reads quant and prot table
                let quantTable: Frame<int,string> =
                    Frame.ReadCsv(path=quantFile ,separators="\t", schema=quantColumnTypes)
                let protTable : Frame<string,string> =
                    Frame.ReadCsv(path=protFile ,separators="\t")
                    |> Frame.indexRowsString "ProteinID"
                let quantTableFiltered: Frame<int,string> =
                    quantTable
                    |> filterFrame param.QuantFieldsToFilterOn
                let protTableFiltered: Frame<string,string> =
                    protTable
                    |> filterFrame param.ProtFieldsToFilterOn
                let quantTableAggregated =
                    // calculate 14N/15N ratios based on the corresponding fields from the quant table
                    if param.Labeled then
                        let ratios =
                            let n14 = quantTableFiltered.GetColumn<float>"N14Quant"
                            let n15 = quantTableFiltered.GetColumn<float>"N15Quant"
                            n14/n15
                        quantTableFiltered.AddColumn ("Ratio", ratios)

                    quantTableFiltered
                    // group all rows based on the peptide sequence. This groups different charges and global modifications of the same peptide
                    |> Frame.groupRowsUsing (fun k ser -> (ser.GetAs<string>("StringSequence")))
                    // makes sure only columns with values are present (i.e. charge/global mod are removed)
                    |> Frame.filterCols (fun ck _ -> 
                        (
                            quantColumnsOfInterest
                            |> Set.ofArray
                        ).Contains(ck)
                    )
                    // aggregates the columns over the peptide sequence with a defined method (i.e. average,median,...)
                    |> Frame.applyLevel (fun (sequence,index) -> sequence) Stats.mean
                    //|> Frame.dropCol "PEPValue"
                // set of every peptide present in the quant table based on the row keys
                let peptidesPresent = quantTableAggregated.RowKeys |> Set.ofSeq
                // table containing every peptide that is mapped to a protein group. The row keys are the protein groups and the PeptideSequences column contains a string array with all peptides mapping to them.
                let peptidesMapped =
                    protTableFiltered
                    |> Frame.getCol "PeptideSequences"
                    |> Series.mapValues (fun (s : string)-> s.Split(';'))
                // calculates count of peptides mapped to each protein group
                let distinctPeptideCount =
                    protTableFiltered
                    |> Frame.getCol "PeptideSequences"
                    |> Series.map (fun _ x ->
                        x
                        |> String.split ';'
                        |> Array.length
                    )
                // an "empty" table based on the prot table. It only contains the row keys with the protein groups.
                let baseTable =
                    protTableFiltered
                    |> Frame.addCol "DistinctPeptideCount" distinctPeptideCount
                    |> Frame.filterCols (fun ck _ -> 
                        (
                            param.ProtColumnsOfInterest
                            |> Set.ofArray
                        ).Contains(ck)
                    )
                // the columns in the array are selected in the filtered quant table and compared to the protein groups. If the peptide is mapping to a protein group, then
                // an entry is added to the base table at the corresponding protein group.
                quantColumnsOfInterest
                |> Array.map (fun name ->
                    baseTable.AddColumn (name, getAggregatedPeptidesVals peptidesPresent peptidesMapped quantTableAggregated name)
                ) |> ignore
                baseTable
                |> Frame.mapRowKeys (fun key -> key, System.IO.Path.GetFileNameWithoutExtension protFile)

            ) quantFiles protFiles
        tables
        |> Frame.mergeAll
        |> Frame.pivotTable
            (fun (proteinGroup, experiment) os -> proteinGroup)
            (fun (proteinGroup, experiment) os -> experiment)
            (fun f ->
                f .GetRowAt<float> 0
            )
            |> Frame.expandAllCols 1
        |> fun frame -> frame.SaveCsv (path=(outDirectory+(@"\TableSort.tab")), separator='\t')