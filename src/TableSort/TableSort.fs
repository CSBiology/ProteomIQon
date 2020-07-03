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

    let aggregationMethod (agMethod: Domain.AggregationMethod): Series<'K0,'V1> -> float =
        match agMethod with
        | Domain.AggregationMethod.Sum    -> Stats.sum
        | Domain.AggregationMethod.Mean   -> Stats.mean
        | Domain.AggregationMethod.Median -> Stats.median

    let aggregationMethodPepToProt (agMethod: Domain.AggregationMethod): float[] -> float =
        match agMethod with
        | Domain.AggregationMethod.Sum ->
            fun (x: float[]) ->
                if x.Length = 0 then nan
                else Array.sum x
        | Domain.AggregationMethod.Mean ->
            fun (x: float[]) ->
                if x.Length = 0 then nan
                else Array.average x
        | Domain.AggregationMethod.Median ->
            fun (x: float[]) ->
                if x.Length = 0 then nan
                else Array.median x

    let getAggregatedPeptidesVals (peptidesPresent: Set<string>) (peptidesMapped:Series<string,string[]>) (data: Frame<string,string>) (columnName:string) (agMethod: Domain.AggregationMethod): Series<string,float> =
        peptidesMapped
        |> Series.mapValues (fun peptides ->
            peptides
            |> Array.map (fun peptide ->
                if peptidesPresent.Contains(peptide) then
                    Some (try Some (data.Item(columnName).Item(peptide)) with _ -> None)
                else
                    None
            )
        )
        |> Series.mapValues (fun x ->
            x
            |> Array.choose id
            |> Array.choose id
            |> (aggregationMethodPepToProt agMethod))

    let applyLevelWithException (levelSel:_ -> 'K) (ex: 'C[]) (op:_ -> 'T) (exOp:_ -> 'T) (frame:Frame<'R, 'C>) =
        let indexBuilder = Deedle.Indices.Linear.LinearIndexBuilder.Instance
        let vectorBuilder = Deedle.Vectors.ArrayVector.ArrayVectorBuilder.Instance
        frame.GetColumns<'T>()
        |> Series.map (fun c s ->
                if (Array.contains c ex) then
                    Series.applyLevel levelSel exOp s
                else
                    Series.applyLevel levelSel op s)
        |> FrameUtils.fromColumns indexBuilder vectorBuilder

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
                let currentField = filterOn.[i]
                let frame' =
                    printfn "1"
                    frameAcc
                    |> Frame.filterRowValues (fun s ->
                        if currentField.FieldName = "Class" then
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

    let filterFrameTukey (tukeyFields: (string*float)[]) (frame: Frame<'a,string>) =
        let rec tukeyLoop (i: int) (frameAcc: Frame<'a,string>) =
            if i = tukeyFields.Length then
                frameAcc
            else
                let tukeyC = snd tukeyFields.[i]
                let field = fst tukeyFields.[i]
                let vals = Frame.getCol field frame
                let borders =
                    vals
                    |> Series.values
                    |> Array.ofSeq
                    |> fun x -> FSharp.Stats.Testing.Outliers.tukey tukeyC x
                let frame' =
                    frameAcc
                    |> Frame.filterRowValues (fun s ->
                        match s.TryGetAs<float>(field) with
                        | OptionalValue.Missing -> true
                        | OptionalValue.Present v -> v <  borders.Upper && v > borders.Lower)
                tukeyLoop (i+1) frame'
        tukeyLoop 0 frame

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

        let isQuantColumnOfInterest = 
            let s = quantColumnsOfInterest |> Set.ofArray
            fun x -> s.Contains x

        let tables =
            Array.map2 ( fun quantFile protFile ->
                // reads quant and prot table
                let quantTable: Frame<int,string> =
                    Frame.ReadCsv(path=quantFile ,separators=param.SeparatorIn, schema=quantColumnTypes)
                let protTable : Frame<string,string> =
                    Frame.ReadCsv(path=protFile ,separators=param.SeparatorIn)
                    |> Frame.indexRowsString "GroupOfProteinIDs"
                let quantTableFiltered: Frame<int,string> =
                    // calculate 14N/15N ratios based on the corresponding fields from the quant table
                    // Ratio gets added before filtering, so that it can also be filtered on
                    if param.Labeled then
                        let ratios =
                            let n14 = quantTable.GetColumn<float>"Quant_Light"
                            let n15 = quantTable.GetColumn<float>"Quant_Heavy"
                            n14/n15
                        quantTable.AddColumn ("Ratio", ratios)

                    quantTable
                    |> filterFrame param.QuantFieldsToFilterOn
                    // is it okay to apply tukey outlier detection after filtering?
                    |> filterFrameTukey param.TukeyQuant
                let protTableFiltered: Frame<string,string> =
                    protTable
                    |> filterFrame param.ProtFieldsToFilterOn
                    |> filterFrameTukey param.TukeyProt
                let quantTableAggregated =
                    quantTableFiltered
                    // group all rows based on the peptide sequence. This groups different charges and global modifications of the same peptide
                    //|> Frame.groupRowsUsing (fun k ser -> (ser.GetAs<string>("StringSequence")))
                    |> Frame.groupRowsBy "StringSequence"
                    // makes sure only columns with values are present (i.e. charge/global mod are removed)
                    |> Frame.filterCols (fun ck _ -> isQuantColumnOfInterest ck)
                    // aggregates the columns over the peptide sequence with a defined method (i.e. average,median,...)
                    |> applyLevelWithException (fun (sequence,index) -> sequence) [|"Quant_Light";"Quant_Heavy"|]
                        (aggregationMethod param.AggregatorFunction) (aggregationMethod param.AggregatorFunctionIntensity)
                    // drop sparse rows to prevent missing 14n/15n quant values, which will later complicate aggregation of peptides to proteins
                    |> Frame.dropSparseRows
                    //|> Frame.dropCol "PEPValue"
                // set of every peptide present in the quant table based on the row keys
                let peptidesPresent = quantTableAggregated.RowKeys |> Set.ofSeq
                // table containing every peptide that is mapped to a protein group. The row keys are the protein groups and the PeptideSequences column contains a string array with all peptides mapping to them.
                let peptidesMapped =
                    protTableFiltered
                    |> Frame.getCol "PeptideSequence"
                    |> Series.mapValues (fun (s : string)-> s.Split(';'))
                // calculates count of peptides mapped to each protein group
                let distinctPeptideCount =
                    protTableFiltered
                    |> Frame.getCol "PeptideSequence"
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
                    baseTable.AddColumn (name, getAggregatedPeptidesVals peptidesPresent peptidesMapped quantTableAggregated name param.AggregatorPepToProt)
                ) |> ignore
                baseTable
                |> Frame.mapRowKeys (fun key -> key, System.IO.Path.GetFileNameWithoutExtension protFile)

            ) quantFiles protFiles
        tables
        |> Frame.mergeAll
        |> fun frame ->
            frame.SaveCsv (path=(outDirectory+(@"\TableSort.tab")), separator=(param.SeparatorOut))
            frame
            |> Frame.pivotTable
                (fun (proteinGroup, experiment) os -> proteinGroup)
                (fun (proteinGroup, experiment) os -> experiment)
                (fun f ->
                    f .GetRowAt<float> 0
                )
                |> Frame.expandAllCols 1
            |> fun framePiv -> framePiv.SaveCsv (path=(outDirectory+(@"\TableSort_horizontal.tab")), separator=param.SeparatorOut)