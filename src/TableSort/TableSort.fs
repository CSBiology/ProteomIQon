namespace ProteomIQon

open Deedle
open FSharpAux.IO
open FSharpAux
open FSharp.Stats

module TableSort =

    let neededQuantColsWithValues = [|"PEPValue";"N14Quant";"N15Quant"|]

    let neededQuantColsIdentifiers = [|"StringSequence";"GlobalMod";"Charge"|]

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

    let sortTables (quantFiles: string[]) (protFiles: string[]) outDirectory (param: Domain.TableSortParams) =

        // creates a schema that sets the column type of every column containing values to float
        let quantColumnTypes =
            createSchema "float" neededQuantColsWithValues
        let tables =
            Array.map2 ( fun quantFile protFile ->
                // reads quant and prot table
                let quantTable: Frame<int,string> =
                    Frame.ReadCsv(path=quantFile ,separators="\t", schema=quantColumnTypes)
                let protTable : Frame<string,string> =
                    Frame.ReadCsv(path=protFile ,separators="\t")
                    |> Frame.indexRowsString "ProteinID"
                // set of all columns from the quant table that are not needed for the final table
                let quantColsNotNeeded =
                    Set.difference (quantTable.ColumnKeys |> Set.ofSeq) ((Array.append neededQuantColsIdentifiers neededQuantColsWithValues) |> Set.ofSeq)
                    |> Set.toArray
                // drops all in 'quantColsNotNeeded' saved columns from the quant table
                let filteredQuantTable =
                    quantTable
                    |> dropCols quantColsNotNeeded
                let quantTableWithRatios =
                    // calculate 14N/15N ratios based on the corresponding fields from the quant table
                    let ratios =
                        let n14 = filteredQuantTable.GetColumn<float>"N14Quant"
                        let n15 = filteredQuantTable.GetColumn<float>"N15Quant"
                        n14/n15

                    filteredQuantTable
                    // add the calculated ratios to the quant table
                    |> Frame.addCol "Ratio" ratios
                    // group all rows based on the peptide sequence. This groups different charges and global modifications of the same peptide
                    |> Frame.groupRowsUsing (fun k ser -> (ser.GetAs<string>("StringSequence")))
                    // makes sure only columns with values are present (i.e. charge/global mod are removed)
                    |> Frame.filterCols (fun ck _ -> 
                        (
                            (
                                neededQuantColsWithValues
                                // appends "Ratio" column since it's not included in the original columns with values
                                |> Array.append [|"Ratio"|]
                            )
                            |> Set.ofArray
                        ).Contains(ck)
                    )
                    // aggregates the columns over the peptide sequence with a defined method (i.e. average,median,...)
                    |> Frame.applyLevel (fun (sequence,index) -> sequence) Stats.mean
                    |> Frame.dropCol "PEPValue"
                // set of every peptide present in the quant table based on the row keys
                let peptidesPresent = quantTableWithRatios.RowKeys |> Set.ofSeq
                // table containing every peptide that is mapped to a protein group. The row keys are the protein groups and the PeptideSequences column contains a string array with all peptides mapping to them.
                let peptidesMapped =
                    protTable
                    |> Frame.getCol "PeptideSequences"
                    |> Series.mapValues (fun (s : string)-> s.Split(';'))
                // calculates count of peptides mapped to each protein group
                let distinctPeptideCount =
                    protTable
                    |> Frame.getCol "PeptideSequences"
                    |> Series.map (fun _ x ->
                        x
                        |> String.split ';'
                        |> Array.length
                    )
                // an "empty" table based on the prot table. It only contains the row keys with the protein groups.
                let baseTable =
                    protTable
                    |> dropCols [|"EvidenceClass"; "PeptideSequences"|]
                    |> Frame.addCol "DistinctPeptideCount" distinctPeptideCount
                // the columns in the array are selected in the filtered quant table and compared to the protein groups. If the peptide is mapping to a protein group, then
                // an entry is added to the base table at the corresponding protein group.
                [|"Ratio";"N14Quant";"N15Quant"|]
                |> Array.map (fun name ->
                    baseTable.AddColumn (name, getAggregatedPeptidesVals peptidesPresent peptidesMapped quantTableWithRatios name)
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