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

    let sortTables (quantFiles: string[]) (protFiles: string[]) =

        let quantColumnTypes =
            createSchema "float" neededQuantColsWithValues
        let tables =
            Array.map2 ( fun quantFile protFile ->
                let quantTable: Frame<int,string> =
                    Frame.ReadCsv(path=quantFile ,separators="\t", schema=quantColumnTypes)
                let protTable : Frame<string,string> =
                    Frame.ReadCsv(path=protFile ,separators="\t")
                    |> Frame.indexRowsString "ProteinID"
                let quantColsNotNeeded =
                    Set.difference (quantTable.ColumnKeys |> Set.ofSeq) ((Array.append neededQuantColsIdentifiers neededQuantColsWithValues) |> Set.ofSeq)
                    |> Set.toArray
                let filteredQuantTable =
                    quantTable
                    |> dropCols quantColsNotNeeded
                let quantTableWithRatios =
                    let ratios =
                        let n14 = filteredQuantTable.GetColumn<float>"N14Quant"
                        let n15 = filteredQuantTable.GetColumn<float>"N15Quant"
                        n14/n15
                    filteredQuantTable
                    |> Frame.addCol "Ratio" ratios
                    |> Frame.groupRowsUsing (fun k ser -> (ser.GetAs<string>("StringSequence")))
                    // makes sure only columns with values are present
                    |> Frame.filterCols (fun ck _ -> ((neededQuantColsWithValues |> Array.append [|"Ratio"|])|> Set.ofArray).Contains(ck) )
                    |> Frame.applyLevel (fun (sequence,index) -> sequence) Stats.mean
                    |> Frame.dropCol "PEPValue"
                let peptidesPresent = quantTableWithRatios.RowKeys |> Set.ofSeq
                let peptidesMapped =
                    protTable
                    |> Frame.getCol "PeptideSequences"
                    |> Series.mapValues (fun (s : string)-> s.Split(';'))
                let distinctPeptideCount =
                    protTable
                    |> Frame.getCol "PeptideSequences"
                    |> Series.map (fun _ x -> 
                        x
                        |> String.split ';'
                        |> Array.length
                    )
                let baseTable =
                    protTable
                    |> dropCols [|"EvidenceClass"; "PeptideSequences"|]
                    |> Frame.addCol "DistinctPeptideCount" distinctPeptideCount
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
