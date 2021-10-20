namespace ProteomIQon

open System.IO
open ProteomIQon.Dto
open FSharpAux
open FSharpAux.IO
open FSharpAux.IO.SchemaReader

module AddDeducedPeptides = 
  
    let addDeducedPeptides (protFiles: string []) (quantFiles: string []) (outDirectory: string) = 
        let prot =
            protFiles
            |> Array.map (fun path ->
                Csv.CsvReader<ProteinInferenceResult>(SchemaMode=Csv.Fill).ReadFile(path,'\t',false,1)
                |> Array.ofSeq
            )
        let quant =
            quantFiles
            |> Array.map (fun path ->
                Csv.CsvReader<QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(path,'\t',false,1)
                |> Array.ofSeq, path
            )
        let ensureEqualProtDetermination =
            prot
            |> Array.concat
            |> Array.groupBy (fun x -> x.ProteinGroup)
            |> Array.map snd
            |> Array.map (fun prots -> 
                prots
                |> Array.distinctBy (fun x -> x.DecoyScore, x.TargetScore, x.QValue)
            )
            |> Array.exists (fun x -> x.Length <> 1)
        if ensureEqualProtDetermination then
            failwith "If you are running this tool please infer proteins over all files"
        let pepProt =
            prot
            |> Array.collect (fun protFile ->
                protFile
                |> Array.collect (fun prot ->
                    prot.PeptideSequence
                    |> String.split ';'
                    |> Array.map (fun pep ->
                        pep, prot
                    )
                )
            )
            |> Array.distinctBy (fun (pep,prot) -> pep, prot.ProteinGroup)
        quant
        |> Array.map (fun (quantFile, filePath) ->
            let outFilePath =
                let filename = Path.GetFileNameWithoutExtension filePath
                Path.Combine[|outDirectory; $"{filename}.prot"|]
            let presentPeptides =
                quantFile
                |> Array.map (fun x -> x.StringSequence)
                |> Set.ofArray
            let protInfResultExpanded =
                pepProt
                |> Array.choose (fun (pep, prot) ->
                    if presentPeptides.Contains pep then
                        Some (pep,prot)
                    else
                        None
                )
            let protInfResultGrouped =
                protInfResultExpanded
                |> Array.groupBy (fun (pep,protInf) -> protInf.ProteinGroup)
                |> Array.map snd
            let protInfResultsCombined =
                protInfResultGrouped
                |> Array.map (fun group ->
                    let peps =
                        group
                        |> Array.map fst
                        |> String.concat ";"
                    group
                    |> Array.head
                    |> snd
                    |> fun protInf ->
                        {protInf with PeptideSequence = peps}

                )
            protInfResultsCombined
            |> SeqIO.Seq.CSV "\t" true false
            |> SeqIO.Seq.writeOrAppend outFilePath
        )
