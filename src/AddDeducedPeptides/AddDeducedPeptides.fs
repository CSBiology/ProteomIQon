namespace ProteomIQon

open System.IO
open ProteomIQon.Dto
open FSharpAux
open FSharpAux.IO

module AddDeducedPeptides = 
  
    let addDeducedPeptides (protFiles: string []) (quantFiles: string []) (outFolder: string) = 
        let prot =
            protFiles
            |> Array.map (fun path ->
                SeqIO.Seq.fromFileWithCsvSchema<ProteinInferenceResult>(path,'\t',false, skipLines=1)
                |> Array.ofSeq
            )
        let quant =
            quantFiles
            |> Array.map (fun path ->
                SeqIO.Seq.fromFileWithCsvSchema<QuantificationResult>(path,'\t',false, skipLines=1)
                |> Array.ofSeq
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
        let pepProtMap =
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
            |> Array.distinctBy fst
            |> Map.ofArray
        let assignProts =
            quant
            |> Array.map (fun quantFile ->
                let presentPeptides =
                    quantFile
                    |> Array.map (fun x -> x.StringSequence)
                    |> Array.distinct
                let protInfResultExpanded =
                    presentPeptides
                    |> Array.choose (fun pep ->
                        let foundProtInf = pepProtMap.TryFind pep
                        match foundProtInf with
                        | Some protInf -> Some (pep,protInf)
                        | None -> None
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
            )
        assignProts
