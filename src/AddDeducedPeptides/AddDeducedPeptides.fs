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
        // Fail when the protein inference was not performed on combined files
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
        // Expand all inferred proteins by their peptides and distinct by peptide/protein pairs
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
        // Map over quant files and keep only peptide/protein pairs of peptides present in the file
        quant
        |> Array.map (fun (quantFile, filePath) ->
            let outFilePath =
                let filename = Path.GetFileNameWithoutExtension filePath
                Path.Combine[|outDirectory; $"{filename}.prot"|]
            // Obtain sequences of present peptides and filter for proteingroups containing those peptides
            let presentPeptides =
                quantFile
                |> Array.map (fun x -> x.StringSequence)
                |> Set.ofArray
            let protInfResultExpanded =
                pepProt
                |> Array.filter (fun (pep, prot) ->
                    presentPeptides.Contains pep
                )
            // Group by Protein Group to combine them
            let protInfResultGrouped =
                protInfResultExpanded
                |> Array.groupBy (fun (pep,protInf) -> protInf.ProteinGroup)
                |> Array.map snd
            // Combine the Protein Groups by creating a new PeptideSequence field based on the peptides present in the file for this group
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
