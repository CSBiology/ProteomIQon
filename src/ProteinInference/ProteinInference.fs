namespace ProteomIQon

open BioFSharp
open System.Data.SQLite
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open BioFSharp.Mz.ProteinInference
open BioFSharp.IO
open FSharpAux
open PeptideClassification
open ProteomIQon
open FSharpAux.IO
open FSharpAux.IO.SchemaReader.Csv
open Domain
open Plotly.NET
open FSharp.Stats
open BioFSharp.IO.GFF3
module ProteinInference =

    /// Represents one peptide-entry with its evidence class and the proteins it points to.
    type ClassInfo =
        {
            Sequence : string
            Class    : PeptideClassification.PeptideEvidenceClass
            Proteins : string []
        }

    /// Calculates the fdr with the chosen fdr method
    let initFDR (fdrMethod: FDRMethod) (data: ProteinInference'.InferredProteinClassItemScored[]) (proteinDB: (string*string)[]) =
        match fdrMethod with
        |Conservative     -> 1.
        |MAYU             -> FDRControl'.calculateFDRwithMAYU data proteinDB
        |DecoyTargetRatio -> FDRControl'.calculateFDRwithDecoyTargetRatio data

    /// Calculates the q value with the chosen method
    let initQValue (qValueMethod: QValueMethod) bandwidth (data: ProteinInference'.InferredProteinClassItemScored[]) (proteinDB: (string*string)[]) =
        match qValueMethod with
        |Storey-> FDRControl'.calculateQValueStorey data
        |LogisticRegression fdrMethod ->
            let fdr = initFDR fdrMethod data proteinDB
            FDRControl'.calculateQValueLogReg fdr bandwidth data 

    /// Given a ggf3 and a fasta file, creates a collection of all theoretically possible peptides and the proteins they might
    /// originate from
    ///
    /// ProteinClassItem: For a single peptide Sequence, contains information about all proteins it might originate from and
    /// its evidence class.
    ///
    /// No experimental data
    let createClassItemCollection gff3Path (memoryDB: SQLiteConnection) regexPattern rawFolderPath=

        let logger = Logging.createLogger "ProteinInference_createClassItemCollection"

        logger.Trace (sprintf "Regex pattern: %s" regexPattern)

        let rawFilePaths = System.IO.Directory.GetFiles (rawFolderPath, "*.qpsm")
                           |> List.ofArray

        // retrieves peptide sequences and IDs from input files
        let psmInputs =
            rawFilePaths
            |> List.map (fun filePath ->
                Seq.fromFileWithCsvSchema<ProteomIQon.ProteinInference'.PSMInput>(filePath, '\t', true,schemaMode = SchemaModes.Fill)
                |> Seq.toList
                )

         //gets distinct IDs of all peptides represented in the raw files
        let inputPepSeqIDs =
            psmInputs
            |> List.collect (fun psmArray ->
                psmArray
                |> List.map (fun psm -> psm.PepSequenceID)
                )
            |> List.distinct

        //list of proteins tupled with list of possible peptides found in psm
        let accessionSequencePairs =
            let preparedProtPepFunc = ProteomIQon.SearchDB'.getProteinPeptideLookUpFromFileBy memoryDB
            inputPepSeqIDs
            |> List.map (fun pepID -> preparedProtPepFunc pepID)
            |> List.concat
            |> List.groupBy (fun (protein, _)-> protein)
            |> List.map (fun (protein, pepList) -> protein, (pepList |> List.map (fun (_,pep)-> BioArray.ofAminoAcidString pep)))

        /// Create proteinModelInfos: Group all genevariants (RNA) for the gene loci (gene)
        ///
        /// The proteinModelInfos
        logger.Trace "Reading GFF3 file"
        let proteinModelInfos =
            match gff3Path with 
            | Some gff3Path -> 
                try
                    GFF3.fromFileWithoutFasta gff3Path
                    |> ProteomIQon.ProteinInference'.assignTranscriptsToGenes regexPattern
                with
                | err ->
                    printfn "ERROR: Could not read gff3 file %s" gff3Path
                    failwithf "%s" err.Message
            | None ->                 
                ProteomIQon.SearchDB'.selectProteins memoryDB
                |> List.mapi (fun i (id,_) -> 
                    [
                    GFFEntryLine (createGFFEntry "Unknown" "." "gene" "." "." "." "." "." (sprintf "ID=%i" i) "");
                    GFFEntryLine (createGFFEntry "Unknown" "." "mRNA" "." "." "." "." "." (sprintf "ID=%s;Parent=%i" id i) "")
                    ]
                )
                |> Seq.concat
                |> ProteomIQon.ProteinInference'.assignTranscriptsToGenes regexPattern

        //reads from file to an array of FastaItems.
        logger.Trace "Assigning FastA sequences to protein model info"
        /// Assigned fasta sequences to model Infos
        let proteinModels =
            try
                accessionSequencePairs
                |> Seq.mapi (fun i (header,sequence) ->
                    let regMatch = System.Text.RegularExpressions.Regex.Match(header,regexPattern)
                    if regMatch.Success then
                        match Map.tryFind regMatch.Value proteinModelInfos with
                        | Some modelInfo ->
                            Some (createProteinModel modelInfo sequence)
                        | None ->
                            failwithf "Could not find protein with id %s in gff3 File" regMatch.Value
                    else
                        logger.Trace (sprintf "Could not extract id of header \"%s\" with regexPattern %s" header regexPattern)
                        createProteinModelInfo header "Placeholder" StrandDirection.Forward (sprintf "Placeholder%i" i) 0 Seq.empty Seq.empty
                        |> fun x -> createProteinModel x sequence
                        |> Some
                    )
            with
            | err ->
                printfn "Could not assign FastA sequences to RNAs"
                failwithf "%s" err.Message
        logger.Trace "Build classification map"
        try
            let ppRelationModel = ProteomIQon.ProteinInference'.createPeptideProteinRelation proteinModels

            let spliceVariantCount = PeptideClassification.createLocusSpliceVariantCount ppRelationModel

            let classified =
                ppRelationModel.GetArrayOfKeys
                |> Array.map (fun peptide ->
                    let proteinInfo = (ppRelationModel.TryGetByKey peptide).Value
                    let proteinIds = Seq.map (fun (x : PeptideClassification.ProteinModelInfo<string,string,string>) -> x.Id) proteinInfo
                    let (c,x) = PeptideClassification.classify spliceVariantCount (peptide, (ppRelationModel.TryGetByKey peptide).Value)
                    {
                        Sequence = BioArray.toString x;
                        Class = c;
                        Proteins = (proteinIds |> Seq.toArray)
                    }
                    )

            classified |> Array.map (fun ci -> ci.Sequence,(createProteinClassItem ci.Proteins ci.Class ci.Sequence)) |> Map.ofArray, psmInputs
        with
        | err ->
            printfn "\nERROR: Could not build classification map"
            failwithf "\t%s" err.Message

    let readAndInferFile classItemCollection protein peptide groupFiles outDirectory rawFolderPath psmInputs (dbConnection: SQLiteConnection) (qValMethod: Domain.QValueMethod) =

        let logger = Logging.createLogger "ProteinInference_readAndInferFile"

        let rawFilePaths = 
            System.IO.Directory.GetFiles (rawFolderPath, "*.qpsm")
            |> Array.toList

        let outFiles: string list =
            rawFilePaths
            |> List.map (fun filePath ->
                //let foldername = (rawFolderPath.Split ([|"\\"|], System.StringSplitOptions.None))
                outDirectory + @"\" (*+ foldername.[foldername.Length - 1]*) + "\\" + (System.IO.Path.GetFileNameWithoutExtension filePath) + ".prot"
                )

        let dbParams = ProteomIQon.SearchDB'.getSDBParams dbConnection

        logger.Trace "Getting proteins with peptide sequence from db"

        let proteinsDB =
            ProteomIQon.SearchDB'.selectProteins dbConnection
            |> Array.ofList
            |> Array.map (fun (protein, peptideSequence) ->
                //is this the best way to get the protein names?
                (protein |> String.split ' ').[0], peptideSequence
            )

        logger.Trace "Perfom digestion on reversed proteins for decoy set"

        // Array of prtoein Accessions tupled with their reverse digested peptides
        let reverseProteins =
            proteinsDB
            |> Array.map (fun (name, sequence) ->
                name,
                Digestion.BioArray.digest dbParams.Protease 0 ((sequence |> String.rev) |> BioArray.ofAminoAcidString)
                |> Digestion.BioArray.concernMissCleavages dbParams.MinMissedCleavages dbParams.MaxMissedCleavages
                |> fun x -> x |> Array.map (fun y -> y.PepSequence |> List.toArray |> BioArray.toString)
                )

        logger.Trace "Map peptide sequences to proteins"

        let classifiedProteins =
            List.map2 (fun psmInput (outFile: string) ->
                try
                    psmInput
                    |> List.map (fun (s: ProteomIQon.ProteinInference'.PSMInput) ->
                        let s = s.Seq
                        match Map.tryFind (ProteinInference'.removeModification s) classItemCollection with
                        | Some (x:ProteinClassItem<'sequence>) -> createProteinClassItem x.GroupOfProteinIDs x.Class s
                        | None ->
                            match Map.tryFind (ProteinInference'.removeModification ("M" + s)) classItemCollection with
                            | Some (x:ProteinClassItem<'sequence>) -> createProteinClassItem x.GroupOfProteinIDs x.Class s
                            | None ->
                                failwithf "Could not find sequence %s in classItemCollection" s
                        ),
                        outFile
                with
                | err ->
                    printfn "Could not map sequences of file %s to proteins:" (System.IO.Path.GetFileNameWithoutExtension outFile)
                    failwithf "%s" err.Message
                ) psmInputs outFiles

        if groupFiles then

            logger.Trace "Create combined list"

            let combinedClasses =
                List.collect (fun (pepSeq,_) -> pepSeq) classifiedProteins
                |> BioFSharp.Mz.ProteinInference.inferSequences protein peptide

            logger.Trace "Create a map with peptide scores"
            // Peptide score Map
            let peptideScoreMap = ProteomIQon.ProteinInference'.createPeptideScoreMap psmInputs

            logger.Trace "Assign scores to decoy set"

            // Assigns scores to reverse digested Proteins using the peptideScoreMap
            let reverseProteinScores = ProteomIQon.ProteinInference'.createReverseProteinScores reverseProteins peptideScoreMap

            logger.Trace "Assign target and decoy scores inferred proteins"

            // Scores each inferred protein and assigns each protein where a reverted peptide was hit its score
            let combinedScoredClasses =
                combinedClasses
                |> Array.ofSeq
                |> Array.map (fun inferredPCI ->
                    // Looks up all peptides assigned to the protein and sums up their score
                    let peptideScore = (ProteomIQon.ProteinInference'.assignPeptideScores inferredPCI.PeptideSequence peptideScoreMap)
                    // Looks if the reverse protein has been randomly matched and assigns the score
                    let decoyScore = (ProteomIQon.ProteinInference'.assignDecoyScoreToTargetScore inferredPCI.GroupOfProteinIDs reverseProteinScores)

                    ProteomIQon.ProteinInference'.createInferredProteinClassItemScored
                        inferredPCI.GroupOfProteinIDs inferredPCI.Class inferredPCI.PeptideSequence
                        peptideScore
                        decoyScore
                        false
                        // for the same score target "wins"
                        (decoyScore > peptideScore)
                        true
                    )

            // creates InferredProteinClassItemScored type for decoy proteins that have no match
            let reverseNoMatch =
                let proteinsPresent =
                    combinedScoredClasses
                    |> Array.collect (fun prots -> prots.GroupOfProteinIDs |> String.split ';')
                reverseProteinScores
                |> Map.toArray
                |> Array.choose (fun (protein, (score, peptides)) ->
                    if (proteinsPresent |> Array.contains protein) then
                        None
                    else
                        // peptides are the peptides which point to the reverse digested protein. This info is currently unused, since in cases where a partner was found this field still contains
                        // the peptides that point to the forward digested protein.
                        Some (ProteomIQon.ProteinInference'.createInferredProteinClassItemScored protein PeptideEvidenceClass.Unknown peptides 0. score true true true)
                )

            logger.Trace "Calculate q values for all proteins"

            let combWithReverse = Array.append combinedScoredClasses reverseNoMatch

            let bandwidth =
                let data = 
                    combWithReverse
                    |> Array.map (fun x ->
                        if x.DecoyBigger then
                            x.DecoyScore
                        else
                            x.TargetScore
                    )
                FSharp.Stats.Distributions.Bandwidth.nrd0 data

            // Assign q values to each protein
            let combinedScoredClassesQVal =
                let decoyBiggerF = (fun (item: ProteinInference'.InferredProteinClassItemScored) -> item.DecoyBigger)
                let targetScoreF = (fun (item: ProteinInference'.InferredProteinClassItemScored) -> item.TargetScore)
                let decoyScoreF  = (fun (item: ProteinInference'.InferredProteinClassItemScored) -> item.DecoyScore)
                let qValueFunction = initQValue qValMethod bandwidth combWithReverse proteinsDB decoyBiggerF decoyScoreF targetScoreF
                let qValuesAssigned =
                    combWithReverse
                    |> Array.map (FDRControl'.assignQValueToIPCIS qValueFunction)
                qValuesAssigned

            ProteinInference'.qValueHitsVisualization bandwidth combinedScoredClassesQVal outDirectory groupFiles

            // Assign results to files in which they can be found
            classifiedProteins
            |> List.iter (fun (prots,outFile) ->
                let pepSeqSet = prots |> List.map (fun x -> x.PeptideSequence) |> Set.ofList
                logger.Trace (sprintf "Writing now %s"(System.IO.Path.GetFileNameWithoutExtension outFile))
                let combinedInferenceresult =
                    combinedScoredClassesQVal |> Array.filter (fun inferredPCIS -> not inferredPCIS.Decoy)
                    |> Seq.choose (fun ic ->
                        let filteredPepSet =
                            ic.PeptideSequence
                            |> Array.filter (fun pep -> Set.contains pep pepSeqSet)
                        if filteredPepSet = [||] then
                            None
                        else
                            Some (ProteomIQon.ProteinInference'.createInferredProteinClassItemOut
                                ic.GroupOfProteinIDs ic.Class (ProteinInference'.proteinGroupToString filteredPepSet) ic.TargetScore ic.DecoyScore ic.QValue)
                        )

                combinedInferenceresult
                |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
                |> Seq.write outFile
                logger.Trace (sprintf "File written to %s" outFile)
            )
        else
            classifiedProteins
            |> List.iter2 (fun psm (sequences,outFile) ->
                logger.Trace (sprintf "start inferring %s" (System.IO.Path.GetFileNameWithoutExtension outFile))
                let inferenceResult = BioFSharp.Mz.ProteinInference.inferSequences protein peptide sequences

                logger.Trace "Create a map with peptide scores"

                // Peptide Score Map
                let peptideScoreMap = ProteomIQon.ProteinInference'.createPeptideScoreMap [psm]

                logger.Trace "Assign scores to decoy set"

                // Assigns scores to reverse digested Proteins using the peptideScoreMap
                let reverseProteinScores = ProteomIQon.ProteinInference'.createReverseProteinScores reverseProteins peptideScoreMap

                logger.Trace "Assign target and decoy scores inferred proteins"

                // Scores each inferred protein and assigns each protein where a reverted peptide was hit its score
                let inferenceResultScored =
                    inferenceResult
                    |> Array.ofSeq
                    |> Array.map (fun inferredPCI ->
                        // Looks up all peptides assigned to the protein and sums up their score
                        let peptideScore = (ProteomIQon.ProteinInference'.assignPeptideScores inferredPCI.PeptideSequence peptideScoreMap)
                        // Looks if the reverse protein has been randomly matched and assigns the score
                        let decoyScore = (ProteomIQon.ProteinInference'.assignDecoyScoreToTargetScore inferredPCI.GroupOfProteinIDs reverseProteinScores)

                        ProteomIQon.ProteinInference'.createInferredProteinClassItemScored
                            inferredPCI.GroupOfProteinIDs inferredPCI.Class inferredPCI.PeptideSequence
                            peptideScore
                            decoyScore
                            false
                            (decoyScore > peptideScore)
                            true
                        )

                // creates InferredProteinClassItemScored type for decoy proteins that have no match
                let reverseNoMatch =
                    let proteinsPresent =
                        inferenceResultScored
                        |> Array.collect (fun prots -> prots.GroupOfProteinIDs |> String.split ';')
                    reverseProteinScores
                    |> Map.toArray
                    |> Array.choose (fun (protein, (score, peptides)) ->
                        if (proteinsPresent |> Array.contains protein) then
                            None
                        else
                            Some (ProteomIQon.ProteinInference'.createInferredProteinClassItemScored protein PeptideEvidenceClass.Unknown peptides 0. score true true true)
                    )

                logger.Trace "Calculate q values for all proteins"

                let combWithReverse = Array.append inferenceResultScored reverseNoMatch

                let bandwidth =
                    let data = 
                        combWithReverse
                        |> Array.map (fun x ->
                            if x.DecoyBigger then
                                x.DecoyScore
                            else
                                x.TargetScore
                        )
                    FSharp.Stats.Distributions.Bandwidth.nrd0 data

                // Assign q values to each protein
                let inferenceResultScoredQVal =
                    let decoyBiggerF = (fun (item: ProteinInference'.InferredProteinClassItemScored) -> item.DecoyBigger)
                    let targetScoreF = (fun (item: ProteinInference'.InferredProteinClassItemScored) -> item.TargetScore)
                    let decoyScoreF  = (fun (item: ProteinInference'.InferredProteinClassItemScored) -> item.DecoyScore)
                    let qValueFunction = initQValue qValMethod bandwidth combWithReverse proteinsDB decoyBiggerF decoyScoreF targetScoreF
                    let qValuesAssigned =
                        combWithReverse
                        |> Array.map (FDRControl'.assignQValueToIPCIS qValueFunction)
                    qValuesAssigned

                ProteinInference'.qValueHitsVisualization bandwidth inferenceResultScoredQVal outFile groupFiles

                inferenceResultScoredQVal
                |> Array.filter (fun inferredPCIQ -> not inferredPCIQ.Decoy)
                |> Array.map (fun ic ->
                    ProteomIQon.ProteinInference'.createInferredProteinClassItemOut 
                        ic.GroupOfProteinIDs ic.Class (ProteinInference'.proteinGroupToString ic.PeptideSequence) ic.TargetScore ic.DecoyScore ic.QValue
                )
                |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
                |> Seq.write outFile
            ) psmInputs

    let inferProteins gff3Location dbConnection (proteinInferenceParams: ProteinInferenceParams) outDirectory rawFolderPath =

        let logger = Logging.createLogger "ProteinInference_inferProteins"

        logger.Trace (sprintf "InputFilePath = %s" rawFolderPath)
        logger.Trace (sprintf "InputGFF3Path = %A" gff3Location)
        logger.Trace (sprintf "OutputFilePath = %s" outDirectory)
        logger.Trace (sprintf "Protein inference parameters = %A" proteinInferenceParams)

        logger.Trace "Copy peptide DB into Memory."
        let memoryDB = SearchDB.copyDBIntoMemory dbConnection
        logger.Trace "Copy peptide DB into Memory: finished."

        logger.Trace "Start building ClassItemCollection"
        let classItemCollection, psmInputs = createClassItemCollection gff3Location memoryDB proteinInferenceParams.ProteinIdentifierRegex rawFolderPath
        logger.Trace "Classify and Infer Proteins"
        readAndInferFile classItemCollection proteinInferenceParams.Protein proteinInferenceParams.Peptide
                         proteinInferenceParams.GroupFiles outDirectory rawFolderPath psmInputs dbConnection proteinInferenceParams.GetQValue