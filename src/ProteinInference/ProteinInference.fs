namespace ProteomIQon

open BioFSharp
open System.Data
open System.Data.SQLite
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open BioFSharp.Mz.ProteinInference
open BioFSharp.IO
open FSharpAux
open System.Text.RegularExpressions
open PeptideClassification
open ProteomIQon
open FSharpAux.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Csv
open FSharpAux.IO.SchemaReader.Attribute
open GFF3
open Domain
open FSharp.Plotly
open FSharp.Stats

module ProteinInference =

    /// This type represents one element of the final output. It's the protein, with all the peptides that were measured and are pointing to it.
    type IntermediateResult =
        {
            ProteinID        : string
            EvidenceClass    : PeptideEvidenceClass
            PeptideSequences : string
            SumOfScores      : float
            DecoyScore       : float
        }

    /// Packages info into one element of the intermediate output. It's the protein, with all the peptides that were measured and are pointing to it.
        static member createIntermediateResult protID evidenceClass sequences score decoyScore =
            {
                ProteinID        = protID
                EvidenceClass    = evidenceClass
                PeptideSequences = sequences
                SumOfScores      = score
                DecoyScore       = decoyScore
            }

    /// Represents one peptide-entry with its evidence class and the proteins it points to.
    type ClassInfo =
        {
            Sequence : string
            Class    : PeptideClassification.PeptideEvidenceClass
            Proteins : string []
        }

    type PSMInput =
        {
            [<FieldAttribute("PepSequenceID")>]
            PepSequenceID   : int
            [<FieldAttribute("StringSequence")>]
            Seq             :string
            [<FieldAttribute("PercolatorScore")>]
            PercolatorScore : float
        }

    // Input for QValue calulation
    type QValueInput =
        {
            Score    : float
            IsDecoy  : bool
        }

    let createQValueInput score isDecoy =
        {
            Score     = score
            IsDecoy   = isDecoy
        }

    /// For a group of proteins, contains information about all peptides that might be used for its quantification and score / q-value calculated for it.
    type InferredProteinClassItemScored<'sequence> =
        {
            GroupOfProteinIDs: string
            PeptideSequence  : 'sequence []
            Class            : PeptideEvidenceClass
            TargetScore      : float
            DecoyScore       : float
            QValue           : float
            Decoy            : bool
            DecoyBigger      : bool
        }

    let createInferredProteinClassItemScored proteinIDs evidenceClass peptideSequences targetScore decoyScore qValue isDecoy decoyBigger =
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence   = peptideSequences
            Class             = evidenceClass
            TargetScore       = targetScore
            DecoyScore        = decoyScore
            QValue            = qValue
            Decoy             = isDecoy
            DecoyBigger       = decoyBigger
        }

    ///checks if GFF line describes gene
    let isGene (item: GFFLine<seq<char>>) =
        match item with
        | GFFEntryLine x -> x.Feature = "gene"
        | _ -> false

    ///checks if GFF line describes rna
    let isRNA (item: GFFLine<seq<char>>) =
        match item with
        | GFFEntryLine x -> if x.Feature = "mRNA" then Some x else None
        | _ -> None

    /// Reads geographical information about protein from gff entry and builds the modelinfo of it
    /// This function takes an RNA gff3 entry and therefore will contain the splice variant id of the gene in its result.
    /// This splice variant id should be the same in the given FastA-file.
    let createProteinModelInfoFromEntry i locus (entry:GFFEntry) =
        let attributes = entry.Attributes
        /// Same as in FastA file
        let spliceVariantID =
            match Map.tryFind "Name" attributes with
            | Some res ->
                res.Head
            | None ->
                failwithf "could not find spliceVariantId for locus %s" locus
        let chromosomeID = entry.Seqid

        let direction =
            match entry.Strand with
            |'+' -> PeptideClassification.StrandDirection.Forward
            |'-' -> PeptideClassification.StrandDirection.Reverse

        PeptideClassification.createProteinModelInfo spliceVariantID chromosomeID direction locus i Seq.empty Seq.empty

    /// By reading GFF creates the protein models (relationships of proteins to each other) which basically means grouping the rnas over the gene loci
    /// TODO: Don't group over order but rather group over id
    let assignTranscriptsToGenes regexPattern (gffLines: seq<GFF3.GFFLine<seq<char>>>)  =
        gffLines
        // transcripts are grouped by the gene they originate from
        |> Seq.groupWhen isGene
        |> Seq.map (fun group ->
            match Seq.head group with
            | GFFEntryLine x ->
                let locus = x.Attributes.["ID"].Head // =genename, this value is used to assign mRNAs of the same gene together
                group
                |> Seq.choose isRNA //the transcripts of the gene are chosen
                |> Seq.mapi (fun i element ->
                    // every transcript of gene gets its own number i and other info is collected from element and used for info of protein
                    let modelInfo = createProteinModelInfoFromEntry i locus element
                    let r = System.Text.RegularExpressions.Regex.Match(modelInfo.Id,regexPattern)
                    // the gff3 id has to be matched with the sequence in the fasta file. therefore the regexpattern is used
                    if r.Success then
                        r.Value,
                        modelInfo
                    else
                        failwithf "could not match gff3 entry id %s with regexpattern %s. Either gff3 file is corrupt or regexpattern is not correct" modelInfo.Id regexPattern
                    )

            | _ -> Seq.empty
            )
        |> Seq.concat
        |> Map.ofSeq

    /// Returns SearchDbParams of a existing database by filePath
    let getSDBParamsBy (cn :SQLiteConnection)=
        let cn =
            match cn.State with
            | ConnectionState.Open ->
                cn
            | ConnectionState.Closed ->
                cn.Open()
                cn
            | _ as x -> failwith "Data base is busy."
        match Db.SQLiteQuery.selectSearchDbParams cn with
        | Some (iD,name,fo,fp,pr,minmscl,maxmscl,mass,minpL,maxpL,isoL,mMode,fMods,vMods,vThr) ->
            createSearchDbParams
                name fo fp id (Digestion.Table.getProteaseBy pr) minmscl maxmscl mass minpL maxpL
                    (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchInfoIsotopic list>(isoL)) (Newtonsoft.Json.JsonConvert.DeserializeObject<MassMode>(mMode)) (massFBy (Newtonsoft.Json.JsonConvert.DeserializeObject<MassMode>(mMode)))
                        (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchModification list>(fMods)) (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchModification list>(vMods)) vThr
        | None ->
            failwith "This database does not contain any SearchParameters. It is not recommended to work with this file."

    /// Prepares statement to select a Protein Accession entry by ID
    let prepareSelectProteinAccessionByID (cn:SQLiteConnection) (tr) =
        let querystring = "SELECT Accession FROM Protein WHERE ID=@id "
        let cmd = new SQLiteCommand(querystring, cn, tr)
        cmd.Parameters.Add("@id", DbType.Int32) |> ignore
        (fun (id:int32)  ->
            cmd.Parameters.["@id"].Value <- id
            use reader = cmd.ExecuteReader()
            match reader.Read() with
            | true  -> (reader.GetString(0))
            | false -> ""
        )

    /// Prepares statement to select a Peptide Sequence entry by ID
    let prepareSelectPepSequenceByPepSequenceID (cn:SQLiteConnection) (tr) =
        let querystring = "SELECT Sequence FROM PepSequence WHERE ID=@pepSequenceID"
        let cmd = new SQLiteCommand(querystring, cn, tr)
        cmd.Parameters.Add("@pepSequenceID", DbType.Int32) |> ignore
        (fun (pepSequenceID:int)  ->
            cmd.Parameters.["@pepSequenceID"].Value <- pepSequenceID
            use reader = cmd.ExecuteReader()
            reader.Read() |> ignore
            reader.GetString(0)
            )

    /// Prepares a function which returns a list of protein Accessions tupled with the peptide sequence whose ID they were retrieved by
    let getProteinPeptideLookUpFromFileBy (memoryDB: SQLiteConnection) =
        let tr = memoryDB.BeginTransaction()
        let selectCleavageIdxByPepSeqID   = Db.SQLiteQuery.prepareSelectCleavageIndexByPepSequenceID memoryDB tr
        let selectProteinByProtID         = prepareSelectProteinAccessionByID memoryDB tr
        let selectPeptideByPepSeqID       = prepareSelectPepSequenceByPepSequenceID memoryDB tr
        (fun pepSequenceID ->
                selectCleavageIdxByPepSeqID pepSequenceID
                |> List.map (fun (_,protID,pepID,_,_,_) -> selectProteinByProtID protID, selectPeptideByPepSeqID pepID )
        )

    /// Returns Accession and Sequence of Proteins from SearchDB
    let selectProteins (cn:SQLiteConnection) =
        let selectProteins =
            let querystring = "SELECT Accession, Sequence FROM Protein"
            let cmd = new SQLiteCommand(querystring, cn)
            use reader = cmd.ExecuteReader()
            (let rec loop (list: (string*string) list) =
                match reader.Read() with
                | false -> list |> List.rev
                | true  -> loop ((reader.GetString(0), reader.GetString(1))::list)
            loop []
            )
        selectProteins

    /// Creates a lookup data base to assign peptides to the proteins they are contained in
    let createPeptideProteinRelation (protModels:seq<ProteinModel<'id,'chromosomeId,'geneLocus,'sequence list> option>) =
        let ppRelation = BidirectionalDictionary<'sequence,ProteinModelInfo<'id,'chromosomeId,'geneLocus>>()
        protModels
        |> Seq.iter (fun prot ->
                        // insert peptide-protein relationship
                        // Todo: change type of proteinID in digest
                        match prot with
                        | Some proteinModel ->
                            proteinModel.Sequence
                            |> Seq.iter (fun pepSequence -> ppRelation.Add pepSequence proteinModel.ProteinModelInfo)
                        | None                   -> ()
                    )
        ppRelation

    // Creates a Map of peptides with their highest found score
    let createPeptideScoreMap (psmInputs: PSMInput list list) =
        psmInputs
        |> List.concat
        |> List.groupBy (fun psm -> psm.Seq)
        |> List.map (fun (sequence, psmList) ->
            sequence,
            psmList
            |> List.sortByDescending (fun psm -> psm.PercolatorScore)
            |> List.head
            |> fun psm -> psm.PercolatorScore
            )
        |> Map.ofList

    // Assigns a score to each protein with reverse digested peptides based on the peptides obtained in psm.
    let createReverseProteinScores (reverseProteins: (string*string[])[]) (peptideScoreMap: Map<string,float>) =
        reverseProteins
        |> Array.map (fun (protein, peptides) ->
        (protein |> String.split ' ').[0],
        peptides
        |> Array.map (fun pep ->
        peptideScoreMap.TryFind pep
        |> (fun x ->
            match x with
            | Some score -> score
            | None -> 0.
            )
        )
        |> Array.sum
        )
        |> Array.filter (fun (protein, score) -> score <> 0.)
        |> Map.ofArray

    /// Sums up score of all peptide sequences
    let assignPeptideScores (peptideSequences : string []) (peptideScoreMap : Map<string,float>) =
        peptideSequences
        |> Array.map (fun sequence -> peptideScoreMap.Item sequence)
        |> Array.sum

    // Looks if the given protein accession is present in a map of identified decoy proteins and assigns its score when found.
    let assignDecoyScoreToTargetScore (proteins: string) (decoyScores: Map<string,float>) =
        let prots = proteins |> String.split ';'
        prots
        |> Array.map (fun protein ->
            decoyScores.TryFind protein
            |> fun protOpt ->
                match protOpt with
                | Some score -> score
                | None -> 0.
        )
        |> Array.max

    /// for given data, creates a logistic regression model and returns a mapping function for this model
    let getLogisticRegressionFunction (x:vector) (y:vector) epsilon =
        let alpha =
            match FSharp.Stats.Fitting.LogisticRegression.Univariable.estimateAlpha epsilon x y with
            | Some a -> a
            | None -> failwith "Could not find an alpha for logistic regression of fdr data"
        let weight = FSharp.Stats.Fitting.LogisticRegression.Univariable.coefficient epsilon alpha x y
        FSharp.Stats.Fitting.LogisticRegression.Univariable.fit weight

    /// returns scores, pep, q
    let binningFunction bandwidth pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[])  =
        let totalDecoyProportion =
            let decoyCount = Array.filter isDecoyF data |> Array.length |> float
            let totalCount = data |> Array.length  |> float
            (2. * decoyCount / totalCount)
        data
        |> Array.groupBy (fun s -> floor (scoreF s / bandwidth))
        |> Array.sortBy fst
        |> Array.map (fun (k,values)->
            let median     = values |> Array.map scoreF |> Array.average
            let totalCount = values |> Array.length |> float
            let decoyCount = values |> Array.filter isDecoyF |> Array.length |> float |> (*) totalDecoyProportion
            //(median |> float,(decoyCount * pi0  / totalCount))
            median,totalCount,decoyCount
                //(median, totalCount )
        )
        |> fun a ->
            a
            |> Array.mapi (fun i (median,totalCountBin,decoyCountBin) ->
                            /// TODO: Accumulate totalCount + totalDecoyCount beforeHand and skip the time intensive mapping accross the array in each iteration.
                            let _,totalCountRight,decoyCountRight = a.[i..a.Length-1] |> Array.reduce (fun (x,y,z) (x',y',z') -> x+x',y+y',z+z')
                            (median,(pi0 * 2. * decoyCountBin / totalCountBin),(pi0 * 2. * decoyCountRight / totalCountRight))
                          )
        |> Array.sortBy (fun (score,pep,q) -> score)
        |> Array.unzip3
        |> fun (score,pep,q) -> vector score, vector pep, vector q

    /// Calculates q value mapping funtion for target/decoy dataset
    let getQValueFunc pi0 bw (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) =
        let (scores,_,q) = binningFunction bw pi0 scoreF isDecoyF data
        getLogisticRegressionFunction scores q 0.0000001

    /// Calculates q values for target/decoy dataset
    let getQValues pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) =
        let f = getQValueFunc pi0 0.01 scoreF isDecoyF data
        Array.map (scoreF >> f) data

    // Calculates a q-value for the indentified proteins
    let calculateQValueLogReg fdrEstimate (targetDecoyMatch: InferredProteinClassItemScored<'sequence>[]) (decoyNoMatch: InferredProteinClassItemScored<'sequence>[]) =
        // Input for q value calculation
        let createTargetDecoyInput =
            targetDecoyMatch
            |> Array.map (fun inferredProteinCIS ->
                if inferredProteinCIS.DecoyBigger then
                    createQValueInput inferredProteinCIS.DecoyScore true
                else
                    createQValueInput inferredProteinCIS.TargetScore false
                )
        let createDecoyNoMatchInput =
            decoyNoMatch
            |> Array.map (fun intermediateResult ->
                createQValueInput intermediateResult.DecoyScore true
                )
        // Combined input for q value calculation
        let combinedInput = Array.append createTargetDecoyInput createDecoyNoMatchInput

        // Combined InferredProteinClassItemScored input, which gets assigned its corresponding q value
        let combinedIPCISInput = Array.append targetDecoyMatch decoyNoMatch

        let qValues = getQValues fdrEstimate (fun (x: QValueInput) -> x.Score) (fun (x: QValueInput) -> x.IsDecoy) combinedInput

        // Create a new instance of InferredProteinClassItemScored with q values assigned
        Array.map2 (fun (qValue: float) (input: InferredProteinClassItemScored<'sequence>) ->
            createInferredProteinClassItemScored input.GroupOfProteinIDs input.Class input.PeptideSequence input.TargetScore input.DecoyScore qValue input.Decoy input.DecoyBigger
            ) qValues combinedIPCISInput

    let calculateQValueStorey (targetDecoyMatch: InferredProteinClassItemScored<'sequence>[]) (decoyNoMatch: InferredProteinClassItemScored<'sequence>[]) =
        // Combined input for q value calculation
        let combinedInput = Array.append targetDecoyMatch decoyNoMatch
                            |> Array.sortByDescending (fun x -> if x.DecoyBigger then
                                                                    x.DecoyScore
                                                                else
                                                                    x.TargetScore
                                                      )
        let rec traverseTopBottom (i: int) (target: float) (decoy: float) (multiplier: float) (qValues: float list)=
            if i >= combinedInput.Length then
                qValues
            else
                match combinedInput.[i].DecoyBigger with
                | true  -> traverseTopBottom (i + 1) target (decoy + 1. * multiplier) multiplier (((decoy + 1. * multiplier) / (if target > 0. then target else 1.))::qValues)
                | false -> traverseTopBottom (i + 1) (target + 1.) decoy multiplier ((decoy / (target + 1.))::qValues)
        let reverseQVal = traverseTopBottom 0 0. 0. 1. []

        let rec traverseBottomTop (i: int) (monotonizedQValues: float list) =
            if i >= reverseQVal.Length then
                monotonizedQValues
            elif i = 0 then
                traverseBottomTop (i + 1) (reverseQVal.[i]::monotonizedQValues)
            else
                if reverseQVal.[i] > monotonizedQValues.[0] then
                    traverseBottomTop (i + 1) (monotonizedQValues.[0]::monotonizedQValues)
                else
                    traverseBottomTop (i + 1) (reverseQVal.[i]::monotonizedQValues)
        let monotoneQVal = traverseBottomTop 0 []
                           |> Array.ofList

        // Create a new instance of InferredProteinClassItemScored with q values assigned
        Array.map2 (fun (qValue: float) (input: InferredProteinClassItemScored<'sequence>) ->
            createInferredProteinClassItemScored input.GroupOfProteinIDs input.Class input.PeptideSequence input.TargetScore input.DecoyScore qValue input.Decoy input.DecoyBigger
            ) monotoneQVal combinedInput

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
                Seq.fromFileWithCsvSchema<PSMInput>(filePath, '\t', true,schemaMode = SchemaModes.Fill)
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
            let preparedProtPepFunc = getProteinPeptideLookUpFromFileBy memoryDB
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
           try
                GFF3.fromFileWithoutFasta gff3Path
                |> assignTranscriptsToGenes regexPattern
            with
            | err ->
                printfn "ERROR: Could not read gff3 file %s" gff3Path
                failwithf "%s" err.Message
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
            let ppRelationModel = createPeptideProteinRelation proteinModels

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

    let removeModification pepSeq =
        String.filter (fun c -> System.Char.IsLower c |> not && c <> '[' && c <> ']') pepSeq

    let proteinGroupToString (proteinGroup:string[]) =
        Array.reduce (fun x y ->  x + ";" + y) proteinGroup

    let readAndInferFile classItemCollection protein peptide groupFiles outDirectory rawFolderPath psmInputs (dbConnection: SQLiteConnection) (qValMethod: QValueMethod) =

        let logger = Logging.createLogger "ProteinInference_readAndInferFile"

        let rawFilePaths = System.IO.Directory.GetFiles (rawFolderPath, "*.qpsm")
                           |> Array.toList

        let outFiles: string list =
            rawFilePaths
            |> List.map (fun filePath ->
                let foldername = (rawFolderPath.Split ([|"\\"|], System.StringSplitOptions.None))
                outDirectory + @"\" + foldername.[foldername.Length - 1] + "\\" + (System.IO.Path.GetFileNameWithoutExtension filePath) + ".prot"
                )

        let dbParams = getSDBParamsBy dbConnection

        // Array of prtoein Accessions tupled with their reverse digested peptides
        let reverseProteins =
            (selectProteins dbConnection)
            |> Array.ofList
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
                    |> List.map (fun s ->
                        let s =
                            s.Seq

                        match Map.tryFind (removeModification s) classItemCollection with
                        | Some (x:ProteinClassItem<'sequence>) -> createProteinClassItem x.GroupOfProteinIDs x.Class s
                        | None -> failwithf "Could not find sequence %s in classItemCollection" s
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

            // Peptide score Map
            let peptideScoreMap = createPeptideScoreMap psmInputs

            // Assigns scores to reverse digested Proteins using the peptideScoreMap
            let reverseProteinScores = createReverseProteinScores reverseProteins peptideScoreMap

            // Scores each inferred protein and assigns each protein where a reverted peptide was hit its score
            let combinedScoredClasses =
                combinedClasses
                |> Array.ofSeq
                |> Array.map (fun inferredPCI ->
                    // Looks up all peptides assigned to the protein and sums up their score
                    let peptideScore = (assignPeptideScores inferredPCI.PeptideSequence peptideScoreMap)
                    // Looks if the reverse protein has been randomly matched and assigns the score
                    let decoyScore = (assignDecoyScoreToTargetScore inferredPCI.GroupOfProteinIDs reverseProteinScores)

                    createInferredProteinClassItemScored
                        inferredPCI.GroupOfProteinIDs inferredPCI.Class inferredPCI.PeptideSequence
                        peptideScore
                        decoyScore
                        // Placeholder for q value
                        -1.
                        false
                        (decoyScore > peptideScore)
                    )

            // creates InferredProteinClassItemScored type for decoy proteins that have no match
            let reverseNoMatch =
                let proteinsPresent =
                    combinedScoredClasses
                    |> Array.collect (fun prots -> prots.GroupOfProteinIDs |> String.split ';')
                reverseProteinScores
                |> Map.toArray
                |> Array.map (fun (protein, score) ->
                    if (proteinsPresent |> Array.contains protein) then
                        createInferredProteinClassItemScored "HasMatch" PeptideEvidenceClass.Unknown [|""|] (-1.) (-1.) (-1.) false false
                    else
                        createInferredProteinClassItemScored protein PeptideEvidenceClass.Unknown [|""|] 0. score (-1.) true true
                )
                |> Array.filter (fun out -> out.Decoy <> false)

            // Assign q values to each protein (now also includes decoy only hits)
            let combinedScoredClassesQVal =
                if qValMethod = Domain.QValueMethod.LogisticRegression then
                    let fdr = (combinedScoredClasses |> Array.filter (fun x -> x.DecoyBigger) |> Array.length |> float |> (*)2.)/(combinedScoredClasses |> Array.filter (fun x -> not x.DecoyBigger) |> Array.length |> float)
                    calculateQValueLogReg fdr combinedScoredClasses reverseNoMatch
                else
                    calculateQValueStorey combinedScoredClasses reverseNoMatch

            let graph =
                let decoy, target = combinedScoredClassesQVal |> Array.partition (fun x -> x.DecoyBigger)
                // Histogram with relative abundance
                let freqTarget = FSharp.Stats.Distributions.Frequency.create 0.01 (target |> Array.map (fun x -> x.TargetScore))
                                 |> Map.toArray
                                 |> Array.map (fun x -> fst x, (float (snd x)) / (float target.Length))
                let freqDecoy  = FSharp.Stats.Distributions.Frequency.create 0.01 (decoy |> Array.map (fun x -> x.DecoyScore))
                                 |> Map.toArray
                                 |> Array.map (fun x -> fst x, (float (snd x)) / (float target.Length))
                // Histogram with absolute values
                let freqTarget1 = FSharp.Stats.Distributions.Frequency.create 0.01 (target |> Array.map (fun x -> x.TargetScore))
                                 |> Map.toArray
                let freqDecoy1  = FSharp.Stats.Distributions.Frequency.create 0.01 (decoy |> Array.map (fun x -> x.DecoyScore))
                                 |> Map.toArray
                let histogram =
                    [
                        Chart.Column freqTarget |> Chart.withTraceName "Target"
                            |> Chart.withAxisAnchor(Y=1);
                        Chart.Column freqDecoy |> Chart.withTraceName "Decoy"
                            |> Chart.withAxisAnchor(Y=1);
                        Chart.Column freqTarget1
                            |> Chart.withAxisAnchor(Y=2)
                            |> Chart.withMarkerStyle (Opacity = 0.)
                            |> Chart.withTraceName (Showlegend = false);
                        Chart.Column freqDecoy1
                            |> Chart.withAxisAnchor(Y=2)
                            |> Chart.withMarkerStyle (Opacity = 0.)
                            |> Chart.withTraceName (Showlegend = false)
                    ]
                    |> Chart.Combine

                let sortedQValues = 
                    combinedScoredClassesQVal 
                    |> Array.map 
                        (fun x -> if x.DecoyBigger then
                                    x.DecoyScore, x.QValue
                                  else
                                    x.TargetScore, x.QValue
                        )
                    |> Array.sortBy (fun (score, qVal) -> score)

                [
                    Chart.Line sortedQValues |> Chart.withTraceName "Q-Values";
                    histogram
                ]
                |> Chart.Combine
                |> Chart.withY_AxisStyle("Relative Frequency / Q-Value",Side=StyleParam.Side.Left,Id=1, MinMax = (0., 1.))
                |> Chart.withY_AxisStyle("Absolute Frequency",Side=StyleParam.Side.Right,Id=2,Overlaying=StyleParam.AxisAnchorId.Y 1, MinMax = (0., float target.Length))
                |> Chart.withX_AxisStyle "Score"
                |> Chart.withSize (900., 900.)
                |> Chart.SaveHtmlAs (outDirectory + @"\QValueGraph")

            // Assign results to files in which they can be found
            classifiedProteins
            |> List.iter (fun (prots,outFile) ->
                let pepSeqSet = prots |> List.map (fun x -> x.PeptideSequence) |> Set.ofList
                logger.Trace (sprintf "Start with %s"(System.IO.Path.GetFileNameWithoutExtension outFile))
                let combinedInferenceresult =
                    combinedScoredClassesQVal |> Array.filter (fun inferredPCIS -> not inferredPCIS.Decoy)
                    |> Seq.choose (fun ic ->
                        let filteredPepSet =
                            ic.PeptideSequence
                            |> Array.filter (fun pep -> Set.contains pep pepSeqSet)
                        if filteredPepSet = [||] then
                            None
                        else
                            Some (createInferredProteinClassItemScored ic.GroupOfProteinIDs ic.Class [|proteinGroupToString filteredPepSet|] ic.TargetScore ic.DecoyScore ic.QValue ic.Decoy ic.DecoyBigger)
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

                // Peptide Score Map
                let peptideScoreMap = createPeptideScoreMap [psm]

                // Assigns scores to reverse digested Proteins using the peptideScoreMap
                let reverseProteinScores = createReverseProteinScores reverseProteins peptideScoreMap

                // Scores each inferred protein and assigns each protein where a reverted peptide was hit its score
                let inferenceResultScored =
                    inferenceResult
                    |> Array.ofSeq
                    |> Array.map (fun inferredPCI ->
                        // Looks up all peptides assigned to the protein and sums up their score
                        let peptideScore = (assignPeptideScores inferredPCI.PeptideSequence peptideScoreMap)
                        // Looks if the reverse protein has been randomly matched and assigns the score
                        let decoyScore = (assignDecoyScoreToTargetScore inferredPCI.GroupOfProteinIDs reverseProteinScores)

                        createInferredProteinClassItemScored
                            inferredPCI.GroupOfProteinIDs inferredPCI.Class inferredPCI.PeptideSequence
                            peptideScore
                            decoyScore
                            // Placeholder for q value
                            -1.
                            false
                            (decoyScore > peptideScore)
                        )

                // creates InferredProteinClassItemScored type for decoy proteins that have no match
                let reverseNoMatch =
                    let proteinsPresent =
                        inferenceResultScored
                        |> Array.collect (fun prots -> prots.GroupOfProteinIDs |> String.split ';')
                    reverseProteinScores
                    |> Map.toArray
                    |> Array.map (fun (protein, score) ->
                        if (proteinsPresent |> Array.contains protein) then
                            createInferredProteinClassItemScored "HasMatch" PeptideEvidenceClass.Unknown [|""|] (-1.) (-1.) (-1.) false false
                        else
                            createInferredProteinClassItemScored protein PeptideEvidenceClass.Unknown [|""|] 0. score (-1.) true true
                    )
                    |> Array.filter (fun out -> out.Decoy <> false)

                // Assign q values to each protein (now also includes decoy only hits)
                let inferenceResultScoredQVal =
                    if qValMethod = Domain.QValueMethod.LogisticRegression then
                        let fdr = (inferenceResultScored |> Array.filter (fun x -> x.DecoyBigger) |> Array.length |> float |> (*)2.)/(inferenceResultScored |> Array.filter (fun x -> not x.DecoyBigger) |> Array.length |> float)
                        calculateQValueLogReg fdr inferenceResultScored reverseNoMatch
                    else
                        calculateQValueStorey inferenceResultScored reverseNoMatch

                let graph =
                    let decoy, target = inferenceResultScoredQVal |> Array.partition (fun x -> x.DecoyBigger)
                    // Histogram with relative abundance
                    let freqTarget = FSharp.Stats.Distributions.Frequency.create 0.01 (target |> Array.map (fun x -> x.TargetScore))
                                     |> Map.toArray
                                     |> Array.map (fun x -> fst x, (float (snd x)) / (float target.Length))
                    let freqDecoy  = FSharp.Stats.Distributions.Frequency.create 0.01 (decoy |> Array.map (fun x -> x.DecoyScore))
                                     |> Map.toArray
                                     |> Array.map (fun x -> fst x, (float (snd x)) / (float target.Length))
                    // Histogram with absolute values
                    let freqTarget1 = FSharp.Stats.Distributions.Frequency.create 0.01 (target |> Array.map (fun x -> x.TargetScore))
                                     |> Map.toArray
                    let freqDecoy1  = FSharp.Stats.Distributions.Frequency.create 0.01 (decoy |> Array.map (fun x -> x.DecoyScore))
                                     |> Map.toArray
                    let histogram =
                        [
                            Chart.Column freqTarget |> Chart.withTraceName "Target"
                                |> Chart.withAxisAnchor(Y=1);
                            Chart.Column freqDecoy |> Chart.withTraceName "Decoy"
                                |> Chart.withAxisAnchor(Y=1);
                            Chart.Column freqTarget1
                                |> Chart.withAxisAnchor(Y=2)
                                |> Chart.withMarkerStyle (Opacity = 0.)
                                |> Chart.withTraceName (Showlegend = false);
                            Chart.Column freqDecoy1
                                |> Chart.withAxisAnchor(Y=2)
                                |> Chart.withMarkerStyle (Opacity = 0.)
                                |> Chart.withTraceName (Showlegend = false)
                        ]
                        |> Chart.Combine

                    let sortedQValues = 
                        inferenceResultScoredQVal 
                        |> Array.map 
                            (fun x -> if x.DecoyBigger then
                                        x.DecoyScore, x.QValue
                                      else
                                        x.TargetScore, x.QValue
                            )
                        |> Array.sortBy (fun (score, qVal) -> score)

                    [
                        Chart.Line sortedQValues |> Chart.withTraceName "Q-Values";
                        histogram
                    ]
                    |> Chart.Combine
                    |> Chart.withY_AxisStyle("Relative Frequency / Q-Value",Side=StyleParam.Side.Left,Id=1, MinMax = (0., 1.))
                    |> Chart.withY_AxisStyle("Absolute Frequency",Side=StyleParam.Side.Right,Id=2,Overlaying=StyleParam.AxisAnchorId.Y 1, MinMax = (0., float target.Length))
                    |> Chart.withX_AxisStyle "Score"
                    |> Chart.withSize (900., 900.)
                    |> Chart.SaveHtmlAs (outFile + @"_QValueGraph")


                inferenceResultScoredQVal
                |> Array.filter (fun inferredPCIS -> not inferredPCIS.Decoy)
                |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
                |> Seq.write outFile
            ) psmInputs

    let inferProteins gff3Location dbConnection (proteinInferenceParams: ProteinInferenceParams) outDirectory rawFolderPath =

        let logger = Logging.createLogger "ProteinInference_inferProteins"

        logger.Trace (sprintf "InputFilePath = %s" rawFolderPath)
        logger.Trace (sprintf "InputGFF3Path = %s" gff3Location)
        logger.Trace (sprintf "OutputFilePath = %s" outDirectory)
        logger.Trace (sprintf "Protein inference parameters = %A" proteinInferenceParams)

        logger.Trace "Copy peptide DB into Memory."
        let memoryDB = SearchDB.copyDBIntoMemory dbConnection
        logger.Trace "Copy peptide DB into Memory: finished."

        logger.Trace "Start building ClassItemCollection"
        let classItemCollection, psmInputs = createClassItemCollection gff3Location memoryDB proteinInferenceParams.ProteinIdentifierRegex rawFolderPath
        logger.Trace "Classify and Infer Proteins"
        readAndInferFile classItemCollection proteinInferenceParams.Protein proteinInferenceParams.Peptide 
                         proteinInferenceParams.GroupFiles outDirectory rawFolderPath psmInputs dbConnection proteinInferenceParams.QValueMethod