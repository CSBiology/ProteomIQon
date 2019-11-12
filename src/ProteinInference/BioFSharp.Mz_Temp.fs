namespace ProteomIQon

open BioFSharp
open BioFSharp.IO
open System.Data
open System.Data.SQLite
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open FSharp.Stats
open FSharp.Stats.Fitting.NonLinearRegression
open System
open FSharpAux
open AminoAcids
open ModificationInfo
open ProteomIQon



module Fitting' =

    module NonLinearRegression' =
        
        module LevenbergMarquardtConstrained' =
            
            /// Returns an estimate for an initial parameter for the linear least square estimator for a given dataset (xData, yData).
            /// The initial estimation is intended for a logistic function.
            let initialParam (xData: float[]) (yData: float[]) =
                let xRange = abs ((xData |> Array.max) - (xData |> Array.min))
                let yRange = abs ((yData |> Array.max) - (yData |> Array.min))
                let maxY = yData |> Array.max
                let combined = Array.map2 (fun x y -> x,y) xData yData
                // Looks for the real point that is closest to the given point
                let rec findClosest i (point: float) (data: float []) (distance: float) (range: float)=
                    let newDistance = abs (data.[i] - point)
                    // Checks if point is located in the middle of the range of possible points
                    if distance < newDistance && (data.[i - 1] < (point + 0.25 * range) && data.[i - 1] > (point - 0.25 * range)) then
                        data.[i - 1]
                    else
                        findClosest (i + 1) point data newDistance range
                let midX,midY =
                    let point = maxY - yRange / 2.
                    let middleYData = findClosest 1 point yData (abs (yData.[0] - point)) yRange
                    Array.filter (fun (x,y) -> y = middleYData) combined
                    |> Array.averageBy fst, middleYData
                let rightSlopeX,rightSlopeY =
                    combined
                    |> Array.filter (fun (x, y) -> (maxY - y) < 0.001 * yRange)
                    |> Array.head
                let leftSlopeX, leftSlopeY =
                    let leftX = midX - (rightSlopeX - midX)
                    leftX, maxY
                let slope =
                    ((rightSlopeY - leftSlopeY)/yRange) / ((rightSlopeX - leftSlopeX)/xRange)
                let steepness = abs slope
                Table.lineSolverOptions [|maxY; steepness; midX|]

            /// Returns an estimate for an initial parameter for the linear least square estimator for a given dataset (xData, yData).
            /// The steepness is given as an array and not estimated. An initial estimate is returned for every given steepness.
            /// The initial estimation is intended for a logistic function.
            let initialParamsOverRange (xData: float[]) (yData: float[]) (steepnessRange: float []) =
                let yRange = abs ((yData |> Array.max) - (yData |> Array.min))
                let maxY = yData |> Array.max
                let combined = Array.map2 (fun x y -> x,y) xData yData
                // Looks for the real point that is closest to the given point
                let rec findClosest i (point: float) (data: float []) (distance: float) (range: float)=
                    let newDistance = abs (data.[i] - point)
                    // Checks if point is located in the middle of the range of possible points
                    if distance < newDistance && (data.[i - 1] < (point + 0.25 * range) && data.[i - 1] > (point - 0.25 * range)) then
                        data.[i - 1]
                    else
                        findClosest (i + 1) point data newDistance range
                let midX,midY =
                    let point = maxY - yRange / 2.
                    let middleYData = findClosest 1 point yData (abs (yData.[0] - point)) yRange
                    Array.filter (fun (x,y) -> y = middleYData) combined
                    |> Array.averageBy fst, middleYData
                steepnessRange
                |> Array.map (fun steepness ->
                    Table.lineSolverOptions [|maxY; steepness; midX|]
                )

            /// Returns a parameter vector tupled with its residual sum of squares (RSS) as a possible solution for linear least square based nonlinear fitting of a given dataset (xData, yData) with a given 
            /// model function. 
            let estimatedParamsWithRSS (model: Model) (solverOptions: SolverOptions) lambdaInitial lambdaFactor (lowerBound: vector) (upperBound: vector) (xData: float[]) (yData: float []) =
                let estParams = LevenbergMarquardtConstrained.estimatedParamsVerbose model solverOptions lambdaInitial lambdaFactor lowerBound upperBound xData yData
                estParams
                |> fun estParams ->
                    let paramGuess = Vector.ofArray solverOptions.InitialParamGuess
                    let rss = getRSS model xData yData paramGuess
                    estParams.[estParams.Count-1], rss
            
    

module SearchDB' =


    module DB' =

        module SQLiteQuery' =
            
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
        let selectProteinByProtID         = DB'.SQLiteQuery'.prepareSelectProteinAccessionByID memoryDB tr
        let selectPeptideByPepSeqID       = DB'.SQLiteQuery'.prepareSelectPepSequenceByPepSequenceID memoryDB tr
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

module FDRControl' =

    open FSharp.Stats

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
            1. / (2. * decoyCount / totalCount)
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
        let (scores,_,q) = FDRControl.binningFunction bw pi0 scoreF isDecoyF data
        FDRControl.getLogisticRegressionFunction scores q 0.0000001

    /// Calculates q values for target/decoy dataset
    let getQValues pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) =
        let f = getQValueFunc pi0 0.01 scoreF isDecoyF data
        Array.map (scoreF >> f) data

module ProteinInference' =

    open BioFSharp.PeptideClassification
    open BioFSharp.IO.GFF3
    open FSharpAux.IO.SchemaReader.Attribute
    open FSharp.Stats.Interpolation
    open FSharp.Plotly

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

    type ScoreTargetDecoyCount =
        {
            Score      : float
            DecoyCount : float
            TargetCount: float
        }

    let createScoreTargetDecoyCount score decoyCount targetCount =
        {
            Score       = score
            DecoyCount  = decoyCount
            TargetCount = targetCount
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

        let qValues = FDRControl'.getQValueFunc fdrEstimate 0.01 (fun (x: QValueInput) -> x.Score) (fun (x: QValueInput) -> x.IsDecoy) combinedInput

        // Create a new instance of InferredProteinClassItemScored with q values assigned
        let combinedInputQVal =
            combinedIPCISInput
            |> Array.map (fun item ->
                if item.Decoy then
                    createInferredProteinClassItemScored item.GroupOfProteinIDs item.Class item.PeptideSequence item.TargetScore item.DecoyScore (qValues item.DecoyScore) item.Decoy item.DecoyBigger
                else
                    createInferredProteinClassItemScored item.GroupOfProteinIDs item.Class item.PeptideSequence item.TargetScore item.DecoyScore (qValues item.TargetScore) item.Decoy item.DecoyBigger
            )
        combinedInputQVal

    /// Calculates the q value using Storeys method and the paired target-decoy approach to determine target or decoy hits.
    let calculateQValueStorey (targetDecoyMatch: InferredProteinClassItemScored<'sequence>[]) (decoyNoMatch: InferredProteinClassItemScored<'sequence>[]) =
        // Combined input for q value calculation
        // Gives an array of scores with the frequency of decoy and target hits at that score
        let combinedInput = Array.append targetDecoyMatch decoyNoMatch
        let scoreFrequencies =
            combinedInput
            |> Array.map (fun x -> if x.DecoyBigger then
                                    x.DecoyScore, true
                                    else
                                    x.TargetScore, false
            )
            |> Array.groupBy fst
            |> Array.map (fun (score,scoreDecoyInfo) ->
                let decoyCount = 
                    Array.filter (fun (score, decoyInfo) -> decoyInfo = true) scoreDecoyInfo
                    |> Array.length
                let targetCount = 
                    Array.filter (fun (score, decoyInfo) -> decoyInfo = false) scoreDecoyInfo
                    |> Array.length
                createScoreTargetDecoyCount score (float decoyCount) (float targetCount)
            )
            |> Array.sortByDescending (fun x -> x.Score)

        // Goes through the list and assigns each protein a "q value" by dividing total decoy hits so far through total target hits so far
        let reverseQVal =
            scoreFrequencies
            |> Array.fold (fun (acc: (float*float*float*float) list) scoreCounts ->
                let _,_,decoyCount,targetCount = acc.Head
                // Decoy hits are doubled
                let newDecoyCount  = decoyCount + scoreCounts.DecoyCount * 2.
                let newTargetCount = targetCount + scoreCounts.TargetCount
                let newQVal = 
                    let nominator = 
                        if newTargetCount > 0. then
                            newTargetCount 
                        else 1.
                    newDecoyCount / nominator
                (scoreCounts.Score, newQVal, newDecoyCount, newTargetCount):: acc
            ) [0., 0., 0., 0.]
            |> fun list -> list.[.. list.Length-2]
            |> List.map (fun (score, qVal, decoyC, targetC) -> score, qVal)

        //Assures monotonicity by going through the list from the bottom to top and assigning the previous q value if it is smaller than the current one
        let score, monotoneQVal =
            let head::tail = reverseQVal
            tail
            |> List.fold (fun (acc: (float*float) list) (score, newQValue) ->
                let _,qValue = acc.Head
                if newQValue > qValue then
                    (score, qValue)::acc
                else
                    (score, newQValue)::acc

            )[head]
            |> Array.ofList
            |> Array.sortBy fst
            |> Array.unzip
        // Linear Interpolation
        let linearSplineCoeff = LinearSpline.initInterpolateSorted score monotoneQVal
        let interpolation = LinearSpline.interpolate linearSplineCoeff
        // Assigns q values based on scores with interpolation function
        let combinedInputQVal =
            combinedInput
            |> Array.map (fun item ->
                if item.Decoy then
                    createInferredProteinClassItemScored item.GroupOfProteinIDs item.Class item.PeptideSequence item.TargetScore item.DecoyScore (interpolation item.DecoyScore) item.Decoy item.DecoyBigger
                else
                    createInferredProteinClassItemScored item.GroupOfProteinIDs item.Class item.PeptideSequence item.TargetScore item.DecoyScore (interpolation item.TargetScore) item.Decoy item.DecoyBigger
            )
        combinedInputQVal

    let qValueHitsVisualization inferredProteinClassItemScored path=
        let decoy, target = inferredProteinClassItemScored |> Array.partition (fun x -> x.DecoyBigger)
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
            inferredProteinClassItemScored 
            |> Array.map 
                (fun x -> if x.Decoy then
                            x.DecoyScore, x.QValue
                            else
                            x.TargetScore, x.QValue
                )
            |> Array.sortBy (fun (score, qVal) -> score)

        [
            Chart.Point sortedQValues |> Chart.withTraceName "Q-Values";
            histogram
        ]
        |> Chart.Combine
        |> Chart.withY_AxisStyle("Relative Frequency / Q-Value",Side=StyleParam.Side.Left,Id=1, MinMax = (0., 1.))
        |> Chart.withY_AxisStyle("Absolute Frequency",Side=StyleParam.Side.Right,Id=2,Overlaying=StyleParam.AxisAnchorId.Y 1, MinMax = (0., float target.Length))
        |> Chart.withX_AxisStyle "Score"
        |> Chart.withSize (900., 900.)
        |> Chart.SaveHtmlAs (path + @"\QValueGraph")


