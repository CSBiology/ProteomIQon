namespace ProteomIQon

open BioFSharp
open BioFSharp.IO
open System.Data
open System.Data.SQLite
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open FSharp.Stats
open FSharp.Stats.Fitting.NonLinearRegression
open FSharpAux
open ModificationInfo

module Fitting' =

    module NonLinearRegression' =

        module LevenbergMarquardtConstrained' =

            /// Logistic Function
            /// Line model of the form "y = a * x + b"
            let LogisticFunction = {
                ParameterNames= [|"L - curve maximum";"k - Steepness"; "x0 xValue of midpoint"|]
                GetFunctionValue = (fun (parameterVector:Vector<float>) xValue -> parameterVector.[0] / (1. + exp(parameterVector.[1]*(xValue-parameterVector.[2]))))
                GetGradientValue = (fun (parameterVector:Vector<float>) (gradientVector: Vector<float>) xValue ->
                                    gradientVector.[0] <- 1. / (1. + exp(parameterVector.[1]*(xValue-parameterVector.[2])))
                                    gradientVector.[1] <- (parameterVector.[0] * (xValue-parameterVector.[2]) * exp(parameterVector.[1]*(xValue-parameterVector.[2])) ) / (exp(parameterVector.[1]*(xValue-parameterVector.[2])) + 1.)**2.
                                    gradientVector.[2] <- (parameterVector.[0] * parameterVector.[1] * exp(parameterVector.[1]*(xValue-parameterVector.[2])) ) / (exp(parameterVector.[1]*(xValue-parameterVector.[2])) + 1.)**2.
                                    gradientVector)
                }

            /// Returns an estimate for an initial parameter for the linear least square estimator for a given dataset (xData, yData).
            /// The initial estimation is intended for a logistic function.
            let initialParam (xData: float[]) (yData: float[]) =
                let xRange = abs ((xData |> Array.max) - (xData |> Array.min))
                let yRange = abs ((yData |> Array.max) - (yData |> Array.min))
                let maxY = yData |> Array.max
                let combined = Array.map2 (fun x y -> x,y) xData yData
                // Looks for the real point that is closest to the given point
                let findClosest (point: float) (data: float []) =
                    let distance =
                        data
                        |> Array.map (fun x ->
                            abs (point - x)
                        )
                    let indexSmallest =
                        distance
                        |> Array.findIndex (fun x ->
                            x = (distance |> Array.min)
                        )
                    data.[indexSmallest]
                let midX,midY =
                    let point = maxY - yRange / 2.
                    let middleYData = findClosest point yData
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
                let findClosest (point: float) (data: float []) =
                    let distance =
                        data
                        |> Array.map (fun x ->
                            abs (point - x)
                        )
                    let indexSmallest =
                        distance
                        |> Array.findIndex (fun x ->
                            x = (distance |> Array.min)
                        )
                    data.[indexSmallest]
                let midX,midY =
                    let point = maxY - yRange / 2.
                    let middleYData = findClosest point yData
                    printfn "%f" middleYData
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


module ProteinInference' =

    open BioFSharp.PeptideClassification
    open BioFSharp.IO.GFF3
    open FSharpAux.IO.SchemaReader.Attribute
    open FSharp.Stats.Interpolation
    open FSharp.Plotly
    open Fitting'.NonLinearRegression'.LevenbergMarquardtConstrained'

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
            FoundInDB        : bool
        }

    let createInferredProteinClassItemScored proteinIDs evidenceClass peptideSequences targetScore decoyScore qValue isDecoy decoyBigger foundInDB =
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence   = peptideSequences
            Class             = evidenceClass
            TargetScore       = targetScore
            DecoyScore        = decoyScore
            QValue            = qValue
            Decoy             = isDecoy
            DecoyBigger       = decoyBigger
            FoundInDB         = foundInDB
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

module FDRControl' =

    open FSharp.Plotly
    open Fitting'.NonLinearRegression'.LevenbergMarquardtConstrained'
    open FSharp.Stats.Interpolation

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
            //let totalCount = values |> Array.length |> float
            let decoyCount = values |> Array.filter isDecoyF |> Array.length |> float |> (*) totalDecoyProportion
            // Include modified decoy count in total count?
            let totalCount = values |> Array.filter (isDecoyF >> not) |> Array.length |> float |> (+) decoyCount
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

// FDR estimation using MAYU

    let stirling_log_factorial (n: float) =
        log(sqrt(2. * pi  *n)) + n * log(n) - n

    let exact_log_factorial (n: float) =
        let rec loop i log_fact =
            if i > n then
                log_fact
            else
                let new_log_fact = log_fact + (log i)
                loop (i + 1.) new_log_fact
        loop 2. 0.

    let log_factorial (n: float) =
        if n < 1000. then
            exact_log_factorial n
        else
            stirling_log_factorial n

    let log_binomial (n: float) (k: float) =
        (log_factorial n) - (log_factorial k) - (log_factorial (n - k))

    let hypergeometric (x: float) (n: float) (w: float) (d: float) =
        //natural logarithm of the probability
        if(d > 0.) then
            exp((log_binomial w x) + (log_binomial (n-w) (d-x)) - (log_binomial n d))
        else 0.

    let estimatePi0HG (n: float) (targets: float) (cf: float) =
        let rec loop (fp: float) (logprob: float list) =
            if fp > cf then
                logprob |> List.rev
            else
                let tp = targets - fp
                let w = n - tp
                let prob = hypergeometric fp n w cf
                loop (fp + 1.) (prob::logprob)
        let logprob = loop 0. []
        let sum = logprob |> List.sum
        let logprob_Norm =
            logprob
            |> List.map (fun x ->
                x / sum
            )
        //// MAYU rounds here to first decimal
        let expectation_value_FP_PID =
            logprob_Norm
            |> List.foldi (fun i acc x ->
                acc + x * (float i)
            ) 0.
        if (isNan expectation_value_FP_PID) || (isInf expectation_value_FP_PID) then
            0.
        else 
            expectation_value_FP_PID

    let binProteinsLength (proteins: ProteinInference'.InferredProteinClassItemScored<'sequence> []) (proteinsFromDB: (string*string)[]) binCount =
        // need to bin all proteins in the db, not only those with hit?
        let proteinsNotMatchedDB =
            let proteinsMatched =
                proteins
                |> Array.collect (fun protein ->
                    protein.GroupOfProteinIDs
                    |> String.split ';'
                )
                |> Set.ofArray
            let proteinsInDB =
                proteinsFromDB
                |> Array.map fst
                |> Set.ofArray
            Set.difference proteinsMatched proteinsInDB
            |> Set.toArray
            |> Array.map (fun protein ->
                ProteinInference'.createInferredProteinClassItemScored protein BioFSharp.PeptideClassification.PeptideEvidenceClass.Unknown [|""|] (-1.) (-1.) (-1.) false false false
            )

        let combined = Array.append proteins proteinsNotMatchedDB
        // sorting treats protein groups as one protein with average length. They are also treated as one protein for total and target counts.
        let sortLength =
            combined
            |> Array.sortBy (fun protein ->
                let proteinsInGroup =
                    protein.GroupOfProteinIDs
                    |> String.split ';'
                let averageLength =
                    proteinsInGroup
                    |> Array.map (fun string ->
                        float string.Length
                    )
                    |> Array.average
                averageLength
            )
        let bins =
            let binSize =
                ceil (float sortLength.Length / float binCount)
            sortLength
            |> Array.chunkBySize (int binSize)
        bins

    // The original paper of Mayu describes a protein as:
    // FP = if only if all its peptides with q <= threshold are decoy
    // TP = at least one of its peptides with q <= threshold is target
    // However, the way mayu estimates it on the program is like this:
    // FP = any protein that contains a decoy psm with q <= threshold
    // TP = any protein that contains a target psm with q <= threshold
    // They do not consider protein containing both decoy and target psms.
    // ProteomIQon currently uses the picked target decoy approach with the following definitions:
    // FP = protein where the score of decoy hits is higher than score of target hits
    // TP = protein where the score of target hits is higher than score of decoy hits
    // Also, mayu estimates q as the empirical (target-decoy) q value.
    // Percolator estimates q as the empirical (target-decoy) q value and adjusted by pi0
    // Mayu extracts the list of TP and FP proteins from PSM level whereas percolator
    // extract the list of TP and FP proteins from peptide level, this avoids redundancy and
    // gives a better calibration since peptide level q values are re-adjusted in percolator.
    // ProteomIQon extracts TP and FP proteins from the result of the picked target decoy approach.
    // This creates sometimes a difference in the number of TP and FP proteins between percolator, Mayu, and ProteomIQon,
    // which causes a slight difference in the estimated protein FDR.

    let expectedFP (proteinBins: ProteinInference'.InferredProteinClassItemScored<'sequence> [] []) =
        proteinBins
        |> Array.map (fun proteinBin ->
            let numberTarget =
                proteinBin
                |> Array.filter (fun protein -> not protein.DecoyBigger && protein.FoundInDB) 
                |> Array.length
                |> float
            let numberDecoy = 
                proteinBin
                |> Array.filter (fun protein -> protein.DecoyBigger && protein.FoundInDB) 
                |> Array.length
                |> float
            let total = 
                let notFound = 
                    proteinBin
                    |> Array.filter (fun protein -> not protein.FoundInDB)
                    |> Array.length 
                    |> float
                notFound + numberTarget + numberDecoy
            let fpInBin = estimatePi0HG total numberTarget numberDecoy
            printfn "target: %f; decoy: %f; total: %f" numberTarget numberDecoy total
            fpInBin
        )
        |> Array.sum


    // Calculates a q-value for the indentified proteins
    let calculateQValueLogReg fdrEstimate (targetDecoyMatch: ProteinInference'.InferredProteinClassItemScored<'sequence>[]) (decoyNoMatch: ProteinInference'.InferredProteinClassItemScored<'sequence>[]) =
        // Input for q value calculation
        let createTargetDecoyInput =
            targetDecoyMatch
            |> Array.map (fun inferredProteinCIS ->
                if inferredProteinCIS.DecoyBigger then
                    ProteinInference'.createQValueInput inferredProteinCIS.DecoyScore true
                else
                    ProteinInference'.createQValueInput inferredProteinCIS.TargetScore false
            )
        let createDecoyNoMatchInput =
            decoyNoMatch
            |> Array.map (fun intermediateResult ->
                ProteinInference'.createQValueInput intermediateResult.DecoyScore true
            )
        // Combined input for q value calculation
        let combinedInput = Array.append createTargetDecoyInput createDecoyNoMatchInput

        // Combined InferredProteinClassItemScored input, which gets assigned its corresponding q value
        let combinedIPCISInput = Array.append targetDecoyMatch decoyNoMatch

        let scores,pep,qVal =
            binningFunction 0.01 fdrEstimate (fun (x: ProteinInference'.QValueInput) -> x.Score) (fun (x: ProteinInference'.QValueInput) -> x.IsDecoy) combinedInput
            |> fun (scores,pep,qVal) -> scores.ToArray(), pep.ToArray(), qVal.ToArray()

        Chart.Point (scores,qVal)
        |> Chart.Show

        let initialGuess =
            initialParamsOverRange scores qVal [|1. .. 30.|]

        let estimate =
            initialGuess
            |> Array.map (fun initial ->
                if initial.InitialParamGuess.Length > 3 then failwith "Invalid initial param guess for Logistic Function"
                let lowerBound =
                    initial.InitialParamGuess
                    |> Array.map (fun param -> param - param * 0.1)
                    |> vector
                let upperBound =
                    initial.InitialParamGuess
                    |> Array.map (fun param -> param + param * 0.1)
                    |> vector
                estimatedParamsWithRSS LogisticFunction initial 0.001 10.0 lowerBound upperBound scores qVal
            )
            |> Array.minBy snd
            |> fst

        let logisticFunction = LogisticFunction.GetFunctionValue estimate

        //let qValues = FDRControl'.getQValueFunc fdrEstimate 0.01 (fun (x: QValueInput) -> x.Score) (fun (x: QValueInput) -> x.IsDecoy) combinedInput

        // Create a new instance of InferredProteinClassItemScored with q values assigned
        let combinedInputQVal =
            combinedIPCISInput
            |> Array.map (fun item ->
                if item.Decoy then
                    ProteinInference'.createInferredProteinClassItemScored item.GroupOfProteinIDs item.Class item.PeptideSequence item.TargetScore item.DecoyScore (logisticFunction item.DecoyScore) item.Decoy item.DecoyBigger true
                else
                    ProteinInference'.createInferredProteinClassItemScored item.GroupOfProteinIDs item.Class item.PeptideSequence item.TargetScore item.DecoyScore (logisticFunction item.TargetScore) item.Decoy item.DecoyBigger true
            )
        combinedInputQVal

    /// Calculates the q value using Storeys method and the paired target-decoy approach to determine target or decoy hits.
    let calculateQValueStorey (targetDecoyMatch: ProteinInference'.InferredProteinClassItemScored<'sequence>[]) (decoyNoMatch: ProteinInference'.InferredProteinClassItemScored<'sequence>[]) =
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
                ProteinInference'.createScoreTargetDecoyCount score (float decoyCount) (float targetCount)
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
                    ProteinInference'.createInferredProteinClassItemScored item.GroupOfProteinIDs item.Class item.PeptideSequence item.TargetScore item.DecoyScore (interpolation item.DecoyScore) item.Decoy item.DecoyBigger true
                else
                    ProteinInference'.createInferredProteinClassItemScored item.GroupOfProteinIDs item.Class item.PeptideSequence item.TargetScore item.DecoyScore (interpolation item.TargetScore) item.Decoy item.DecoyBigger true
            )
        combinedInputQVal
