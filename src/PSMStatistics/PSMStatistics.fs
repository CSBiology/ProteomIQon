namespace ProteomIQon

open System
open System.Data.SQLite
open BioFSharp
open BioFSharp.Mz.SearchDB
open Domain
open Core 
open System.IO
open BioFSharp.Mz
open FSharpAux.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Csv
open FSharpAux.IO.SchemaReader.Attribute
open System.Reflection
open Plotly.NET
open Microsoft
open Microsoft.ML
open Microsoft.ML.Data 

module PSMStatistics = 
    open System.IO
    ///
    let downcastPipeline (x : IEstimator<_>) = 
        match x with 
        | :? IEstimator<ITransformer> as y -> y
        | _ -> failwith "downcastPipeline: expecting a IEstimator<ITransformer>"

    [<CLIMutable>]
    type PSMToLearn = {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
        [<ColumnName("PSMId")>]
        PSMId                        : string
        [<ColumnName("Label")>]
        Label                        : bool
        [<ColumnName("ScanNr")>]
        ScanNr                       : int
        [<ColumnName("Charge")>]
        Charge                       : int//float32 []
        [<ColumnName("PrecursorMZ")>]
        PrecursorMZ                  : float32
        [<ColumnName("TheoMass")>]
        TheoMass                     : float32
        [<ColumnName("AbsDeltaMass")>]
        AbsDeltaMass                 : float32
        [<ColumnName("PeptideLength")>]
        PeptideLength                : float32
        [<ColumnName("MissCleavages")>]
        MissCleavages                : float32
        [<ColumnName("SequestScore")>]
        SequestScore                 : float32
        [<ColumnName("SequestNormDeltaBestToRest")>]
        SequestNormDeltaBestToRest   : float32
        [<ColumnName("SequestNormDeltaNext")>]
        SequestNormDeltaNext         : float32
        [<ColumnName("AndroScore")>]
        AndroScore                   : float32
        [<ColumnName("AndroNormDeltaBestToRest")>]
        AndroNormDeltaBestToRest     : float32
        [<ColumnName("AndroNormDeltaNext")>]
        AndroNormDeltaNext           : float32
        [<ColumnName("XtandemScore")>]
        XtandemScore                 : float32
        [<ColumnName("XtandemNormDeltaBestToRest")>]
        XtandemNormDeltaBestToRest   : float32
        [<ColumnName("XtandemNormDeltaNext")>]
        XtandemNormDeltaNext         : float32
        /// Hyperscore of PeptideSpectrumMatchingTIMs, 0 when the input has none.
        [<ColumnName("Hyperscore")>]
        Hyperscore                   : float32
        /// Negative log10 of the expectation value of PeptideSpectrumMatchingTIMs, 0 when the input has none.
        [<ColumnName("LogExpectscore")>]
        LogExpectscore               : float32
        /// Scan time minus the predicted scan time of the peptide, in robust standard deviations, 0 without RetentionTimeFeature.
        [<ColumnName("RtDeviation")>]
        RtDeviation                  : float32
        /// Absolute value of RtDeviation.
        [<ColumnName("AbsRtDeviation")>]
        AbsRtDeviation               : float32
        [<ColumnName("Peptide")>]
        Peptide                      : string
        [<ColumnName("Protein")>]
        Protein                      : string
        }

    [<CLIMutable>]
    type PSMPrediction = 
        { 
            // ColumnName attribute is used to change the column name from
            // its default value, which is the name of the field.
            [<ColumnName("PredictedLabel")>]
            PredictedLabel : bool; 

            // No need to specify ColumnName attribute, because the field
            // name "Probability" is the column name we want.
            Probability : float32; 
            Score : float32 
        }

    type AppliedModel = {
        NPositivesAtFDR : int
        CalcQValue      : float -> float
        CalcPepValue    : float -> float
        Model           : PSMToLearn -> PSMPrediction
        }
        

    type PSMStats = {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
        [<FieldAttribute(0)>]
        PSMId                        : string;
        [<FieldAttribute(1)>]
        ModelScore              : float
        [<FieldAttribute(2)>]
        QValue                       : float
        [<FieldAttribute(3)>]
        PosteriorErrorProbability    : float
        [<FieldAttribute(4)>]
        StringSequence               : string
        }
    
    /// Returns a list of proteins retrieved by PepsequenceID
    let initProteinAndClvIdxLookUp (cn:SQLiteConnection) tr =
        let selectCleavageIdxByPepSeqID   = Db.SQLiteQuery.prepareSelectCleavageIndexByPepSequenceID cn tr
        let selectProteinByProtID         = Db.SQLiteQuery.prepareSelectProteinByID cn tr
        (fun pepSequenceID ->
                selectCleavageIdxByPepSeqID pepSequenceID
                |> List.map (fun (clid,protID,peps,clvgS,clvgE,missClv) -> selectProteinByProtID protID,clvgS,clvgE,missClv )
        )
    
    /// Retrieves the flanking amino acids of a peptide Sequence originating of a given Protein Sequence
    let getFlankingAminoAcids missClvStart missClvEnd (protSequence:string) =
        let Nterminal = if missClvStart = 0 then "-" else protSequence.[missClvStart-1].ToString()
        let CTerminal = if missClvEnd = protSequence.Length-1 then "-" else  protSequence.[missClvEnd+1].ToString()
        Nterminal,CTerminal

    /// The residues of a sequence, each with the modification tokens written before it, such as "[ox]M". The terminator "*" and the gap "-" are dropped.
    let sequenceTokens (sequence: string) =
        Text.RegularExpressions.Regex.Matches(sequence, @"(\[[^\]]*\])*[^\[\]]")
        |> Seq.map (fun m -> m.Value)
        |> List.ofSeq
        |> List.filter (fun token -> token <> "*" && token <> "-")

    /// The decoy of a peptide sequence as the searches build it with
    /// Fragmentation.Series.pseudoReverse: the C-terminal residue stays in place and the rest is
    /// reversed, or the whole sequence if that gives the sequence itself. The terminator "*" and the
    /// gap "-" are dropped first, as the searches drop them. A residue keeps the
    /// modification tokens written before it, such as "[ox]M".
    let decoySequence (sequence: string) =
        sequenceTokens sequence
        |> Fragmentation.Series.pseudoReverse
        |> String.concat ""

    /// The score that starts the training: the negative log10 of the expectation value if the input has it, otherwise the SEQUEST-like score.
    let seedScoreOf useTimsScores (logExpectscore: float32) (sequestScore: float32) =
        if useTimsScores then logExpectscore else sequestScore

    /// Negative log10 of the expectation value of PeptideSpectrumMatchingTIMs, 0 when the input has none.
    let logExpectscoreOf (psm: Dto.PeptideSpectrumMatchingResult) =
        if Double.IsFinite psm.Expectscore && psm.Expectscore > 0. then float32 (-log10 psm.Expectscore) else 0.f

    /// Predicts the scan time of a peptide from its sequence. psmStats adds the deviation of the observed from the predicted time as model feature.
    module RetentionTime =

        let nFolds = 5
        /// Fewer aggregated peptides than this, or a lower out-of-fold R2, and the feature is not used.
        let minPeptides = 500
        let minR2 = 0.3

        let private residueIndex = "ACDEFGHIKLMNPQRSTVWY" |> Seq.mapi (fun i c -> c, i) |> dict
        let private nPairs = 210

        /// Index of the unordered pair of two residue indices.
        let private pairIndex a b =
            let lo, hi = min a b, max a b
            20 * lo - lo * (lo - 1) / 2 + hi - lo

        let private residueOf (token: string) = token.[token.Length - 1]

        let private modificationsOf (token: string) =
            if token.Length = 1 then [||]
            else Text.RegularExpressions.Regex.Matches(token, @"\[[^\]]*\]") |> Seq.map (fun m -> m.Value) |> Array.ofSeq

        /// Names the peptide group of a sequence: the smaller of the sequence without modifications and its pseudo-reversal.
        /// A target, its decoy and all their modified forms share the group.
        let canonicalSequence (sequence: string) =
            let residues = sequenceTokens sequence |> List.map residueOf
            let toString (residues: char list) = residues |> Array.ofList |> System.String
            let own = toString residues
            let reversed = toString (Fragmentation.Series.pseudoReverse residues)
            if String.CompareOrdinal(own, reversed) <= 0 then own else reversed

        /// The modifications of the sequences, numbered in order.
        let modificationIndex (sequences: seq<string>) =
            sequences
            |> Seq.collect (fun sequence -> sequenceTokens sequence |> Seq.collect modificationsOf)
            |> Seq.distinct
            |> Seq.sort
            |> Seq.mapi (fun i modification -> modification, i)
            |> dict

        let featureWidth (modIndex: Collections.Generic.IDictionary<string,int>) = 21 + modIndex.Count + 2 * nPairs + 20

        /// Amino acid counts, length, count per modification, the unordered residue pairs at the positions 1 and n-1 and at 2 and n-2,
        /// and the C-terminal residue. A sequence and its pseudo-reversal give the same vector.
        let featureVector (modIndex: Collections.Generic.IDictionary<string,int>) (sequence: string) =
            let tokens = sequenceTokens sequence
            let features = Array.zeroCreate<float32> (featureWidth modIndex)
            let residues = tokens |> List.map residueOf |> Array.ofList
            let n = residues.Length
            let index residue = match residueIndex.TryGetValue residue with | true, i -> i | _ -> -1
            for token in tokens do
                let i = index (residueOf token)
                if i >= 0 then features.[i] <- features.[i] + 1.f
                for modification in modificationsOf token do
                    match modIndex.TryGetValue modification with
                    | true, j -> features.[21 + j] <- features.[21 + j] + 1.f
                    | _ -> ()
            features.[20] <- float32 n
            let pairs1 = 21 + modIndex.Count
            let pairs2 = pairs1 + nPairs
            let cTerminal = pairs2 + nPairs
            let setPair offset a b =
                if a >= 0 && b >= 0 then features.[offset + pairIndex a b] <- 1.f
            if n >= 3 then
                setPair pairs1 (index residues.[0]) (index residues.[n - 2])
                setPair pairs2 (index residues.[1]) (index residues.[n - 3])
            if n >= 1 && index residues.[n - 1] >= 0 then features.[cTerminal + index residues.[n - 1]] <- 1.f
            features

        [<CLIMutable>]
        type Example = { Features: float32 []; Label: float32 }

        [<CLIMutable>]
        type Scored = { Score: float32 }
        type Predictor =
            {
                /// 1.4826 times the median absolute deviation of the out-of-fold residuals of the training rows
                Sigma          : float
                R2             : float
                MedianAbsError : float
                /// Predicts sequences with the predictor of the fold of their peptide group, which was trained without that fold
                Predict        : string [] -> float32 []
            }

        /// Fold of a peptide group: a seeded FNV-1a hash of the canonical sequence, finished with the murmur3 mixer. It is stable between runs,
        /// and a group needs no table.
        let foldOfGroup (seed: int) (group: string) =
            let mutable hash = 2166136261u ^^^ uint32 seed
            for c in group do
                hash <- (hash ^^^ uint32 c) * 16777619u
            hash <- hash ^^^ (hash >>> 16)
            hash <- hash * 0x85ebca6bu
            hash <- hash ^^^ (hash >>> 13)
            hash <- hash * 0xc2b2ae35u
            hash <- hash ^^^ (hash >>> 16)
            int (hash % uint32 nFolds)

        let private median (values: float []) = FSharp.Stats.Array.median values

        /// Fits one FastTree regression per fold on the rows (sequence, scan time) of the other folds. The folds are groups of peptides,
        /// so every sequence is predicted by a predictor that has not seen its group.
        let fit (seed: int) (rows: (string * float) []) =
            // an own context keeps the random stream of the classifier untouched
            let ctx = new MLContext(seed)
            let modIndex = modificationIndex (rows |> Seq.map fst)
            let schema =
                let schema = SchemaDefinition.Create(typeof<Example>)
                schema.["Features"].ColumnType <- VectorDataViewType(NumberDataViewType.Single, featureWidth modIndex)
                schema
            let load (examples: seq<Example>) = ctx.Data.LoadFromEnumerable(examples, schema)
            let foldOf = canonicalSequence >> foldOfGroup seed
            let rowFolds = rows |> Array.map (fst >> foldOf)
            let trainer =
                ctx.Regression.Trainers.FastTree(
                    Trainers.FastTree.FastTreeRegressionTrainer.Options(
                        LabelColumnName = "Label",
                        FeatureColumnName = "Features",
                        NumberOfTrees = 200,
                        NumberOfLeaves = 64,
                        MinimumExampleCountPerLeaf = 20,
                        LearningRate = 0.1,
                        NumberOfThreads = 1
                    )
                )
            let models =
                Array.init nFolds (fun k ->
                    rows
                    |> Seq.indexed
                    |> Seq.filter (fun (i, _) -> rowFolds.[i] <> k)
                    |> Seq.map (fun (_, (sequence, time)) -> { Features = featureVector modIndex sequence; Label = float32 time })
                    |> load
                    |> trainer.Fit
                )
            let predict (sequences: string []) =
                let folds = sequences |> Array.map foldOf
                let predictions = Array.zeroCreate<float32> sequences.Length
                for k in 0 .. nFolds - 1 do
                    let indices = [| 0 .. sequences.Length - 1 |] |> Array.filter (fun i -> folds.[i] = k)
                    if indices.Length > 0 then
                        let scores =
                            indices
                            |> Seq.map (fun i -> { Features = featureVector modIndex sequences.[i]; Label = 0.f })
                            |> load
                            |> models.[k].Transform
                            |> fun scored -> ctx.Data.CreateEnumerable<Scored>(scored, reuseRowObject = false)
                            |> Seq.map (fun x -> x.Score)
                            |> Array.ofSeq
                        indices |> Array.iteri (fun j i -> predictions.[i] <- scores.[j])
                predictions
            let outOfFold = predict (rows |> Array.map fst)
            let times = rows |> Array.map snd
            let residuals = Array.map2 (fun time (prediction: float32) -> time - float prediction) times outOfFold
            let medianResidual = median residuals
            let mean = Array.average times
            {
                Sigma          = 1.4826 * median (residuals |> Array.map (fun r -> abs (r - medianResidual))) + 1e-9
                R2             = 1. - (residuals |> Array.sumBy (fun r -> r * r)) / (times |> Array.sumBy (fun t -> (t - mean) * (t - mean)))
                MedianAbsError = median (residuals |> Array.map abs)
                Predict        = predict
            }

        /// The deviation of the scan time of a PSM from the time predicted for its peptide, in robust standard deviations, or None if the
        /// file has too few confident peptides or the predictor is too inaccurate. The predictor is trained on the targets with q < 0.01
        /// among the best candidates per scan by the starting score, one row per modified sequence and global modification with the median scan time.
        /// Each ModSequenceID is predicted once.
        let tryInitDeviation (logger: NLog.Logger) (seed: int) (seedScore: Dto.PeptideSpectrumMatchingResult -> float32) (psms: Dto.PeptideSpectrumMatchingResult []) =
            let watch = Diagnostics.Stopwatch.StartNew()
            let bestPSMPerScan =
                psms
                |> Array.groupBy (fun x -> x.ScanNr)
                |> Array.map (fun (_, scan) -> scan |> Array.maxBy seedScore)
            let q = BioFSharp.Mz.FDRControl.calculateQValueStorey bestPSMPerScan (fun x -> x.Label <> 1) (fun x -> float (seedScore x)) (fun x -> float (seedScore x))
            let confident = bestPSMPerScan |> Array.filter (fun x -> x.Label = 1 && q (float (seedScore x)) < 0.01)
            let rows =
                confident
                |> Array.groupBy (fun x -> x.ModSequenceID, x.GlobalMod)
                |> Array.sortBy fst
                |> Array.map (fun (_, group) -> group.[0].StringSequence, group |> Array.map (fun x -> x.ScanTime) |> median)
            let skip reason =
                logger.Warn (sprintf "Retention time feature not used (%s): %i confident PSMs, %i unique peptides, %.1f s" reason confident.Length rows.Length watch.Elapsed.TotalSeconds)
                None
            if rows.Length < minPeptides then
                skip (sprintf "fewer than %i confident unique peptides" minPeptides)
            else
                let predictor = fit seed rows
                // written as a negated comparison so that a missing R2 skips as well
                if not (predictor.R2 >= minR2) then
                    skip (sprintf "out-of-fold R2 %.3f below %.1f, median absolute error %.2f min, sigma %.2f min" predictor.R2 minR2 predictor.MedianAbsError predictor.Sigma)
                else
                    let predictions = Collections.Generic.Dictionary<int, float32>()
                    // The features are the same for a sequence and its pseudo-reversal, so a decoy is predicted from the sequence of its target.
                    psms
                    |> Array.distinctBy (fun x -> x.ModSequenceID)
                    |> Array.chunkBySize 50000
                    |> Array.iter (fun chunk ->
                        let predicted = predictor.Predict (chunk |> Array.map (fun x -> x.StringSequence))
                        Array.iter2 (fun (x: Dto.PeptideSpectrumMatchingResult) p -> predictions.[x.ModSequenceID] <- p) chunk predicted
                        )
                    logger.Trace (sprintf "Retention time feature used: %i confident PSMs, %i unique peptides, out-of-fold R2 %.3f, median absolute error %.2f min, sigma %.2f min, %.1f s"
                                    confident.Length rows.Length predictor.R2 predictor.MedianAbsError predictor.Sigma watch.Elapsed.TotalSeconds)
                    Some (fun (x: Dto.PeptideSpectrumMatchingResult) -> float32 ((x.ScanTime - float predictions.[x.ModSequenceID]) / predictor.Sigma))

    let restorePSMID psmID =
        FSharpAux.String.split '_' psmID  
        |> Array.item 0
        |> FSharpAux.String.split '-'
        |> String.concat " "

    ///
    let initToPSMToLearn maxCharge fastaHeaderToName proteinAndClvIdxLookUp (rtDeviation: Dto.PeptideSpectrumMatchingResult -> float32) (psm:Dto.PeptideSpectrumMatchingResult) =
        let lookUp = proteinAndClvIdxLookUp psm.PepSequenceID
        let missCleavages,nTerminal,cTerminal =
            match lookUp with
            | [] -> 0,"",""
            | items ->
                items
                |> List.minBy (fun ((protID,name,seq),clvgS,clvgE,missClv) -> missClv)
                |> fun ((protID,name,seq),clvgS,clvgE,missClv) ->
                    let nTerm,cTerm = getFlankingAminoAcids clvgS clvgE seq
                    missClv,nTerm + ".", "." + cTerm
        let proteinNames =
            let parseLabel label = if label = 1 then "" else "decoy_"
            match lookUp with
            | [] -> sprintf "%snoProt" (parseLabel psm.Label)
            | items ->
                items
                |> List.map (fun ((protID,name,seq),clvgS,clvgE,missClv) -> sprintf "%s%s" (parseLabel psm.Label) (fastaHeaderToName name))
                |> String.concat ";"
        
        let flankedPepSequence =
            let sequence = if psm.Label = 1 then psm.StringSequence else decoySequence psm.StringSequence
            nTerminal + sequence + cTerminal
        let deviation = rtDeviation psm
        {
            PSMId                        = psm.PSMId
            Label                        = if psm.Label = 1 then true else false
            ScanNr                       = psm.ScanNr
            Charge                       = psm.Charge 
            PrecursorMZ                  = float32 psm.PrecursorMZ
            TheoMass                     = float32 psm.TheoMass
            AbsDeltaMass                 = float32 psm.AbsDeltaMass
            PeptideLength                = float32 psm.PeptideLength
            MissCleavages                = float32 missCleavages
            SequestScore                 = float32 psm.SequestScore
            SequestNormDeltaBestToRest   = float32 psm.SequestNormDeltaBestToRest
            SequestNormDeltaNext         = float32 psm.SequestNormDeltaNext
            AndroScore                   = float32 psm.AndroScore
            AndroNormDeltaBestToRest     = float32 psm.AndroNormDeltaBestToRest
            AndroNormDeltaNext           = float32 psm.AndroNormDeltaNext
            XtandemScore                 = float32 psm.XtandemScore
            XtandemNormDeltaBestToRest   = float32 psm.XtandemNormDeltaBestToRest
            XtandemNormDeltaNext         = float32 psm.XtandemNormDeltaNext
            Hyperscore                   = if Double.IsFinite psm.Hyperscore then float32 psm.Hyperscore else 0.f
            LogExpectscore               = logExpectscoreOf psm
            RtDeviation                  = deviation
            AbsRtDeviation               = abs deviation
            Peptide                      = flankedPepSequence
            Protein                      = proteinNames
        }
             
    /// Writes the results as a tab separated table, or only the header if there are none.
    let writeResults (outFilePath: string) (result: Dto.PSMStatisticsResult []) =
        if result.Length = 0 then
            // Seq.CSV takes the header from the first row and fails on an empty sequence.
            Reflection.FSharpType.GetRecordFields(typeof<Dto.PSMStatisticsResult>)
            |> Seq.map (fun field -> field.Name)
            |> String.concat "\t"
            |> FSharpAux.IO.FileIO.writeStringToFile false outFilePath
        else
            result
            |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
            |> FSharpAux.IO.FileIO.writeToFile false outFilePath

    /// The distinct (sequence, global modification, charge) ions of a file of results.
    let readIons (qpsmPath: string) =
        FSharpAux.IO.SchemaReader.Csv.CsvReader<Dto.PSMStatisticsResult>(SchemaMode=FSharpAux.IO.SchemaReader.Csv.Fill).ReadFile(qpsmPath,'\t',false,1)
        |> Seq.map (fun x -> x.StringSequence, x.GlobalMod, x.Charge)
        |> Seq.distinct
        |> Array.ofSeq

    ///
    let psmStats diagCharts (processParams:PSMStatisticsParams) (outputDir:string) (cn:SQLiteConnection) (psms:string) =

        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension psms)

        logger.Trace (sprintf "Input file: %s" psms)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)

        logger.Trace (sprintf "Now performing peptide spectrum matching: %s Results will be written to: %s" psms outputDir)
        let outFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension psms) + ".qpsm"
            Path.Combine [|outputDir;fileName|]
        let plotDirectory =
            let fileName = sprintf "%s_plots" (Path.GetFileNameWithoutExtension psms)
            let path = Path.Combine [|outputDir;fileName|]
            if System.IO.Directory.Exists path then
                path
            else
                System.IO.Directory.CreateDirectory path |> ignore
                path
        logger.Trace (sprintf "outFilePath:%s" outFilePath)
        logger.Trace (sprintf "plotDirectory:%s" plotDirectory)
        logger.Trace "Copy peptide DB into Memory"
        let memoryDB = SearchDB.copyDBIntoMemory cn 
        let pepDBTr = memoryDB.BeginTransaction()
        logger.Trace "Copy peptide DB into Memory: finished"
        
        logger.Trace "Read scored PSMs."
        let psms =
            FSharpAux.IO.SchemaReader.Csv.CsvReader<Dto.PeptideSpectrumMatchingResult>().ReadFile(psms,'\t',false,1)
            |> Array.ofSeq
        logger.Trace "Read scored PSMs: finished"

        logger.Trace "Prepare processing functions."
        // PeptideSpectrumMatchingTIMs fills Hyperscore and Expectscore, PeptideSpectrumMatching writes NaN.
        // With finite values the two become features and the expectation value seeds the training.
        let useTimsScores = psms |> Array.forall (fun x -> Double.IsFinite x.Hyperscore && Double.IsFinite x.Expectscore)
        logger.Trace (sprintf "Hyperscore and expectation value present: %b" useTimsScores)
        let seedScore (x: PSMToLearn) = seedScoreOf useTimsScores x.LogExpectscore x.SequestScore
        let maxCharge = psms |> Array.map (fun x -> x.Charge) |> Array.max
        let proteinAndClvIdxLookUp = initProteinAndClvIdxLookUp memoryDB pepDBTr
        let toPSMToLearn = initToPSMToLearn maxCharge processParams.FastaHeaderToName proteinAndClvIdxLookUp
        logger.Trace "Finished preparing processing functions."
        
        match processParams.Threshold with 
        | Threshold.Estimate estParams -> 
            try
            // The predictor is fitted once per file, before any model is trained.
            let rtDeviation =
                if estParams.RetentionTimeFeature then
                    RetentionTime.tryInitDeviation logger 1024 (fun x -> seedScoreOf useTimsScores (logExpectscoreOf x) (float32 x.SequestScore)) psms
                else
                    None
            let useRtFeature = rtDeviation.IsSome
            let toPSMToLearn = toPSMToLearn (defaultArg rtDeviation (fun _ -> 0.f))
            logger.Trace "Prepare training pipeline"
            let ctx = new ML.MLContext(1024)
            // prefix: "" or "Fold<k>_" with cross-validation, added to log lines and chart names
            let trainModel prefix positives' negatives' =
                let data = ctx.Data.LoadFromEnumerable(positives' + negatives')
                let split = ctx.Data.TrainTestSplit(data, testFraction= 0.1)
                let featureColumns =
                    [|
                        "OneHotCharge"
                        "PrecursorMZ"
                        "TheoMass"
                        "AbsDeltaMass"
                        "PeptideLength"
                        "MissCleavages"
                        "SequestScore"
                        "SequestNormDeltaBestToRest"
                        "SequestNormDeltaNext"
                        "AndroScore"
                        "AndroNormDeltaBestToRest"
                        "AndroNormDeltaNext"
                        "XtandemScore"
                        "XtandemNormDeltaBestToRest"
                        "XtandemNormDeltaNext"
                        if useTimsScores then
                            "Hyperscore"
                            "LogExpectscore"
                        if useRtFeature then
                            "RtDeviation"
                            "AbsRtDeviation"
                    |]
                // A leaf output of up to 100 (the default) lets a few leaves push near-empty spectra and decoys to the top.
                let fastTreeOptions =
                    Trainers.FastTree.FastTreeBinaryTrainer.Options(
                        FeatureColumnName = "featuresNorm",
                        LabelColumnName = "Label",
                        MaximumTreeOutput = 10.,
                        NumberOfThreads = 1
                    )
                let pipeline =
                    (ctx.Transforms.Categorical.OneHotEncoding("OneHotCharge","Charge") |> downcastPipeline)
                        .Append(ctx.Transforms.Concatenate("Features", featureColumns))
                        .Append(ctx.Transforms.NormalizeMeanVariance("featuresNorm","Features"))
                        .Append(ctx.BinaryClassification.Trainers.FastTree(fastTreeOptions))
             
                let model = pipeline.Fit(split.TrainSet)    
                try
                    let metrics = ctx.BinaryClassification.Evaluate(model.Transform(split.TestSet),labelColumnName="Label")
                    logger.Trace (sprintf "%sMetrics on the test split: accuracy %.4f, AUC %.4f, F1 %.4f" prefix metrics.Accuracy metrics.AreaUnderRocCurve metrics.F1Score)
                    if diagCharts then 
                        Chart.Column(
                            [
                            "Accuracy",metrics.Accuracy
                            "PositivePrecision",metrics.PositivePrecision
                            "PositiveRecall",metrics.PositiveRecall
                            "NegativePrecision",metrics.NegativePrecision
                            "NegativeRecall",metrics.NegativeRecall
                            "F1Score",metrics.F1Score
                            ]
                            )
                        |> Chart.withTitle "Metrics"
                        |> Chart.saveHtml (Path.Combine [|plotDirectory;prefix + "Metrics"|])
                with
                | ex -> logger.Warn (sprintf "%sMetrics on the test split not available: %s" prefix ex.Message)
                let predF = ctx.Model.CreatePredictionEngine<PSMToLearn,PSMPrediction>(model)
                let cache = System.Collections.Generic.Dictionary<PSMToLearn,PSMPrediction>(HashIdentity.Reference)
                let predict (psm: PSMToLearn) =
                    match cache.TryGetValue psm with
                    | true, prediction -> prediction
                    | _ ->
                        let prediction = predF.Predict psm
                        cache.[psm] <- prediction
                        prediction
                predict
            logger.Trace "Prepare training pipeline:finished"
            logger.Trace "Converting psms to PSMToLearn format."
            let psmsToLearn =
                psms
                |> Array.map toPSMToLearn
            logger.Trace "Converting psms to PSMToLearn format: finished"        

            let applyModel prefix iteration (trainedModel: PSMToLearn -> PSMPrediction) psms =
                let bestPSMPerScan = 
                    psms
                    |> Array.groupBy (fun x -> x.ScanNr)
                    |> Array.map (fun (psmId,psms) -> 
                        psms 
                        |> Array.maxBy (fun x -> 
                            (trainedModel x).Score
                            )
                        )
                let getQ = BioFSharp.Mz.FDRControl.calculateQValueStorey bestPSMPerScan (fun x -> x.Label |> not) (fun x -> float (trainedModel x).Score ) (fun x -> float (trainedModel x).Score)             
                let getPep = 
                    let bw = 
                        bestPSMPerScan
                        |> Array.map (fun x -> (trainedModel x).Score |> float)
                        |> FSharp.Stats.Distributions.Bandwidth.nrd0
                        |> fun x -> x / 4.
                    match estParams.PepValueFittingMethod with
                    //| LinearLogit -> ProteomIQon.FDRControl'.initCalculateLin logger bw (fun x -> x.Label |> not) (fun x -> float (trainedModel x).Score ) (fun x -> float (trainedModel x).Score) bestPSMPerScan
                    | IRLS -> ProteomIQon.PepValueCalculation.initCalculatePEPValueIRLS logger 1. (fun x -> x.Label |> not) (fun x -> float (trainedModel x).Score ) (fun x -> float (trainedModel x).Score) bestPSMPerScan
                let scoreVsQ = 
                    bestPSMPerScan
                    |> Array.map (fun x -> (trainedModel x).Score, getQ (float (trainedModel x).Score))
                let scoreVSPep = 
                    bestPSMPerScan
                    |> Array.map (fun x -> (trainedModel x).Score, getPep (float (trainedModel x).Score))
                let tar = 
                    bestPSMPerScan 
                    |> Array.filter (fun x -> x.Label) 
                    |> Array.map (fun x -> (trainedModel x).Score )
                let decoy = 
                    bestPSMPerScan 
                    |> Array.filter (fun x -> x.Label |> not) 
                    |> Array.map (fun x -> (trainedModel x).Score )
                // With cross-validation, models are compared by their targets with q <= 0.01, independent of the output thresholds.
                let nPosTar = 
                    if estParams.CrossValidation then
                        bestPSMPerScan |> Array.filter (fun x -> x.Label && getQ (float (trainedModel x).Score) <= 0.01) |> Array.length
                    else
                        let posQ = scoreVsQ |> Array.filter (fun x -> snd x <= estParams.QValueThreshold) |> Array.length
                        let posPep = scoreVSPep |> Array.filter (fun x -> snd x <= estParams.PepValueThreshold) |> Array.length
                        System.Math.Min(posQ,posPep)
                if diagCharts then 
                    [
                        [
                        Chart.Histogram(tar)
                        |> Chart.withTraceInfo "positives"

                        Chart.Histogram(decoy)
                        |> Chart.withTraceInfo "negatives"
                        ]
                        |> Chart.combine
                        |> Chart.withAxisAnchor(Y=1)
                        Chart.Line(scoreVsQ |> Array.sortBy fst)
                        |> Chart.withAxisAnchor(Y=2)
                        |> Chart.withTraceInfo "QValue"
                        Chart.Line(scoreVSPep |> Array.sortBy fst)
                        |> Chart.withAxisAnchor(Y=2)
                        |> Chart.withTraceInfo "PepValue"
                    ]
                    |> Chart.combine
                    |> Chart.withXAxisStyle("Score")
                    |> Chart.withYAxisStyle("Count",Side=StyleParam.Side.Left,Id=StyleParam.SubPlotId.YAxis 1,ShowGrid=false)
                    |> Chart.withYAxisStyle("FDR",Side=StyleParam.Side.Right,Id=StyleParam.SubPlotId.YAxis 2,Overlaying=StyleParam.LinearAxisId.Y 1,ShowGrid=false,MinMax=(0.,0.5))
                    |> Chart.withTitle (sprintf "#iteration: %i, %i with q < 0.01" iteration (scoreVsQ |> Array.filter (fun x -> snd x <= 0.01) |> Array.length))
                    |> Chart.saveHtml (Path.Combine [|plotDirectory;prefix + "separationAtIteration_" + iteration.ToString()|])
                {
                    NPositivesAtFDR = nPosTar
                    CalcQValue     = getQ
                    CalcPepValue   = getPep
                    Model          = trainedModel
                }

            /// Trains the model on psmsToTrain. The positives are seeded with seedScore and then taken from the scores of the previous model.
            let train prefix (psmsToTrain: PSMToLearn []) =
                logger.Trace (sprintf "%sTraining on %i PSMs of %i scans" prefix psmsToTrain.Length (psmsToTrain |> Array.distinctBy (fun x -> x.ScanNr) |> Array.length))
                let bestPSMPerScan = 
                    psmsToTrain
                    |> Array.groupBy (fun x -> x.ScanNr)
                    |> Array.map (fun (psmId,psms) -> 
                        psms |> Array.maxBy seedScore
                        )
                let q = BioFSharp.Mz.FDRControl.calculateQValueStorey bestPSMPerScan (fun s -> s.Label |> not) (fun s -> float (seedScore s)) (fun s -> float (seedScore s)) 

                let scoreVsQ = 
                    bestPSMPerScan
                    |> Array.map (fun x -> seedScore x,q (float (seedScore x)))

                let tar = 
                    bestPSMPerScan 
                    |> Array.filter (fun x -> x.Label) 
                    |> Array.map seedScore

                let decoy = 
                    bestPSMPerScan 
                    |> Array.filter (fun x -> x.Label |> not) 
                    |> Array.map seedScore
                if diagCharts then
                    [
                        [
                        Chart.Histogram(tar)
                        |> Chart.withTraceInfo "positives"

                        Chart.Histogram(decoy)
                        |> Chart.withTraceInfo "negatives"
                        ]
                        |> Chart.combine
                        |> Chart.withAxisAnchor(Y=1)
                        Chart.Point(scoreVsQ)
                        |> Chart.withAxisAnchor(Y=2)
                    ]
                    |> Chart.combine
                    |> Chart.withXAxisStyle("Score")
                    |> Chart.withYAxisStyle("Count",Side=StyleParam.Side.Left,Id=StyleParam.SubPlotId.YAxis 1,ShowGrid=false)
                    |> Chart.withYAxisStyle("FDR",Side=StyleParam.Side.Right,Id=StyleParam.SubPlotId.YAxis 2,Overlaying=StyleParam.LinearAxisId.Y 1,ShowGrid=false,MinMax=(0.,0.5))
                    |> Chart.withTitle (sprintf "#%i with q < 0.01" (scoreVsQ |> Array.filter (fun x -> snd x <= 0.01) |> Array.length))
                    |> Chart.saveHtml (Path.Combine [|plotDirectory;prefix + "InitialSeparation"|])
                logger.Trace (prefix + "Selecting positives for training")
                let positives' = 
                    bestPSMPerScan 
                    |> Array.filter (fun x -> q (float (seedScore x)) < 0.001)
                    |> Array.filter (fun x -> x.Label = true)
                    |> Array.map (fun x -> x.ScanNr,x)
                    |> Map.ofArray
                logger.Trace (prefix + "Selecting negatives for training")
                let negatives' =
                    psmsToTrain
                    |> Array.filter (fun x -> positives' |> Map.containsKey x.ScanNr |> not)
                    |> Array.filter (fun x -> x.Label = false)
                    |> Set.ofArray
        
                let predict = trainModel prefix (positives'|> Map.toArray |> Array.map snd |> Set.ofArray) negatives' 
                let initModel = applyModel prefix 0 predict psmsToTrain
                let rec loop (prevModel:AppliedModel) i =
                    if i >= estParams.MaxIterations then 
                        prevModel
                    else
                        let bestPSMPerScan = 
                            psmsToTrain
                            |> Array.groupBy (fun x -> x.ScanNr)
                            |> Array.map (fun (psmId,psms) -> 
                                psms 
                                |> Array.maxBy (fun x -> 
                                    (prevModel.Model x).Score
                                    )
                                )   
                        logger.Trace (prefix + "Selecting positives for training")
                        // With cross-validation, the positives are the targets with q < 0.01, independent of the output thresholds.
                        let positives' = 
                            bestPSMPerScan 
                            |> Array.filter (fun x -> 
                                if estParams.CrossValidation then
                                    (prevModel.Model x).Score |> float |> prevModel.CalcQValue < 0.01
                                else
                                    let passQVal = (prevModel.Model x).Score |> float |> prevModel.CalcQValue  < (estParams.QValueThreshold / 10.)
                                    let passPEPVal = (prevModel.Model x).Score |> float |> prevModel.CalcPepValue  < (estParams.PepValueThreshold)
                                    passQVal && passPEPVal
                                )
                            |> Array.filter (fun x -> x.Label = true)
                            |> Array.map (fun x -> x.ScanNr,x)
                            |> Map.ofArray
                        logger.Trace (prefix + "Selecting negatives for training")
                        let negatives' =
                            psmsToTrain
                            |> Array.filter (fun x -> positives' |> Map.containsKey x.ScanNr |> not)
                            |> Array.filter (fun x -> x.Label = false)
                            |> Set.ofArray
                        logger.Trace (sprintf "%sTraining iteration #%i with %i positives and %i negatives" prefix i positives'.Count negatives'.Count )
                        try
                            let predict = trainModel prefix (positives'|> Map.toArray |> Array.map snd |> Set.ofArray)  negatives' 
                            let newModel = applyModel prefix i predict psmsToTrain                    
                            let ratio = (prevModel.NPositivesAtFDR |> float) /  (newModel.NPositivesAtFDR |> float) 
                            if ratio > 1. then 
                                prevModel 
                            elif 1. - ratio < estParams.MinimumIncreaseBetweenIterations then    
                                newModel
                            else
                                loop newModel (i+1)
                        with
                        | ex -> logger.Trace (sprintf "%sFall back to previous model: %s" prefix ex.Message); prevModel
                let refinedModel = loop initModel 0
                logger.Trace (sprintf "%sFinal model: %i %s in-sample" prefix refinedModel.NPositivesAtFDR (if estParams.CrossValidation then "targets with q <= 0.01" else "scans pass the q-value and PEP thresholds"))
                refinedModel

            // finalScore: the model score, with cross-validation the normalised score of the model that did not see the scan
            let finalScore, calcQValue, calcPepValue =
                if not estParams.CrossValidation then
                    let refinedModel = train "" psmsToLearn
                    (fun x -> (refinedModel.Model x).Score |> float), refinedModel.CalcQValue, refinedModel.CalcPepValue
                else
                    // Percolator-style: each scan is scored by the model trained on the other two of three scan folds.
                    let nFolds = 3
                    let foldOfScan =
                        let rnd = System.Random(1024)
                        psmsToLearn
                        |> Array.map (fun x -> x.ScanNr)
                        |> Array.distinct
                        |> Array.sort
                        |> Array.sortBy (fun _ -> rnd.Next())
                        |> Array.mapi (fun i scanNr -> scanNr, i % nFolds)
                        |> dict
                    let scores = System.Collections.Generic.Dictionary<PSMToLearn,float>(HashIdentity.Reference)
                    // 1 / (score at q = 0.01 - median decoy score) of each fold: one model score unit in normalised units
                    let unitPerFold =
                        Array.init nFolds (fun k ->
                            let prefix = sprintf "Fold%i_" (k + 1)
                            let trainSet, heldOut = psmsToLearn |> Array.partition (fun x -> foldOfScan.[x.ScanNr] <> k)
                            let model = train prefix trainSet
                            let score x = (model.Model x).Score |> float
                            // Normalise as Percolator does: the score at q = 0.01 maps to 0 and the median decoy score to -1, so the folds are comparable.
                            let bestPSMPerScan =
                                heldOut
                                |> Array.groupBy (fun x -> x.ScanNr)
                                |> Array.map (fun (scanNr,psms) -> psms |> Array.maxBy score)
                            let getQ = BioFSharp.Mz.FDRControl.calculateQValueStorey bestPSMPerScan (fun x -> x.Label |> not) score score
                            let passing = bestPSMPerScan |> Array.filter (fun x -> x.Label && getQ (score x) <= 0.01) |> Array.map score
                            let scoreAtQ = if passing.Length > 0 then Array.min passing else bestPSMPerScan |> Array.map score |> Array.max
                            let medianDecoy = bestPSMPerScan |> Array.filter (fun x -> x.Label |> not) |> Array.map score |> FSharp.Stats.Array.median
                            logger.Trace (sprintf "%sHeld-out: %i PSMs of %i scans, %i targets with q <= 0.01, score at q = 0.01 %f, median decoy score %f" prefix heldOut.Length bestPSMPerScan.Length passing.Length scoreAtQ medianDecoy)
                            for x in heldOut do
                                scores.[x] <- (score x - scoreAtQ) / (scoreAtQ - medianDecoy)
                            1. / (scoreAtQ - medianDecoy)
                        )
                    let score x = scores.[x]
                    let bestPSMPerScan =
                        psmsToLearn
                        |> Array.groupBy (fun x -> x.ScanNr)
                        |> Array.map (fun (scanNr,psms) -> psms |> Array.maxBy score)
                    let getQ = BioFSharp.Mz.FDRControl.calculateQValueStorey bestPSMPerScan (fun x -> x.Label |> not) score score
                    let getPep =
                        match estParams.PepValueFittingMethod with
                        // the histogram bin width of 1 model score unit, as without cross-validation
                        | IRLS -> ProteomIQon.PepValueCalculation.initCalculatePEPValueIRLS logger (Array.average unitPerFold) (fun x -> x.Label |> not) score score bestPSMPerScan
                    logger.Trace (sprintf "Cross-validated scores: %i targets with q <= 0.01" (bestPSMPerScan |> Array.filter (fun x -> x.Label && getQ (score x) <= 0.01) |> Array.length))
                    score, getQ, getPep

            let bestPSMPerScan = 
                psmsToLearn
                |> Array.groupBy (fun x -> x.ScanNr)
                |> Array.map (fun (psmId,psms) -> 
                    psms 
                    |> Array.maxBy (fun x -> 
                        finalScore x
                        )
                    )

            let qpsm = 
                bestPSMPerScan 
                |> Array.filter (fun x -> 
                    let qValPass = finalScore x |> calcQValue  < estParams.QValueThreshold
                    let pepValPass = finalScore x |> calcPepValue < estParams.PepValueThreshold
                    qValPass && pepValPass
                    )
                |> Array.filter (fun x -> x.Label = true)
                |> Array.map (fun x -> x.PSMId,x)
                |> Map.ofArray

            let result: Dto.PSMStatisticsResult [] =
                psms
                |> Array.choose (fun candidatePSM ->
                                    match candidatePSM.Label with
                                    | x when x = 1 ->
                                        match Map.tryFind candidatePSM.PSMId qpsm with
                                        | Some validPSM ->
                                            let psmID = restorePSMID validPSM.PSMId 
                                            let score = finalScore validPSM 
                                            let qValue = calcQValue score 
                                            let pepValue = calcPepValue score 
                                            Some {
                                            PSMId                       = psmID
                                            GlobalMod                   = candidatePSM.GlobalMod
                                            PepSequenceID               = candidatePSM.PepSequenceID
                                            ModSequenceID               = candidatePSM.ModSequenceID
                                            Label                       = candidatePSM.Label
                                            ScanNr                      = candidatePSM.ScanNr
                                            ScanTime                    = candidatePSM.ScanTime
                                            Charge                      = candidatePSM.Charge
                                            PrecursorMZ                 = candidatePSM.PrecursorMZ
                                            TheoMass                    = candidatePSM.TheoMass
                                            AbsDeltaMass                = candidatePSM.AbsDeltaMass
                                            PeptideLength               = candidatePSM.PeptideLength
                                            MissCleavages               = int validPSM.MissCleavages
                                            SequestScore                = candidatePSM.SequestScore
                                            SequestNormDeltaBestToRest  = candidatePSM.SequestNormDeltaBestToRest
                                            SequestNormDeltaNext        = candidatePSM.SequestNormDeltaNext
                                            AndroScore                  = candidatePSM.AndroScore
                                            AndroNormDeltaBestToRest    = candidatePSM.AndroNormDeltaBestToRest
                                            AndroNormDeltaNext          = candidatePSM.AndroNormDeltaNext
                                            XtandemScore                = candidatePSM.XtandemScore
                                            XtandemNormDeltaBestToRest  = candidatePSM.XtandemNormDeltaBestToRest
                                            XtandemNormDeltaNext        = candidatePSM.XtandemNormDeltaNext
                                            ModelScore                  = score
                                            QValue                      = qValue
                                            PEPValue                    = pepValue
                                            StringSequence              = candidatePSM.StringSequence
                                            ProteinNames                = validPSM.Protein
                                            }
                                        | None  -> 
                                            None
                                    | _ -> 
                                        None
                                )
            logger.Trace (sprintf "Number of results: %i" result.Length)
            writeResults outFilePath result
            logger.Trace "Done."
            with 
            | ex ->
                logger.Trace (sprintf "%A" ex)
                logger.Trace "Machine learning based psm filtering did not succeed, please try a different PEP-value fitting method. This method was developed for large data sets so it might fail for data sets too small of a size."
                writeResults outFilePath [||]
                logger.Trace "Done."
        | Threshold.Fixed cutoff -> 
            let pepSequenceIDToMissCleavagesAndProt =
                psms
                |> Array.map (toPSMToLearn (fun _ -> 0.f))
                |> Array.filter (fun candidatePSM -> candidatePSM.Label )
                |> Array.map (fun candidatePSM -> candidatePSM.PSMId,(candidatePSM.Protein,candidatePSM.MissCleavages))
                |> Map.ofArray
            let result : Dto.PSMStatisticsResult [] =
                psms
                |> Array.filter (fun x -> x.SequestScore > cutoff.SequestLike && x.AndroScore > cutoff.Andromeda)
                |> Array.groupBy (fun x -> x.ScanNr)
                |> Array.map (fun (scanNr,scans) ->
                    scans
                    |> Array.maxBy (fun x -> x.SequestScore)
                    )
                |> Array.choose (fun candidatePSM  ->
                                    match candidatePSM.Label with
                                    | x when x = 1 ->
                                        let psmID = restorePSMID candidatePSM.PSMId 
                                        let proteins,missCleavages = pepSequenceIDToMissCleavagesAndProt.[candidatePSM.PSMId]
                                        let res : Dto.PSMStatisticsResult = 
                                            {
                                            PSMId                       = psmID
                                            GlobalMod                   = candidatePSM.GlobalMod
                                            PepSequenceID               = candidatePSM.PepSequenceID
                                            ModSequenceID               = candidatePSM.ModSequenceID
                                            Label                       = candidatePSM.Label
                                            ScanNr                      = candidatePSM.ScanNr
                                            ScanTime                    = candidatePSM.ScanTime
                                            Charge                      = candidatePSM.Charge
                                            PrecursorMZ                 = candidatePSM.PrecursorMZ
                                            TheoMass                    = candidatePSM.TheoMass
                                            AbsDeltaMass                = candidatePSM.AbsDeltaMass
                                            PeptideLength               = candidatePSM.PeptideLength
                                            MissCleavages               = int missCleavages
                                            SequestScore                = candidatePSM.SequestScore
                                            SequestNormDeltaBestToRest  = candidatePSM.SequestNormDeltaBestToRest
                                            SequestNormDeltaNext        = candidatePSM.SequestNormDeltaNext
                                            AndroScore                  = candidatePSM.AndroScore
                                            AndroNormDeltaBestToRest    = candidatePSM.AndroNormDeltaBestToRest
                                            AndroNormDeltaNext          = candidatePSM.AndroNormDeltaNext
                                            XtandemScore                = candidatePSM.XtandemScore
                                            XtandemNormDeltaBestToRest  = candidatePSM.XtandemNormDeltaBestToRest
                                            XtandemNormDeltaNext        = candidatePSM.XtandemNormDeltaNext
                                            ModelScore                  = candidatePSM.SequestScore
                                            QValue                      = nan
                                            PEPValue                    = nan
                                            StringSequence              = candidatePSM.StringSequence
                                            ProteinNames                = proteins 
                                            }
                                        res
                                        |> Some
                                    | _ -> 
                                        None
                                )
            logger.Trace (sprintf "Number of results: %i" result.Length)
            writeResults outFilePath result

