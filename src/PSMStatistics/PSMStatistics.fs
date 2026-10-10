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

    /// The decoy of a peptide sequence as the searches build it with
    /// Fragmentation.Series.pseudoReverse: the C-terminal residue stays in place and the rest is
    /// reversed, or the whole sequence if that gives the sequence itself. The terminator "*" and the
    /// gap "-" are dropped first, as the searches drop them. A residue keeps the
    /// modification tokens written before it, such as "[ox]M".
    let decoySequence (sequence: string) =
        Text.RegularExpressions.Regex.Matches(sequence, @"(\[[^\]]*\])*[^\[\]]")
        |> Seq.map (fun m -> m.Value)
        |> List.ofSeq
        |> List.filter (fun token -> token <> "*" && token <> "-")
        |> Fragmentation.Series.pseudoReverse
        |> String.concat ""

    let restorePSMID psmID =
        FSharpAux.String.split '_' psmID  
        |> Array.item 0
        |> FSharpAux.String.split '-'
        |> String.concat " "

    ///
    let initToPSMToLearn maxCharge fastaHeaderToName proteinAndClvIdxLookUp (psm:Dto.PeptideSpectrumMatchingResult) =
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
            LogExpectscore               = if Double.IsFinite psm.Expectscore && psm.Expectscore > 0. then float32 (-log10 psm.Expectscore) else 0.f
            Peptide                      = flankedPepSequence
            Protein                      = proteinNames
        }
             
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
        let seedScore (x: PSMToLearn) = if useTimsScores then x.LogExpectscore else x.SequestScore
        let maxCharge = psms |> Array.map (fun x -> x.Charge) |> Array.max
        let proteinAndClvIdxLookUp = initProteinAndClvIdxLookUp memoryDB pepDBTr
        let toPSMToLearn = initToPSMToLearn maxCharge processParams.FastaHeaderToName proteinAndClvIdxLookUp
        logger.Trace "Finished preparing processing functions."
        
        match processParams.Threshold with 
        | Threshold.Estimate estParams -> 
            try
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
                    |]
                let pipeline =    
                    (ctx.Transforms.Categorical.OneHotEncoding("OneHotCharge","Charge") |> downcastPipeline)
                        .Append(ctx.Transforms.Concatenate("Features", featureColumns))
                        .Append(ctx.Transforms.NormalizeMeanVariance("featuresNorm","Features"))
                        .Append(ctx.BinaryClassification.Trainers.FastTree(featureColumnName="featuresNorm",labelColumnName="Label"))
             
                let model = pipeline.Fit(split.TrainSet)    
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
                let predF = ctx.Model.CreatePredictionEngine<PSMToLearn,PSMPrediction>(model)
                let predict psm = 
                    psm 
                    |> predF.Predict
                predict 
            logger.Trace "Prepare training pipeline:finished"
            logger.Trace "Converting psms to PSMToLearn format."
            let psmsToLearn =
                psms
                |> Array.map toPSMToLearn
                |> Array.filter (fun x -> nan.Equals(x.SequestNormDeltaNext) = false && nan.Equals(x.AndroNormDeltaNext) = false )
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
                let nPosTar = 
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
                        let positives' = 
                            bestPSMPerScan 
                            |> Array.filter (fun x -> 
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
                        | _ -> logger.Trace (prefix + "Fall back to previous model"); prevModel
                let refinedModel = loop initModel 0
                logger.Trace (sprintf "%sFinal model: %i scans pass the q-value and PEP thresholds in-sample" prefix refinedModel.NPositivesAtFDR)
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
            result
            |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
            |> FSharpAux.IO.FileIO.writeToFile false outFilePath
            logger.Trace "Done."
            with 
            | ex ->
                logger.Trace (sprintf "%A" ex)
                logger.Trace "Machine learning based psm filtering did not succeed, please try a different PEP-value fitting method. This method was developed for large data sets so it might fail for data sets too small of a size."
                let fields = Reflection.FSharpType.GetRecordFields(typeof<Dto.PSMStatisticsResult>)
                fields
                |> Seq.map(fun field -> field.Name)
                |> String.concat "\t"
                |> FSharpAux.IO.FileIO.writeStringToFile false outFilePath
                logger.Trace "Done."
        | Threshold.Fixed cutoff -> 
            let pepSequenceIDToMissCleavagesAndProt =
                psms
                |> Array.map toPSMToLearn
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
                                            ModelScore                  = nan
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
            result
            |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
            |> FSharpAux.IO.FileIO.writeToFile false outFilePath

