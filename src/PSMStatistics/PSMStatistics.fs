namespace ProteomIQon

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
                    if psm.Label = 1 then
                        missClv,nTerm + ".", "." + cTerm
                    else
                        missClv, cTerm + ".", "." + nTerm
        let proteinNames =
            let parseLabel label = if label = 1 then "" else "decoy_"
            match lookUp with
            | [] -> sprintf "%snoProt" (parseLabel psm.Label)
            | items ->
                items
                |> List.map (fun ((protID,name,seq),clvgS,clvgE,missClv) -> sprintf "%s%s" (parseLabel psm.Label) (fastaHeaderToName name))
                |> String.concat ";"
        
        let flankedPepSequence =
            let sequence = if psm.Label = 1 then psm.StringSequence else FSharpAux.String.rev psm.StringSequence
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
        let maxCharge = psms |> Array.map (fun x -> x.Charge) |> Array.max
        let proteinAndClvIdxLookUp = initProteinAndClvIdxLookUp memoryDB pepDBTr
        let toPSMToLearn = initToPSMToLearn maxCharge processParams.FastaHeaderToName proteinAndClvIdxLookUp
        logger.Trace "Finished preparing processing functions."
        
        match processParams.Threshold with 
        | Threshold.Estimate estParams -> 
            try
            logger.Trace "Prepare training pipeline"
            let ctx = new ML.MLContext(1024)
            let trainModel positives' negatives' =
                let data = ctx.Data.LoadFromEnumerable(positives' + negatives')
                let split = ctx.Data.TrainTestSplit(data, testFraction= 0.1)
                let pipeline =    
                    (ctx.Transforms.Categorical.OneHotEncoding("OneHotCharge","Charge") |> downcastPipeline)
                        .Append(
                            ctx.Transforms.Concatenate(
                                "Features",
                                "OneHotCharge",
                                "PrecursorMZ",
                                "TheoMass",
                                "AbsDeltaMass",
                                "PeptideLength",
                                "MissCleavages",
                                "SequestScore",
                                "SequestNormDeltaBestToRest",
                                "SequestNormDeltaNext",
                                "AndroScore",
                                "AndroNormDeltaBestToRest",
                                "AndroNormDeltaNext",
                                "XtandemScore",
                                "XtandemNormDeltaBestToRest",
                                "XtandemNormDeltaNext"
                                )
                            )
                        .Append(ctx.Transforms.NormalizeMeanVariance("featuresNorm","Features"))
                        .Append(ctx.BinaryClassification.Trainers.FastTree(featureColumnName="featuresNorm",labelColumnName="Label"))
             
                let model = pipeline.Fit(split.TrainSet)    
                let metrics = ctx.BinaryClassification.Evaluate(model.Transform(split.TestSet),labelColumnName="Label")
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
                    |> Chart.SaveHtmlAs (Path.Combine [|plotDirectory;"Metrics"|])
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

            let bestPSMPerScan = 
                psmsToLearn
                |> Array.groupBy (fun x -> x.ScanNr)
                |> Array.map (fun (psmId,psms) -> 
                    psms |> Array.maxBy (fun x -> x.SequestScore)
                    )
            let q = BioFSharp.Mz.FDRControl.calculateQValueStorey bestPSMPerScan (fun s -> s.Label |> not) (fun s -> float s.SequestScore) (fun s -> float s.SequestScore) 

            let scoreVsQ = 
                bestPSMPerScan
                |> Array.map (fun x -> x.SequestScore,q (float x.SequestScore))

            let tar = 
                bestPSMPerScan 
                |> Array.filter (fun x -> x.Label) 
                |> Array.map (fun x -> x.SequestScore )

            let decoy = 
                bestPSMPerScan 
                |> Array.filter (fun x -> x.Label |> not) 
                |> Array.map (fun x -> x.SequestScore)
            if diagCharts then
                [
                    [
                    Chart.Histogram(tar)
                    |> Chart.withTraceName "positives"

                    Chart.Histogram(decoy)
                    |> Chart.withTraceName "negatives"
                    ]
                    |> Chart.Combine
                    |> Chart.withAxisAnchor(Y=1)
                    Chart.Point(scoreVsQ)
                    |> Chart.withAxisAnchor(Y=2)
                ]
                |> Chart.Combine
                |> Chart.withX_AxisStyle("Score")
                |> Chart.withY_AxisStyle("Count",Side=StyleParam.Side.Left,Id=1,Showgrid=false)
                |> Chart.withY_AxisStyle("FDR",Side=StyleParam.Side.Right,Id=2,Overlaying=StyleParam.AxisAnchorId.Y 1,Showgrid=false,MinMax=(0.,0.5))
                |> Chart.withTitle (sprintf "#%i with q < 0.01" (scoreVsQ |> Array.filter (fun x -> snd x <= 0.01) |> Array.length))
                |> Chart.SaveHtmlAs (Path.Combine [|plotDirectory;"InitialSeparation"|])
            logger.Trace "Selecting positives for training"
            let positives' = 
                bestPSMPerScan 
                |> Array.filter (fun x -> q (float x.SequestScore) < 0.001)
                |> Array.filter (fun x -> x.Label = true)
                |> Array.map (fun x -> x.ScanNr,x)
                |> Map.ofArray
            logger.Trace "Selecting negatives for training"
            let negatives' =
                psmsToLearn
                |> Array.filter (fun x -> positives' |> Map.containsKey x.ScanNr |> not)
                |> Array.filter (fun x -> x.Label = false)
                |> Set.ofArray
        
            let predict = trainModel (positives'|> Map.toArray |> Array.map snd |> Set.ofArray) negatives' 
        
            let applyModel iteration (trainedModel: PSMToLearn -> PSMPrediction) psms =
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
                    | LinearLogit -> ProteomIQon.FDRControl'.initCalculateLin logger bw (fun x -> x.Label |> not) (fun x -> float (trainedModel x).Score ) (fun x -> float (trainedModel x).Score) bestPSMPerScan
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
                        |> Chart.withTraceName "positives"

                        Chart.Histogram(decoy)
                        |> Chart.withTraceName "negatives"
                        ]
                        |> Chart.Combine
                        |> Chart.withAxisAnchor(Y=1)
                        Chart.Point(scoreVsQ)
                        |> Chart.withAxisAnchor(Y=2)
                        |> Chart.withTraceName "QValue"
                        Chart.Point(scoreVSPep)
                        |> Chart.withAxisAnchor(Y=2)
                        |> Chart.withTraceName "PepValue"
                    ]
                    |> Chart.Combine
                    |> Chart.withX_AxisStyle("Score")
                    |> Chart.withY_AxisStyle("Count",Side=StyleParam.Side.Left,Id=1,Showgrid=false)
                    |> Chart.withY_AxisStyle("FDR",Side=StyleParam.Side.Right,Id=2,Overlaying=StyleParam.AxisAnchorId.Y 1,Showgrid=false,MinMax=(0.,0.5))
                    |> Chart.withTitle (sprintf "#iteration: %i, %i with q < 0.01" iteration (scoreVsQ |> Array.filter (fun x -> snd x <= 0.01) |> Array.length))
                    |> Chart.SaveHtmlAs (Path.Combine [|plotDirectory;"separationAtIteration_" + iteration.ToString()|])
                {
                    NPositivesAtFDR = nPosTar
                    CalcQValue     = getQ
                    CalcPepValue   = getPep
                    Model          = trainedModel
                }
        
            let initModel = applyModel 0 predict psmsToLearn
            let refinedModel = 
                let rec loop (prevModel:AppliedModel) i =
                    if i >= estParams.MaxIterations then 
                        prevModel
                    else
                        let bestPSMPerScan = 
                            psmsToLearn
                            |> Array.groupBy (fun x -> x.ScanNr)
                            |> Array.map (fun (psmId,psms) -> 
                                psms 
                                |> Array.maxBy (fun x -> 
                                    (prevModel.Model x).Score
                                    )
                                )   
                        logger.Trace "Selecting positives for training"
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
                        logger.Trace "Selecting negatives for training"
                        let negatives' =
                            psmsToLearn
                            |> Array.filter (fun x -> positives' |> Map.containsKey x.ScanNr |> not)
                            |> Array.filter (fun x -> x.Label = false)
                            |> Set.ofArray
                        logger.Trace (sprintf "Training iteration #%i with %i positives and %i negatives" i positives'.Count negatives'.Count )
                        try
                            let predict = trainModel (positives'|> Map.toArray |> Array.map snd |> Set.ofArray)  negatives' 
                            let newModel = applyModel i predict psmsToLearn                    
                            let ratio = (prevModel.NPositivesAtFDR |> float) /  (newModel.NPositivesAtFDR |> float) 
                            if ratio > 1. then 
                                prevModel 
                            elif 1. - ratio < estParams.MinimumIncreaseBetweenIterations then    
                                newModel
                            else
                                loop newModel (i+1)
                        with
                        | _ -> logger.Trace "Fall back to previous model"; prevModel
                loop initModel 0

            let bestPSMPerScan = 
                psmsToLearn
                |> Array.groupBy (fun x -> x.ScanNr)
                |> Array.map (fun (psmId,psms) -> 
                    psms 
                    |> Array.maxBy (fun x -> 
                        (refinedModel.Model x).Score
                        )
                    )

            let qpsm = 
                bestPSMPerScan 
                |> Array.filter (fun x -> 
                    let qValPass = (refinedModel.Model x).Score |> float |> refinedModel.CalcQValue  < estParams.QValueThreshold
                    let pepValPass = (refinedModel.Model x).Score |> float |> refinedModel.CalcPepValue < estParams.PepValueThreshold
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
                                            let score = (refinedModel.Model validPSM).Score |> float 
                                            let qValue = refinedModel.CalcQValue score 
                                            let pepValue = refinedModel.CalcPepValue score 
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

