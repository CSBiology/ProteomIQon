namespace ProteomIQon


open System.IO
open Argu
open System.Data.SQLite
open System
open ProteomIQon.Core
open Core.MzIO
open Dto
open FSharp.Stats
open BioFSharp.Mz
open FSharpAux.IO.SchemaReader
open FSharp.Plotly
open BioFSharp
open MzIO.Processing
open Microsoft
open Microsoft.ML
open Microsoft.ML.Transforms
open Microsoft.ML.Data   
open FSharpAux.IO.SchemaReader.Attribute
open Microsoft.ML.Transforms.Text

module QuantBasedAlignment = 

    ///
    type AlignmentResult = 
        {
            [<FieldAttribute(0)>]
            GlobalMod                    : int
            [<FieldAttribute(1)>]
            PepSequenceID                : int
            [<FieldAttribute(2)>]
            ModSequenceID                : int
            [<FieldAttribute(3)>]
            ScanTime                     : float
            [<FieldAttribute(4)>]
            Charge                       : int
            [<FieldAttribute(5)>]
            Mz                           : float
            [<FieldAttribute(10)>]
            StringSequence               : string
            [<FieldAttribute(11)>]
            ProteinNames                 : string
        }   


    type PeptideIon = 
        {
            Sequence             : string
            GlobalMod            : int
            Charge               : int        
        }
    
    let toPeptideIon (qp:QuantificationResult) = 
        {
            Sequence             = qp.StringSequence
            GlobalMod            = qp.GlobalMod
            Charge               = qp.Charge
        }

    type AlignmentFile = {
        FileName                   : string
        QuantifiedPeptides         : Map<PeptideIon,QuantificationResult> 
        MissingPeptides            : PeptideIon []
        GainedPeptides             : AlignmentResult []
        }

    type AlignmentParams = {
        Placeholder : bool 
        }

    [<CLIMutable>]
    type PeptideForLearning = 
        {
            [<ColumnName("Sequence")>]
            Sequence : string
            [<ColumnName("SourceScanTime")>]
            SourceScanTime : float32
            [<ColumnName("TargetScanTime")>]
            TargetScanTime : float32
        }

    [<CLIMutable>]
    type ScanTimePrediction = 
        {
            [<ColumnName("Score")>]
            TargetScanTime : float32
        }

    /// Retrieves the scan time based on the fitted parameter values (HULQ output).
    let getScanTime (qp:QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.Params_Light.Split(';').[1] |> float
        else
            qp.Params_Heavy.Split(';').[1] |> float
        with
        | _ -> nan

    /// Retrieves the scan time based on the fitted parameter values (HULQ output).
    let getStabw (qp:QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.Params_Light.Split(';').[2] |> float
        else
            qp.Params_Heavy.Split(';').[2] |> float
        with
        | _ -> nan

    /// Retrieves the scan time based on the fitted parameter values (HULQ output).
    let tryGetScanTime (qp:QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.Params_Light.Split(';').[1] 
            |> float
            |> Some
        else
            qp.Params_Heavy.Split(';').[1] 
            |> float
            |> Some
        with
        | _ -> None

    ///
    let getQuantifiedPeptides (quantFilePath:string) = 
        ///
        let peptides =
            Csv.CsvReader<QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(quantFilePath,'\t',false,1)
            |> Array.ofSeq
        let filteredPeptides =
            peptides
            |> Array.filter (fun qp -> 
                // Filter for peptides where fit assures a good estimation of the ScanTime
                (qp.GlobalMod = 0 && qp.Params_Light <> "") || (qp.GlobalMod = 1 && qp.Params_Heavy <> "")
                )
            |> Array.filter (fun qp -> tryGetScanTime qp |> Option.isSome)
            |> Array.filter (fun qp ->
                (qp.GlobalMod = 0 && qp.Difference_SearchRT_FittedRT_Light |> abs < 0.5) || (qp.GlobalMod = 1 && qp.Difference_SearchRT_FittedRT_Heavy |> abs < 0.5) 
                )
        filteredPeptides 
    
    /// Reads quant files, filters for quality quantifications and creates AlignmentFiles.
    let createAlignmentFiles (quantFilePaths:string []) = 
        let quantifiedPeptides = 
            quantFilePaths
            |> Array.map getQuantifiedPeptides
        let allPeptideIons =
            quantifiedPeptides
            |> Array.concat
            |> Array.map (fun qp -> toPeptideIon qp)
        Array.map2 (fun (qFp:string) qps -> 
                let presentPeptides = 
                    qps 
                    |> Array.map (fun qp -> toPeptideIon qp, qp)
                    |> Array.distinctBy fst
                    |> Map.ofArray
                let missingPeptides = 
                    allPeptideIons
                    |> Array.filter (fun pep -> 
                        let pepUnMod = {pep with GlobalMod = 1}
                        let pepMod = {pep with GlobalMod = 1}
                        presentPeptides.ContainsKey pepMod |> not || presentPeptides.ContainsKey pepUnMod |> not 
                        )
                {
                    FileName                   = Path.GetFileNameWithoutExtension(qFp)
                    QuantifiedPeptides         = presentPeptides 
                    MissingPeptides            = missingPeptides 
                    GainedPeptides             = [||]
                }
            ) quantFilePaths quantifiedPeptides

    /// Calculates the median of absolute scan time differences between shared peptide Ions. This serves as an estimator for overall
    /// File difference
    let calculateFileDifference (a:AlignmentFile) (b:AlignmentFile) =
        a.QuantifiedPeptides
        |> Seq.choose (fun a -> 
            match Map.tryFind a.Key b.QuantifiedPeptides with 
            | Some qp -> abs (getScanTime qp - getScanTime a.Value) |> Some
            | None -> None 
            )
        |> Seq.filter (isNan >> not)
        |> Seq.median

    /// determines order of alignment based on the calculated file difference. 
    // calculateFileDifference can be easily factored out if needed. 
    let findAlignmentOrder (alignmentFiles: AlignmentFile []) =
        alignmentFiles
        |> Array.mapi (fun i af -> 
                af,
                [| 
                    for j = 0 to alignmentFiles.Length-1 do 
                        if j <> i then calculateFileDifference af alignmentFiles.[j], alignmentFiles.[j]
                |]
                |> Array.sortBy fst
            )
    
    let createAlignmentResult (quantifiedPeptide:QuantificationResult) (scanTimePrediction:ScanTimePrediction) = 
        {
            GlobalMod                    = quantifiedPeptide.GlobalMod
            PepSequenceID                = quantifiedPeptide.PepSequenceID
            ModSequenceID                = quantifiedPeptide.ModSequenceID
            ScanTime                     = float scanTimePrediction.TargetScanTime
            Charge                       = quantifiedPeptide.Charge
            Mz                           = quantifiedPeptide.PrecursorMZ
            StringSequence               = quantifiedPeptide.StringSequence
            ProteinNames                 = quantifiedPeptide.ProteinNames
        }   

    let downcastPipeline (x : IEstimator<_>) = 
        match x with 
        | :? IEstimator<ITransformer> as y -> y
        | _ -> failwith "downcastPipeline: expecting a IEstimator<ITransformer>"

    type ModelMetrics = 
        {
        Metrics     : RegressionMetrics
        X_Test      : float []
        Y_Test      : float []
        YHat_Test   : float []
        }

    let formatString s = String.filter (fun x -> Char.IsUpper x) s
    
    ///
    let initAlign (ctx:MLContext) (pepsForLearning: PeptideForLearning []) = 
        let data = ctx.Data.LoadFromEnumerable(pepsForLearning)
        let split = ctx.Data.TrainTestSplit(data, testFraction= 0.1)
        let pipeline =    
            (ctx.Transforms.Concatenate("Features","SourceScanTime") |> downcastPipeline)
              //.Append(ctx.Regression.Trainers.FastTree(featureColumnName="Features",labelColumnName="TargetScanTime",numberOfLeaves=200))
              .Append(ctx.Regression.Trainers.FastTree(featureColumnName="Features",labelColumnName="TargetScanTime"))
              
        //let pipeline =
        //    (ctx.Transforms.Text.TokenizeIntoCharactersAsKeys("MessageChars", "Sequence") |> downcastPipeline)
        //        .Append(ctx.Transforms.Text.ProduceNgrams("BagOfOnechar", "MessageChars", ngramLength= 1))
        //        .Append(ctx.Transforms.Concatenate("Features","SourceScanTime"(*,"BagOfOnechar"*)))
        //        .Append(ctx.Regression.Trainers.FastForest (featureColumnName="Features",labelColumnName="TargetScanTime",numberOfLeaves=200))
        ///minmax
        //let pipeline =    
        //    (ctx.Transforms.Concatenate("ToMinMaxNormalized","SourceScanTime") |> downcastPipeline)
        //        .Append(ctx.Transforms.NormalizeMinMax("Features", "ToMinMaxNormalized"))
        //        .Append(ctx.Transforms.Concatenate("ToMinMaxNormalizedTarget","TargetScanTime"))
        //        .Append(ctx.Transforms.NormalizeMinMax("MinMaxTargetScanTime", "ToMinMaxNormalizedTarget"))

        //        //.Append(ctx.Regression.Trainers.FastTree(featureColumnName="Features",labelColumnName="TargetScanTime",numberOfLeaves=200))
        //      .Append(ctx.Regression.Trainers.FastTree(featureColumnName="Features",labelColumnName="MinMaxTargetScanTime"))
        let model = pipeline.Fit(split.TrainSet)    
        let metrics = 
            
            let metrics = ctx.Regression.Evaluate(model.Transform(split.TestSet),labelColumnName="TargetScanTime")
            //let metrics = ctx.Regression.Evaluate(model.Transform(split.TestSet),labelColumnName="MinMaxTargetScanTime")
            let evalView = model.Transform(split.TestSet)
            let x    = evalView.GetColumn<float32>("SourceScanTime") |> Seq.map float |> Array.ofSeq
            let y    = evalView.GetColumn<float32>("TargetScanTime") |> Seq.map float |> Array.ofSeq
            let yHat = evalView.GetColumn<float32>("Score") |> Seq.map float |> Array.ofSeq
            //let x    = evalView.GetColumn<float32 []>("Features") |> Seq.map  (Array.item 0) |> Seq.map float |> Array.ofSeq
            //let y    = evalView.GetColumn<float32 []>("MinMaxTargetScanTime") |> Seq.map  (Array.item 0) |> Seq.map float |> Array.ofSeq
            //let yHat = evalView.GetColumn<float32>("Score") |> Seq.map float |> Array.ofSeq
            {
                Metrics     = metrics
                X_Test      = x
                Y_Test      = y
                YHat_Test   = yHat
            }

        let predF = ctx.Model.CreatePredictionEngine<PeptideForLearning,ScanTimePrediction>(model)
        let toPeptideForLearning (quantifiedPeptide:QuantificationResult) = {Sequence=formatString quantifiedPeptide.StringSequence; SourceScanTime=getScanTime quantifiedPeptide |> float32;TargetScanTime=nan|> float32}
        let predict quantifiedPeptide = 
            predF.Predict(toPeptideForLearning quantifiedPeptide)
            |> createAlignmentResult quantifiedPeptide
        metrics, predict
            
        //let initPrediction pepsForLearning:= 
    ///
    let performAlignment align (target: AlignmentFile) (source: AlignmentFile) =
        
        let peptidesForLearning = 
            target.QuantifiedPeptides
            |> Seq.choose (fun tarQP -> 
                    match Map.tryFind tarQP.Key source.QuantifiedPeptides with 
                    | Some sourceQP ->
                        {Sequence = formatString sourceQP.StringSequence;SourceScanTime=getScanTime sourceQP|> float32;TargetScanTime=getScanTime tarQP.Value|> float32}
                        |> Some
                    | None -> None
                )
            |> Array.ofSeq
            |> Array.shuffleFisherYates
        
        let metrics,model: ModelMetrics*(QuantificationResult->AlignmentResult) = 
            align peptidesForLearning
        
        let peptideIonsToTransfer,peptideIonsStillMissing = 
            target.MissingPeptides
            |> Array.fold (fun (pepsToTransfer,stillMissingPeps) missingPep -> 
                match Map.tryFind missingPep source.QuantifiedPeptides with 
                | Some pepToTransfer -> (pepToTransfer::pepsToTransfer,stillMissingPeps)
                | None -> (pepsToTransfer,missingPep::stillMissingPeps)
                ) ([],[])

        let stabwMedian = 
            target.QuantifiedPeptides
            |> Seq.map (fun x -> getStabw x.Value)
            |> Array.ofSeq
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Seq.median

        let alignmentResults : AlignmentResult [] = 
            peptideIonsToTransfer 
            |> List.toArray 
            |> Array.map model 
        let yVsYHat = 
            Chart.Point(metrics.Y_Test,metrics.YHat_Test)
            |> Chart.withTraceName (sprintf "#TestPeptides:%i Rsquared:%f RMS:%f, %s " metrics.X_Test.Length metrics.Metrics.RSquared metrics.Metrics.MeanAbsoluteError source.FileName)
            |> Chart.withX_AxisStyle("target ScanTimes (Y_Test)")
            |> Chart.withY_AxisStyle("predicted target ScanTimes (YHat_Test)")
        let xVsDifferenceYandYHat = 
            Chart.Point(metrics.X_Test, Array.map2 (fun y yHat -> y - yHat) metrics.Y_Test metrics.YHat_Test)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("target ScanTimes (Y_Test) - predicted target ScanTimes (YHat_Test)")
            |> Chart.withTraceName (sprintf "#TestPeptides:%i Rsquared:%f RMS:%f, %s " metrics.X_Test.Length metrics.Metrics.RSquared metrics.Metrics.MeanAbsoluteError source.FileName)
        let xVsDifferenceYandYHatNormed = 
            Chart.Point(metrics.X_Test, Array.map2 (fun y yHat -> (y - yHat) / stabwMedian) metrics.Y_Test metrics.YHat_Test)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("target ScanTimes (Y_Test) - predicted target ScanTimes (YHat_Test)")
            |> Chart.withTraceName (sprintf "#TestPeptides:%i Rsquared:%f RMS:%f, %s " metrics.X_Test.Length metrics.Metrics.RSquared metrics.Metrics.MeanAbsoluteError source.FileName)


        [yVsYHat;xVsDifferenceYandYHat;xVsDifferenceYandYHatNormed]
        |> Chart.Stack 2
        |> Chart.withSize(2000.,1500.),
        {target with MissingPeptides = peptideIonsStillMissing |> Array.ofList; GainedPeptides = Array.append target.GainedPeptides alignmentResults}
        

    
    let alignFiles (logger:NLog.Logger) (processParams:AlignmentParams) (outputDir:string) (quantFiles:string) = 

        //let logger = Logging.createLogger (Path.GetFileNameWithoutExtension quantFiles)

        logger.Trace (sprintf "Input directory: %s" quantFiles)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        
        let getPlotFilePathFilePath (plotName:string) (fileName:string) =
            let fileName = (Path.GetFileNameWithoutExtension fileName) + "_" + plotName 
            Path.Combine [|outputDir;fileName|]
        
        logger.Trace "Init align function"
        let ctx = new ML.MLContext()
        let align = initAlign ctx
        logger.Trace "Init align function: finished"
         
        logger.Trace "Reading and preparing .quant file for alignment"
        let quantFiles = 
            System.IO.Directory.GetFiles (quantFiles, "*.quant")
        let alignmentFiles = createAlignmentFiles quantFiles
        logger.Trace "Reading and preparing .quant file for alignment: finished"
               
        logger.Trace "Determining Alignment order"
        let alignmentFilesOrdered = findAlignmentOrder alignmentFiles
        logger.Trace "Determining Alignment order: finished"
        
        
        logger.Trace "Plotting file distances"
        let chart = 
            alignmentFilesOrdered
            |> Array.map (fun (target,sources) ->
                    Chart.Point(sources |> Array.mapi (fun i x -> (snd x).FileName, fst x))
                    |> Chart.withTraceName target.FileName
                    |> Chart.withX_AxisStyle("FileNames")
                    |> Chart.withY_AxisStyle("Median absolute difference of peptide ion scan times")
                    |> Chart.withSize(1000.,1000.)
                    |> Chart.SaveHtmlAs(getPlotFilePathFilePath "differences" target.FileName)
                )
        logger.Trace "Plotting file distances: finished"
        
            
        logger.Trace "Performing Alignment"
        let alignments = 
            alignmentFilesOrdered 
            |> Array.map (fun (target,sources) -> 
                    let chart,result =
                        sources
                        |> Array.fold (fun (charts,tar) (dis,source)  -> 
                            logger.Trace (sprintf "Performing Alignment %s vs %s" tar.FileName source.FileName)
                            let chart',tar' = 
                                performAlignment align tar source
                            logger.Trace (sprintf "Performing Alignment %s vs %s: finished" tar.FileName source.FileName)
                            logger.Trace (sprintf "Missing peptides before: %i, now:%i" tar.MissingPeptides.Length tar'.MissingPeptides.Length)
                            chart'::charts,tar'
                            ) ([],target)
                    chart
                    |> Chart.Combine
                    |> Chart.withTitle(target.FileName)
                    |> Chart.SaveHtmlAs(getPlotFilePathFilePath "whatever" target.FileName)                    
                    result 
                )
                        
        logger.Trace "Performing Alignment: finished"
                
        alignments
        |> ignore
        

        




































