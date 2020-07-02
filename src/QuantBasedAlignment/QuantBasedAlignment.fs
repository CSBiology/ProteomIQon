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
open FSharpAux.Colors.Table.StatisticalGraphics24

module QuantBasedAlignment = 

    let paletteArray =
        [|
            Blue1     
            Blue2     
            Blue3     
            Red1      
            Red2      
            Red3      
            LightBlue1
            LightBlue2
            LightBlue3
            LightRed1 
            LightRed2 
            LightRed3 
            Green1    
            Green2    
            Green3    
            Orange1   
            Orange2   
            Orange3   
            Cyan1     
            Cyan2     
            Cyan3     
            Magenta1  
            Magenta2  
            Magenta3  
        |]

    /// Define the random number generator outside of a potential loop.
    let getRandomColor (rnd: System.Random) =
        let index = rnd.Next(0,23)
        paletteArray.[index]

    /// Retrieves the scan time based on the fitted parameter values (HULQ output).
    let getTargetIntensity (qp:QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.Quant_Light
        else
            qp.Quant_Heavy
        with
        | _ -> nan

    /// Retrieves the scan time based on the fitted parameter values (HULQ output).
    let getTargetScanTime (qp:QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.Params_Light.[1] 
        else
            qp.Params_Heavy.[1] 
        with
        | _ -> nan

    /// Retrieves the scan time based on the fitted parameter values (HULQ output).
    let getTargetStabw (qp:QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.Params_Light.[2] 
        else
            qp.Params_Heavy.[2] 
        with
        | _ -> nan

    /// Retrieves the scan time based on the fitted parameter values (HULQ output).
    let getTargetScanTimeDifference (qp:QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.Difference_SearchRT_FittedRT_Light
        else
            qp.Difference_SearchRT_FittedRT_Heavy
        with
        | _ -> nan

    /// Retrieves the scan time based on the fitted parameter values (HULQ output).
    let tryTargetGetScanTime (qp:QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.Params_Light.[1] 
            |> float
            |> Some
        else
            qp.Params_Heavy.[1] 
            |> float
            |> Some
        with
        | _ -> None

    /// 
    let getTargetRtTrace (qp:QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.RtTrace_Light
        else
            qp.RtTrace_Heavy 
        with
        | _ -> [||]

    /// 
    let getTargetIntensityTrace (qp:QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.IntensityTrace_Corrected_Light
        else
            qp.IntensityTrace_Corrected_Heavy
        with
        | _ -> [||]

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
            Charge                       : int
            [<FieldAttribute(4)>]
            Mz                           : float
            [<FieldAttribute(5)>]
            StringSequence               : string
            [<FieldAttribute(10)>]
            ProteinNames                 : string
            [<FieldAttribute(11)>]
            PredictedScanTime            : float
            [<FieldAttribute(12)>][<TraceConverter>]
            RtTrace_SourceFile           : float []
            [<FieldAttribute(13)>][<TraceConverter>]
            IntensityTrace_SourceFile    : float []
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
            Sequence                     : string
            [<ColumnName("SourceScanTime")>]
            SourceScanTime               : float32
            [<ColumnName("SourceIntensity")>]
            SourceIntensity              : float32
            [<ColumnName("SourceStabw")>]
            SourceStabw                  : float32
            [<ColumnName("TargetScanTime")>]
            TargetScanTime               : float32
            [<ColumnName("RtTrace_SourceFile")>]
            RtTrace_SourceFile           : float []
            [<ColumnName("IntensityTrace_SourceFile")>]
            IntensityTrace_SourceFile    : float []
            [<ColumnName("RtTrace_TargetFile")>]
            RtTrace_TargetFile           : float []
            [<ColumnName("IntensityTrace_TargetFile")>]
            IntensityTrace_TargetFile    : float []
        }
    
    let formatString s = String.filter (fun x -> Char.IsUpper x) s
    
    let toPeptideForLearning (targetPep:QuantificationResult option) (sourcePep:QuantificationResult) = 
        match targetPep with 
        | Some tP -> 
            {
                Sequence                    = formatString sourcePep.StringSequence 
                SourceScanTime              = getTargetScanTime sourcePep  |> float32
                SourceIntensity             = getTargetIntensity sourcePep |> float32
                SourceStabw                 = getTargetStabw sourcePep     |> float32
                TargetScanTime              = getTargetScanTime tP         |> float32
                RtTrace_SourceFile          = getTargetRtTrace sourcePep
                IntensityTrace_SourceFile   = getTargetIntensityTrace sourcePep
                RtTrace_TargetFile          = getTargetRtTrace tP
                IntensityTrace_TargetFile   = getTargetIntensityTrace tP
            }
        | None -> 
            {
                Sequence                    = formatString sourcePep.StringSequence 
                SourceScanTime              = getTargetScanTime sourcePep  |> float32
                SourceIntensity             = getTargetIntensity sourcePep |> float32
                SourceStabw                 = getTargetStabw sourcePep     |> float32
                TargetScanTime              = nan                          |> float32
                RtTrace_SourceFile          = getTargetRtTrace sourcePep
                IntensityTrace_SourceFile   = getTargetIntensityTrace sourcePep
                RtTrace_TargetFile          = [||]
                IntensityTrace_TargetFile   = [||]
            }

    
    [<CLIMutable>]
    type ScanTimePrediction = 
        {
            [<ColumnName("Score")>]
            TargetScanTime : float32
        }


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
                (qp.GlobalMod = 0 && qp.Params_Light |> Array.isEmpty |> not) || (qp.GlobalMod = 1 && qp.Params_Heavy |> Array.isEmpty |> not)
                )
            |> Array.filter (fun qp -> tryTargetGetScanTime qp |> Option.isSome)
            |> Array.filter (fun qp -> (getTargetScanTimeDifference qp |> abs) / (getTargetStabw qp) < 2. )          
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
            |> Array.distinct
        Array.map2 (fun (qFp:string) qps -> 
                let presentPeptides = 
                    qps 
                    |> Array.map (fun qp -> toPeptideIon qp, qp)
                    |> Map.ofArray
                let missingPeptides = 
                    allPeptideIons
                    |> Array.filter (fun pep -> 
                        let pepUnMod = {pep with GlobalMod = 0}
                        let pepMod = {pep with GlobalMod = 1}
                        (presentPeptides.ContainsKey pepMod || presentPeptides.ContainsKey pepUnMod) 
                        |> not
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
            | Some qp -> abs (getTargetScanTime qp - getTargetScanTime a.Value) |> Some
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
            Charge                       = quantifiedPeptide.Charge
            Mz                           = Mass.toMZ (quantifiedPeptide.TheoMass) (float quantifiedPeptide.Charge) 
            StringSequence               = quantifiedPeptide.StringSequence
            ProteinNames                 = quantifiedPeptide.ProteinNames
            PredictedScanTime            = float scanTimePrediction.TargetScanTime
            RtTrace_SourceFile           = getTargetRtTrace quantifiedPeptide
            IntensityTrace_SourceFile    = getTargetIntensityTrace quantifiedPeptide
        }   

    let downcastPipeline (x : IEstimator<_>) = 
        match x with 
        | :? IEstimator<ITransformer> as y -> y
        | _ -> failwith "downcastPipeline: expecting a IEstimator<ITransformer>"

    type ModelMetrics = 
        {
        Metrics             : RegressionMetrics
        X_Intensities       : float []
        X_Stabw             : float []        
        X_Test              : float []
        Y_Test              : float []
        YHat_Test           : float []
        YHat_Refined_Test   : float []
        }

    let createMetricsChart fileName stabwMedian (rnd:System.Random) (metrics:ModelMetrics) = 
        ///
        let traceName = (sprintf "#TestPeptides:%i Rsquared:%f RMS:%f, %s " metrics.X_Test.Length metrics.Metrics.RSquared metrics.Metrics.MeanAbsoluteError fileName)
        let color = getRandomColor rnd |> FSharpAux.Colors.toWebColor        
        let xVsY = 
            [
            Chart.Point(metrics.X_Test,metrics.Y_Test) |> Chart.withMarkerStyle(Color = color)
            Chart.Line(Array.zip metrics.X_Test metrics.YHat_Test |> Array.sortBy fst)
            ]
            |> Chart.Combine
            |> Chart.withX_AxisStyle("source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("measured (dotted) and predicted (Line) target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let yVsYHat = 
            Chart.Point(metrics.Y_Test,metrics.YHat_Test)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("target ScanTimes (Y_Test)")
            |> Chart.withY_AxisStyle("predicted target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let xVsDifferenceYandYHat = 
            Chart.Point(metrics.X_Test, Array.map2 (fun y yHat -> y - yHat) metrics.Y_Test metrics.YHat_Test)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("target ScanTimes (Y_Test) - predicted target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let xVsDifferenceYandYHatNormed = 
            Chart.Point(metrics.X_Test, Array.map2 (fun y yHat -> (y - yHat) / stabwMedian) metrics.Y_Test metrics.YHat_Test)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("normed target ScanTimes (Y_Test) - predicted target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let xVsDifferenceYandYHat_refined = 
            Chart.Point(metrics.X_Test, Array.map2 (fun y yHat -> y - yHat) metrics.Y_Test metrics.YHat_Refined_Test)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("target ScanTimes (Y_Test) - refined target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let xVsDifferenceYandYHatNormed_refined = 
            Chart.Point(metrics.X_Test, Array.map2 (fun y yHat -> (y - yHat) / stabwMedian) metrics.Y_Test metrics.YHat_Refined_Test)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("normed target ScanTimes (Y_Test) - refined target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let makeBar yHat_Test = 
            let metric =  Array.map2 (fun y yHat -> (y - yHat) / stabwMedian) metrics.Y_Test yHat_Test
            let UpToOne   = metric |> Array.filter (fun x -> abs x >= 0. && abs x < 1.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToTwo   = metric |> Array.filter (fun x -> abs x >= 1. && abs x < 2.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToThree = metric |> Array.filter (fun x -> abs x >= 2. && abs x < 3.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToFour  = metric |> Array.filter (fun x -> abs x >= 3. && abs x < 4.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToFive  = metric |> Array.filter (fun x -> abs x >= 4. && abs x < 5.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToSix   = metric |> Array.filter (fun x -> abs x >= 5. && abs x < 6.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToSeven = metric |> Array.filter (fun x -> abs x >= 6. && abs x < 7.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let OutOfSeven  = metric |> Array.filter (fun x -> abs x > 7.)                |> Array.length |> float |> fun x -> x / float metric.Length           
            [
                Chart.StackedColumn(values=[UpToOne],keys=[traceName],Name="UpToOne")       |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Green1)
                Chart.StackedColumn(values=[UpToTwo],keys=[traceName],Name="UpToTwo")       |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Green2)
                Chart.StackedColumn(values=[UpToThree],keys=[traceName],Name="UpToThree")   |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Green3)
                Chart.StackedColumn(values=[UpToFour],keys=[traceName],Name="UpToFour")     |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Orange1)
                Chart.StackedColumn(values=[UpToFive],keys=[traceName],Name="UpToFive")     |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Orange2)
                Chart.StackedColumn(values=[UpToSix],keys=[traceName],Name="UpToSix")       |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Orange3)
                Chart.StackedColumn(values=[UpToSeven],keys=[traceName],Name="UpToSeven")   |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Red1)
                Chart.StackedColumn(values=[OutOfSeven],keys=[traceName],Name="OutOfSeven") |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Red2)
            ]
            |> Chart.Combine
            |> Chart.withY_AxisStyle("")
        let bar = makeBar metrics.YHat_Test
        let barRef = makeBar metrics.YHat_Refined_Test
        [
        xVsY
        yVsYHat
        xVsDifferenceYandYHat
        xVsDifferenceYandYHatNormed
        xVsDifferenceYandYHat_refined
        xVsDifferenceYandYHatNormed_refined        
        bar
        barRef
        ]
        |> Chart.Stack 2
        |> Chart.withSize(2000.,2500.)
    
    type MetricsDTO = 
       {
           X_Intensities       : float 
           X_Stabw             : float         
           X_Test              : float 
           Y_Test              : float 
           YHat_Test           : float 
           YHat_Refined_Test   : float 
        }
    
    let saveMetrics outputDir fileName (metrics:ModelMetrics) = 
        let fileName = (fileName) + ".metric"
        let outFilePath = Path.Combine [|outputDir;fileName|]
        [|
        for i = 1 to metrics.X_Intensities.Length-1 do 
            yield
                {
                    X_Intensities       = metrics.X_Intensities.[i]
                    X_Stabw             = metrics.X_Stabw.[i]
                    X_Test              = metrics.X_Test.[i]
                    Y_Test              = metrics.Y_Test.[i]
                    YHat_Test           = metrics.YHat_Test.[i]
                    YHat_Refined_Test   = metrics.YHat_Refined_Test.[i]
                }
        |]
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
    ///
    let initAlign (ctx:MLContext) (pepsForLearning: PeptideForLearning []) = 
        let data = ctx.Data.LoadFromEnumerable(pepsForLearning)
        let split = ctx.Data.TrainTestSplit(data, testFraction= 0.1)
        let trainer = ctx.Regression.Trainers.Gam(featureColumnName="Features",labelColumnName="TargetScanTime")
        let pipeline =    
            (ctx.Transforms.Concatenate("Features","SourceScanTime")|> downcastPipeline)
              .Append(trainer)
        let model = pipeline.Fit(split.TrainSet)    
        let metrics =             
            let metrics = ctx.Regression.Evaluate(model.Transform(split.TestSet),labelColumnName="TargetScanTime")
            let evalView = model.Transform(split.TestSet)
            let i    = evalView.GetColumn<float32>("SourceIntensity")|> Seq.map float |> Array.ofSeq
            let std  = evalView.GetColumn<float32>("SourceStabw")    |> Seq.map float |> Array.ofSeq
            let x    = evalView.GetColumn<float32>("SourceScanTime") |> Seq.map float |> Array.ofSeq
            let y    = evalView.GetColumn<float32>("TargetScanTime") |> Seq.map float |> Array.ofSeq
            let yHat = evalView.GetColumn<float32>("Score")          |> Seq.map float |> Array.ofSeq
            let yHatAfterRefinement = 
                let xSource = evalView.GetColumn<float[]>("RtTrace_SourceFile")        |> Seq.map (Array.ofSeq) |> Array.ofSeq
                let ySource = evalView.GetColumn<float[]>("IntensityTrace_SourceFile") |> Seq.map (Array.ofSeq) |> Array.ofSeq
                let xTarget = evalView.GetColumn<float[]>("RtTrace_TargetFile")        |> Seq.map (Array.ofSeq) |> Array.ofSeq
                let yTarget = evalView.GetColumn<float[]>("IntensityTrace_TargetFile") |> Seq.map (Array.ofSeq) |> Array.ofSeq
                [|
                    for i = 0 to xSource.Length-1 do                          
                        printfn "%i %i %i %i" xSource.[i].Length ySource.[i].Length xTarget.[i].Length yTarget.[i].Length
                        let target = Array.zip xTarget.[i] (DTW'.zNorm yTarget.[i])
                        let source = Array.zip xSource.[i] (DTW'.zNorm ySource.[i])
                        let tmp = 
                            DTW'.align' target source x.[i] 
                            |> snd
                        yield tmp 
                |]                  
            {
                Metrics             = metrics
                X_Intensities       = i
                X_Stabw             = std      
                X_Test              = x
                Y_Test              = y
                YHat_Test           = yHat
                YHat_Refined_Test   = yHatAfterRefinement
            }
        let predF = ctx.Model.CreatePredictionEngine<PeptideForLearning,ScanTimePrediction>(model)
        let predict quantifiedPeptide = 
            quantifiedPeptide
            |> toPeptideForLearning None
            |> predF.Predict
            |> createAlignmentResult quantifiedPeptide
        metrics, predict
            
    ///
    let performAlignment outDir rnd align (target: AlignmentFile) (source: AlignmentFile) =
        
        ///
        let peptidesForLearning = 
            target.QuantifiedPeptides
            |> Seq.choose (fun tarQP -> 
                match Map.tryFind tarQP.Key source.QuantifiedPeptides with 
                | Some sourceQP ->
                    toPeptideForLearning (Some tarQP.Value) sourceQP
                    |> Some
                | None -> None
                ) 
            |> Array.ofSeq
            |> Array.shuffleFisherYates
        
        ///
        let metrics,model: ModelMetrics*(QuantificationResult->AlignmentResult) = 
            align peptidesForLearning
        
        ///
        let peptideIonsToTransfer,peptideIonsStillMissing = 
            target.MissingPeptides
            |> Array.fold (fun (pepsToTransfer,stillMissingPeps) missingPep -> 
                match Map.tryFind missingPep source.QuantifiedPeptides with 
                | Some pepToTransfer -> (pepToTransfer::pepsToTransfer,stillMissingPeps)
                | None -> (pepsToTransfer,missingPep::stillMissingPeps)
                ) ([],[])
        
        ///
        let alignmentResults : AlignmentResult [] = 
            peptideIonsToTransfer 
            |> List.toArray 
            |> Array.map model 
        
        ///
        let stabwMedian = 
            target.QuantifiedPeptides
            |> Seq.map (fun x -> getTargetStabw x.Value)
            |> Seq.filter (fun x -> nan.Equals x |> not)
            |> Seq.median

        ///
        saveMetrics outDir source.FileName metrics
        createMetricsChart source.FileName stabwMedian rnd metrics,
        {target with MissingPeptides = peptideIonsStillMissing |> Array.ofList; GainedPeptides = Array.append target.GainedPeptides alignmentResults}
        
    
    let alignFiles (logger:NLog.Logger) (processParams:AlignmentParams) (outputDir:string) (quantFiles:string) = 

        logger.Trace (sprintf "Input directory: %s" quantFiles)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        
        let getPlotFilePathFilePath (plotName:string) (fileName:string) =
            let fileName = (Path.GetFileNameWithoutExtension fileName) + "_" + plotName 
            Path.Combine [|outputDir;fileName|]
        
        logger.Trace "Init align function"
        let ctx = new ML.MLContext()
        let rnd = new System.Random()
        let align = initAlign ctx
        logger.Trace "Init align function: finished"
         
        logger.Trace "Reading and preparing .quant file for alignment"
        let quantFiles = System.IO.Directory.GetFiles (quantFiles, "*.quant")
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
                        let chart',tar' = performAlignment outputDir rnd align tar source
                        logger.Trace (sprintf "Performing Alignment %s vs %s: finished" tar.FileName source.FileName)
                        logger.Trace (sprintf "Missing peptides before: %i, now:%i" tar.MissingPeptides.Length tar'.MissingPeptides.Length)
                        chart'::charts,tar'
                        ) ([],target)
                chart
                |> Chart.Combine
                |> Chart.withTitle(target.FileName)
                |> Chart.SaveHtmlAs(getPlotFilePathFilePath "Metrics" target.FileName)                    
                result 
                )
        logger.Trace "Performing Alignment: finished"
        logger.Trace "Writing Results"
        alignments
        |> Array.iter (fun tar -> 
            let outFilePath =
                let fileName = (tar.FileName) + ".align"
                Path.Combine [|outputDir;fileName|]
            logger.Trace (sprintf "outFilePath:%s" outFilePath)
            tar.GainedPeptides
            |> SeqIO'.csv "\t" true false
            |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
            )
        logger.Trace "Writing Results:finished"
                
        

        




































