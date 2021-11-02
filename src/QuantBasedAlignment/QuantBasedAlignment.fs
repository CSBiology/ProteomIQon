namespace ProteomIQon


open System.IO
open System
open FSharpAux.Colors.Table.StatisticalGraphics24
open FSharp.Stats
open FSharpAux.IO.SchemaReader
open Plotly.NET
open BioFSharp
open Microsoft
open Microsoft.ML
open Microsoft.ML.Data
open Microsoft.ML.AutoML   
open Dto
open Dto.QuantificationResult

module QuantBasedAlignment = 

    ///
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

    ///
    let downcastPipeline (x : IEstimator<_>) = 
        match x with 
        | :? IEstimator<ITransformer> as y -> y
        | _ -> failwith "downcastPipeline: expecting a IEstimator<ITransformer>"

    ///
    type PeptideIon = 
        {
            Sequence             : string
            GlobalMod            : int
            Charge               : int        
        }
    
    ///
    let toPeptideIon (qp:QuantificationResult) = 
        {
            Sequence             = qp.StringSequence
            GlobalMod            = qp.GlobalMod
            Charge               = qp.Charge
        }

    ///
    type AlignmentFile = {
        FileName                   : string
        QuantifiedPeptides         : Map<PeptideIon,QuantificationResult> 
        MissingPeptides            : PeptideIon []
        GainedPeptides             : AlignmentResult []
        }
    

    [<CLIMutable>]
    type PeptideForLearning = 
        {
            [<ColumnName("Sequence")>]
            Sequence                     : string
            [<ColumnName("GlobalMod")>]
            GlobalMod                    : string
            [<ColumnName("Charge")>]
            Charge                       : string
            [<ColumnName("PepSequenceID")>]
            PepSequenceID                : string
            [<ColumnName("ModSequenceID")>]
            ModSequenceID                : string
            [<ColumnName("SourceScanTime")>]
            SourceScanTime               : float32
            [<ColumnName("TargetScanTime")>]
            TargetScanTime               : float32
            [<ColumnName("ScanTimeDifference")>]
            ScanTimeDifference           : float32
        }
    
    [<CLIMutable>]
    type PeptideComplete = 
        {
            [<ColumnName("Sequence")>]
            Sequence                     : string
            [<ColumnName("GlobalMod")>]
            GlobalMod                    : string
            [<ColumnName("Charge")>]
            Charge                       : string
            [<ColumnName("PepSequenceID")>]
            PepSequenceID                : string
            [<ColumnName("ModSequenceID")>]
            ModSequenceID                : string
            [<ColumnName("SourceScanTime")>]
            SourceScanTime               : float32
            [<ColumnName("SourceIntensity")>]
            SourceIntensity              : float32
            [<ColumnName("SourceStabw")>]
            SourceStabw                  : float32
            [<ColumnName("TargetScanTime")>]
            TargetScanTime               : float32
            [<ColumnName("TargetIntensity")>]
            TargetIntensity              : float32
            [<ColumnName("RtTrace_SourceFile")>]
            RtTrace_SourceFile                              : float [] 
            [<ColumnName("IntensityTrace_SourceFile")>]
            IntensityTrace_SourceFile                       : float []
            [<ColumnName("RtTrace_TargetFile")>]
            RtTrace_TargetFile                              : float []
            [<ColumnName("IntensityTrace_TargetFile")>]
            IntensityTrace_TargetFile                       : float []          
            [<ColumnName("IsotopicPatternMz_SourceFile")>]
            IsotopicPatternMz_SourceFile                    : float []          
            [<ColumnName("IsotopicPatternIntensity_Observed_SourceFile")>]
            IsotopicPatternIntensity_Observed_SourceFile    : float []         
            [<ColumnName("IsotopicPatternMz_TargetFile")>]
            IsotopicPatternMz_TargetFile                    : float []         
            [<ColumnName("IsotopicPatternIntensity_Observed_TargetFile")>]
            IsotopicPatternIntensity_Observed_TargetFile    : float []
        }

    /////
    //let formatString s = String.filter (fun x -> Char.IsUpper x) s

    ///
    let toPeptideForLearning (targetPep:QuantificationResult option) (sourcePep:QuantificationResult) = 
        match targetPep with 
        | Some tP -> 
            {
                Sequence                                        = (*formatString*) sourcePep.StringSequence
                GlobalMod                                       = sourcePep.GlobalMod          |> string      
                Charge                                          = sourcePep.Charge             |> string
                PepSequenceID                                   = sourcePep.PepSequenceID      |> string
                ModSequenceID                                   = sourcePep.ModSequenceID      |> string
                SourceScanTime                                  = getTargetScanTime sourcePep  |> float32
                TargetScanTime                                  = getTargetScanTime tP         |> float32
                ScanTimeDifference                              = ((getTargetScanTime tP) - (getTargetScanTime sourcePep)) |> float32
            },
            {
                Sequence                                        = (*formatString*) sourcePep.StringSequence
                GlobalMod                                       = sourcePep.GlobalMod          |> string      
                Charge                                          = sourcePep.Charge             |> string
                PepSequenceID                                   = sourcePep.PepSequenceID      |> string
                ModSequenceID                                   = sourcePep.ModSequenceID      |> string
                SourceScanTime                                  = getTargetScanTime sourcePep  |> float32
                SourceIntensity                                 = getTargetIntensity sourcePep |> float32
                SourceStabw                                     = getTargetStabw sourcePep     |> float32
                TargetScanTime                                  = getTargetScanTime tP         |> float32
                TargetIntensity                                 = getTargetIntensity tP        |> float32
                RtTrace_SourceFile                              = getTargetRtTrace sourcePep
                IntensityTrace_SourceFile                       = getTargetIntensityTrace sourcePep
                RtTrace_TargetFile                              = getTargetRtTrace tP
                IntensityTrace_TargetFile                       = getTargetIntensityTrace tP
                IsotopicPatternMz_SourceFile                    = getIsotopicPatternMz sourcePep
                IsotopicPatternIntensity_Observed_SourceFile    = getIsotopicPatternIntensity_Observed sourcePep
                IsotopicPatternMz_TargetFile                    = getIsotopicPatternMz tP  
                IsotopicPatternIntensity_Observed_TargetFile    = getIsotopicPatternIntensity_Observed tP
            }
        | None -> 
            {
                Sequence                                        = (*formatString*) sourcePep.StringSequence
                GlobalMod                                       = sourcePep.GlobalMod          |> string   
                Charge                                          = sourcePep.Charge             |> string
                PepSequenceID                                   = sourcePep.PepSequenceID      |> string
                ModSequenceID                                   = sourcePep.ModSequenceID      |> string
                SourceScanTime                                  = getTargetScanTime sourcePep  |> float32
                TargetScanTime                                  = nan                          |> float32
                ScanTimeDifference                              = nan                          |> float32
            },
            {
                Sequence                                        = (*formatString*) sourcePep.StringSequence
                GlobalMod                                       = sourcePep.GlobalMod          |> string   
                Charge                                          = sourcePep.Charge             |> string
                PepSequenceID                                   = sourcePep.PepSequenceID      |> string
                ModSequenceID                                   = sourcePep.ModSequenceID      |> string
                SourceScanTime                                  = getTargetScanTime sourcePep  |> float32
                SourceIntensity                                 = getTargetIntensity sourcePep |> float32
                SourceStabw                                     = getTargetStabw sourcePep     |> float32
                TargetScanTime                                  = nan                          |> float32
                TargetIntensity                                 = nan                          |> float32
                RtTrace_SourceFile                              = getTargetRtTrace sourcePep
                IntensityTrace_SourceFile                       = getTargetIntensityTrace sourcePep
                RtTrace_TargetFile                              = [||]
                IntensityTrace_TargetFile                       = [||]
                IsotopicPatternMz_SourceFile                    = getIsotopicPatternMz sourcePep
                IsotopicPatternIntensity_Observed_SourceFile    = getIsotopicPatternIntensity_Observed sourcePep
                IsotopicPatternMz_TargetFile                    = [||]
                IsotopicPatternIntensity_Observed_TargetFile    = [||]
            }

    ///
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
    let createAlignmentFiles (sourceFiles :string []) (targetFile:string) = 
        let toAlignmentFile (allPepIons:PeptideIon []) (filePath:string) peptides = 
            let presentPeptides = 
                peptides 
                |> Array.map (fun qp -> toPeptideIon qp, qp)
                |> Map.ofArray
            let missingPeptides = 
                allPepIons
                |> Array.filter (fun pep -> 
                    let pepUnMod = {pep with GlobalMod = 0}
                    let pepMod = {pep with GlobalMod = 1}
                    (presentPeptides.ContainsKey pepMod || presentPeptides.ContainsKey pepUnMod) 
                    |> not
                    )
            {
                FileName                   = Path.GetFileNameWithoutExtension(filePath)
                QuantifiedPeptides         = presentPeptides 
                MissingPeptides            = missingPeptides 
                GainedPeptides             = [||]
            }
        let targetPeptides = 
            targetFile 
            |> getQuantifiedPeptides
        let sourcePeptides = 
            sourceFiles
            |> Array.map getQuantifiedPeptides
        let allPeptideIons =
            sourcePeptides
            |> Array.append [|targetPeptides|] 
            |> Array.concat
            |> Array.map (fun qp -> toPeptideIon qp)
            |> Array.distinct
        let targetAlignmentFile = toAlignmentFile allPeptideIons targetFile targetPeptides
        let sourceAlignmentFiles = Array.map2 (toAlignmentFile allPeptideIons) sourceFiles sourcePeptides
        sourceAlignmentFiles, targetAlignmentFile 

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

    // /// determines order of alignment based on the calculated file difference. 
    // // calculateFileDifference can be easily factored out if needed. 
    // let findAlignmentOrder (alignmentFiles: AlignmentFile []) =
    //     alignmentFiles
    //     |> Array.mapi (fun i af -> 
    //             af,
    //             [| 
    //                 for j = 0 to alignmentFiles.Length-1 do 
    //                     if j <> i then calculateFileDifference af alignmentFiles.[j], alignmentFiles.[j]
    //             |]
    //             |> Array.sortBy fst
    //         )
    
    ///
    let createAlignmentResult (quantifiedPeptide:QuantificationResult) (scanTimePrediction:ScanTimePrediction) = 
        {
            StringSequence                                  = quantifiedPeptide.StringSequence
            GlobalMod                                       = quantifiedPeptide.GlobalMod
            Charge                                          = quantifiedPeptide.Charge
            PepSequenceID                                   = quantifiedPeptide.PepSequenceID
            ModSequenceID                                   = quantifiedPeptide.ModSequenceID
            Mz                                              = Mass.toMZ (quantifiedPeptide.TheoMass) (float quantifiedPeptide.Charge) 
            ProteinNames                                    = quantifiedPeptide.ProteinNames
            PredictedScanTime                               = float scanTimePrediction.TargetScanTime
            ScanTime_SourceFile                             = getTargetScanTime quantifiedPeptide
            RtTrace_SourceFile                              = getTargetRtTrace quantifiedPeptide
            IntensityTrace_SourceFile                       = getTargetIntensityTrace quantifiedPeptide
            IsotopicPatternMz_SourceFile                    = getIsotopicPatternMz quantifiedPeptide       
            IsotopicPatternIntensity_Observed_SourceFile    = getIsotopicPatternIntensity_Observed quantifiedPeptide       
        }   


    ///
    type ModelMetrics = 
        {
        RSquared                             : float
        Sequence                            : string []
        GlobalMod                           : int []
        Charge                              : int []
        PepSequenceID                       : int []
        ModSequenceID                       : int []
        X_Intensities                       : float []
        X_Stabw                             : float []        
        X_Test                              : float []
        X_IsotopicPatternMz                 : float [][]
        X_IsotopicPatternIntensity_Observed : float [][]
        X_RtTrace                           : float [][]
        X_IntensityTrace                    : float [][]   
        Y_Test                              : float []
        YHat_Test                           : float []
        YHat_Refined_Test                   : float []
        Y_Intensities                       : float []
        Y_IsotopicPatternMz                 : float [][]
        Y_IsotopicPatternIntensity_Observed : float [][]
        Y_RtTrace                           : float [][]
        Y_IntensityTrace                    : float [][]
        DtwDistanceBefore                   : float []
        }



    ///
    let createMetricsChart fileName (*stabwMedian*) (rnd:System.Random) (metrics:ModelMetrics) = 
        ///
        let traceName = ""//(sprintf "#TestPeptides:%i Rsquared:%f RMS:%f, %s " metrics.X_Test.Length metrics.Metrics.RSquared metrics.Metrics.MeanAbsoluteError fileName)
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
            Chart.Point(metrics.X_Test, Array.map3 (fun y yHat stabwMedian -> (y - yHat) / stabwMedian) metrics.Y_Test metrics.YHat_Test metrics.X_Stabw)
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
            Chart.Point(metrics.X_Test, Array.map3 (fun y yHat stabwMedian -> (y - yHat) / stabwMedian) metrics.Y_Test metrics.YHat_Refined_Test metrics.X_Stabw)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("normed target ScanTimes (Y_Test) - refined target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let makeBar yHat_Test = 
            let metric =  Array.map3 (fun y yHat stabwMedian -> (y - yHat) / stabwMedian) metrics.Y_Test yHat_Test metrics.X_Stabw
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


    ///
    let saveMetrics outputDir targetFileName sourceFileName (metrics:ModelMetrics) = 
        let fileName = (targetFileName ) + ".alignmetric"
        let outFilePath = Path.Combine [|outputDir;fileName|]
        let data = 
            [|
            for i = 1 to metrics.X_Intensities.Length-1 do 
                yield
                    {           
                        Sequence                             = metrics.Sequence.[i]
                        GlobalMod                            = metrics.GlobalMod.[i]
                        Charge                               = metrics.Charge.[i]
                        PepSequenceID                        = metrics.PepSequenceID.[i]
                        ModSequenceID                        = metrics.ModSequenceID.[i]
                        X_FileName                           = sourceFileName 
                        X_Intensities                        = metrics.X_Intensities.[i] 
                        X_Stabw                              = metrics.X_Stabw.[i]
                        X_Test                               = metrics.X_Test.[i] 
                        X_IsotopicPatternMz                  = metrics.X_IsotopicPatternMz.[i]
                        X_IsotopicPatternIntensity_Observed  = metrics.X_IsotopicPatternIntensity_Observed.[i]
                        X_RtTrace                            = metrics.X_RtTrace.[i]
                        X_IntensityTrace                     = metrics.X_IntensityTrace.[i]
                        Y_Test                               = metrics.Y_Test.[i] 
                        YHat_Test                            = metrics.YHat_Test.[i] 
                        YHat_Refined_Test                    = metrics.YHat_Refined_Test.[i]
                        Y_Intensity                          = metrics.Y_Intensities.[i]
                        Y_IsotopicPatternMz                  = metrics.Y_IsotopicPatternMz.[i]
                        Y_IsotopicPatternIntensity_Observed  = metrics.Y_IsotopicPatternIntensity_Observed.[i]
                        Y_RtTrace                            = metrics.Y_RtTrace.[i]
                        Y_IntensityTrace                     = metrics.Y_IntensityTrace.[i]   
                        DtwDistanceBefore                    = metrics.DtwDistanceBefore.[i]
                    }
            |]
        if System.IO.File.Exists outFilePath then 
            data
            |> SeqIO'.csv "\t" false false
            |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        else
            data
            |> SeqIO'.csv "\t" true false
            |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)

    ///
    let initAlign (logger:NLog.Logger) (ctx:MLContext) (pepsForLearning: (PeptideForLearning*PeptideComplete) []) = 
        logger.Trace ("Sampling test and train data")
        let train,test = 
            let getBinIdx width scantime = int ((scantime / width))    
            let knotX,train,test = 
                pepsForLearning
                |> Array.groupBy (fun (pl,_) -> getBinIdx 10. (pl.SourceScanTime |> float))
                |> Array.map (fun (binIdx,ions) -> 
                    let dataS = ions |> Array.shuffleFisherYates
                    let knotX = 
                        dataS 
                        |> Array.map fst
                        |> Array.maxBy (fun x -> x.SourceScanTime)
                    let train,test =
                        let nTest = 
                            (float dataS.Length) * 0.9
                            |> int
                        dataS.[.. nTest], dataS.[nTest+1 ..] 
                    float knotX.SourceScanTime, train, test
                    )
                |> Array.unzip3
            train |> Array.concat |> Array.sortBy (fun (pepLearning,pepComp) -> pepLearning.SourceScanTime), test |> Array.concat |> Array.sortBy (fun (pepLearning,pepComp) -> pepLearning.SourceScanTime)
        let trainLearn,trainComp = train |> Array.unzip
        let testLearn,testComp = test |> Array.unzip
        logger.Trace ("Sampling test and train data:finished")
        
        logger.Trace ("Training Fast Tree")
        let trainView = ctx.Data.LoadFromEnumerable(trainLearn)
        // let testView = ctx.Data.LoadFromEnumerable(testLearn)
        let trainer = ctx.Regression.Trainers.Gam(featureColumnName="Features",labelColumnName="TargetScanTime")
        let pipeline =    
            (ctx.Transforms.Concatenate("Features","SourceScanTime")|> downcastPipeline)
              .Append(trainer)
        let model = pipeline.Fit(trainView) 
        logger.Trace ("Training Fast Tree:finished")        
        logger.Trace ("Training Spline")
        let knots = 
            let upB = model.LastTransformer.Model.GetBinUpperBounds 0
            let bounds = (upB |> Array.ofSeq |> fun x -> x.[0..x.Length-2])
            bounds 
            |> Array.choose (fun b -> 
                match trainLearn |> Array.tryFind (fun x -> float x.SourceScanTime > b) with 
                | Some x -> x.SourceScanTime |> float |> Some
                | None -> None
                )   
            |> Array.mapi (fun i x -> if i%2 = 0 then Some x else None)
            |> Array.choose id   
        logger.Trace (sprintf "Number of knots selected:%i" knots.Length)
        // let spline = ProteomIQon.Spline.smoothingSpline (train |> Array.map (fun (x,y) -> float x.SourceScanTime, float x.TargetScanTime)) (knots)
        let spline = FSharp.Stats.Fitting.Spline.smoothingSpline (train |> Array.map (fun (x,y) -> float x.SourceScanTime, float x.TargetScanTime)) (knots)        
        let trainer lambda = 
            let fit = spline lambda 
            let x,y,yHat = 
                train
                |> Array.map (fun (x,y) -> 
                    x, x.TargetScanTime, fit (float x.SourceScanTime)
                    )
                |> Array.unzip3
            let rS = FSharp.Stats.Fitting.GoodnessOfFit.calculateDeterminationFromValue yHat (Array.map float y)
            lambda, rS, fit                  
        logger.Trace ("Optimizing Lambdas")
        let models =      
            [|0.001 .. 0.01 .. 1.|]
            |> Array.map trainer
        let lambda,rSquared,model = 
            models
            |> Array.maxBy (fun (x,y,z) -> y )                  
        logger.Trace (sprintf "Optimizing Lambdas:Finished, selected lamda:%f" lambda)
        [
            train
            |> Array.map (fun (x,y) -> 
                    x.SourceScanTime, (float x.TargetScanTime)
                    )
            |> Chart.Point    
            Chart.Point(knots,knots)
            models
            |> Array.map (fun (x,y,z) -> 
                train
                |> Array.map (fun (x,y) -> 
                        x.SourceScanTime, z (float x.SourceScanTime)
                        )
                |> Chart.Line
                |> Chart.withTraceName (sprintf "lambda %f" x)
                )
            |> Chart.Combine
        ]
        |> Chart.Combine
        |> Chart.Show
        models
        |> Array.map (fun (x,y,z) -> x,y)
        |> Chart.Point
        |> Chart.Show
        logger.Trace ("Training Spline:Finished")
        let metrics =             
            let rSquared = rSquared
            let yHat          = testComp |> Seq.map (fun x -> model (float x.SourceScanTime))|> Array.ofSeq
            let sequence      = testComp |> Seq.map (fun x -> x.Sequence)       |> Array.ofSeq
            let globalMod     = testComp |> Seq.map (fun x -> x.GlobalMod)      |> Array.ofSeq
            let charge        = testComp |> Seq.map (fun x -> x.Charge)         |> Array.ofSeq
            let pepSequenceID = testComp |> Seq.map (fun x -> x.PepSequenceID)  |> Array.ofSeq
            let modSequenceID = testComp |> Seq.map (fun x -> x.ModSequenceID)  |> Array.ofSeq
            let i             =  testComp |> Seq.map (fun x -> float x.SourceIntensity)  |> Array.ofSeq
            let std           =  testComp |> Seq.map (fun x -> float x.SourceStabw)      |> Array.ofSeq
            let x             =  testComp |> Seq.map (fun x -> float x.SourceScanTime)   |> Array.ofSeq
            let y             =  testComp |> Seq.map (fun x -> float x.TargetScanTime)   |> Array.ofSeq
            let yIntensities  =  testComp |> Seq.map (fun x -> float x.TargetIntensity)  |> Array.ofSeq
            let xSource = testComp |> Seq.map (fun x -> Array.ofSeq x.RtTrace_SourceFile)         |> Array.ofSeq
            let ySource = testComp |> Seq.map (fun x -> Array.ofSeq x.IntensityTrace_SourceFile)  |> Array.ofSeq
            let xTarget = testComp |> Seq.map (fun x -> Array.ofSeq x.RtTrace_TargetFile)         |> Array.ofSeq
            let yTarget = testComp |> Seq.map (fun x -> Array.ofSeq x.IntensityTrace_TargetFile)  |> Array.ofSeq     
            let yHatAfterRefinement,dtwDistanceBefore,dtwDistanceAfter = 
                [|
                    for i = 0 to xSource.Length-1 do                         
                        printfn "%i %i %i %i" xSource.[i].Length ySource.[i].Length xTarget.[i].Length yTarget.[i].Length
                        let target = Array.zip xTarget.[i] (DTW'.zNorm yTarget.[i])
                        let source = Array.zip xSource.[i] (DTW'.zNorm ySource.[i])
                        let yRefined = 
                            DTW'.align' target source x.[i] 
                            |> snd
                        let alignment = 
                            DTW'.align target source 
                            |> Array.ofList
                        let dtwDistanceBefore = 
                            DTW'.distance None None None None None None (Array.map snd target)  (Array.map snd source)
                        let dtwDistanceAfter = 
                            DTW'.distance None None None None None None (Array.map snd target)  (Array.map snd alignment)
                        yield yRefined,dtwDistanceBefore,dtwDistanceAfter
                |]
                |> Array.unzip3

            let x_IsotopicPatternMz                 = testComp |> Seq.map (fun x -> Array.ofSeq x.IsotopicPatternMz_SourceFile)                 |> Array.ofSeq
            let x_IsotopicPatternIntensity_Observed = testComp |> Seq.map (fun x -> Array.ofSeq x.IsotopicPatternIntensity_Observed_SourceFile) |> Array.ofSeq
            let y_IsotopicPatternMz                 = testComp |> Seq.map (fun x -> Array.ofSeq x.IsotopicPatternMz_TargetFile)                 |> Array.ofSeq
            let y_IsotopicPatternIntensity_Observed = testComp |> Seq.map (fun x -> Array.ofSeq x.IsotopicPatternIntensity_Observed_TargetFile) |> Array.ofSeq       
            {
                RSquared                            = rSquared
                Sequence                            = sequence                           
                GlobalMod                           = globalMod     |> Array.map int                     
                Charge                              = charge        |> Array.map int                     
                PepSequenceID                       = pepSequenceID |> Array.map int                     
                ModSequenceID                       = modSequenceID |> Array.map int                     
                X_Intensities                       = i
                X_Stabw                             = std      
                X_Test                              = x
                X_IsotopicPatternMz                 = x_IsotopicPatternMz
                X_IsotopicPatternIntensity_Observed = x_IsotopicPatternIntensity_Observed
                X_RtTrace                           = xSource
                X_IntensityTrace                    = ySource
                Y_Test                              = y
                YHat_Test                           = yHat
                YHat_Refined_Test                   = yHatAfterRefinement
                Y_Intensities                       = yIntensities
                Y_IsotopicPatternMz                 = y_IsotopicPatternMz
                Y_IsotopicPatternIntensity_Observed = y_IsotopicPatternIntensity_Observed
                Y_RtTrace                           = xTarget
                Y_IntensityTrace                    = yTarget
                DtwDistanceBefore                   = dtwDistanceBefore
            }

        // Spline
        let predict quantifiedPeptide = 
            quantifiedPeptide
            |> toPeptideForLearning None
            |> fun (pepToLearn,pepComp) -> 
                let v = model (float pepToLearn.SourceScanTime)        
                {
                TargetScanTime = float32 v
                }
            |> createAlignmentResult quantifiedPeptide

        metrics, predict

    type Alignment = {
        Metrics      : ModelMetrics
        MetricsChart : GenericChart.GenericChart
        AlignFunc    : (QuantificationResult->AlignmentResult)
        SourceFile   : AlignmentFile
        }
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

        {
        Metrics      = metrics
        MetricsChart = createMetricsChart source.FileName (*stabwMedian*) rnd metrics
        AlignFunc    = model
        SourceFile   = source
        }
        
    ///
    let alignFiles diagCharts (processParams:AlignmentParams) (outputDir:string) (sourceFiles:string []) (targetFile:string) = 
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension targetFile)
        logger.Trace (sprintf "target file: %A" targetFile)
        logger.Trace (sprintf "source files: %A" sourceFiles)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        let getPlotFilePathFilePath (plotName:string) (fileName:string) =
            let fileName = (Path.GetFileNameWithoutExtension fileName) + "_" + plotName 
            Path.Combine [|outputDir;fileName|]
        logger.Trace "Init align function"
        let ctx = ML.MLContext()
        let rnd = System.Random()
        let align = initAlign logger ctx
        logger.Trace "Init align function: finished"
         
        logger.Trace "Reading and preparing .quant files for alignment"
        let sourceAlignmentFiles, targetAlignmentFile  = createAlignmentFiles sourceFiles targetFile
        logger.Trace "Reading and preparing .quant files for alignment: finished"

        logger.Trace "Performing Alignments"
        let alignmentsSortedByQuality = 
            let alignments =
                sourceAlignmentFiles
                |> Array.map (fun (sourceFile)  -> 
                    logger.Trace (sprintf "Performing Alignment %s vs %s" targetAlignmentFile.FileName sourceFile.FileName)
                    let alignResult = performAlignment outputDir rnd align targetAlignmentFile sourceFile
                    saveMetrics outputDir targetAlignmentFile.FileName sourceFile.FileName alignResult.Metrics
                    logger.Trace (sprintf "Performing Alignment %s vs %s: finished" targetAlignmentFile.FileName sourceFile.FileName)
                    alignResult
                    ) 
            let sortedByQuality = 
                alignments
                |> Array.sortByDescending (fun x -> x.Metrics.RSquared)
            if diagCharts then 
                sortedByQuality
                |> Array.map (fun x -> x.MetricsChart)
                |> Chart.Combine
                |> Chart.withTitle(targetAlignmentFile.FileName)
                |> Chart.SaveHtmlAs(getPlotFilePathFilePath "Metrics" targetAlignmentFile.FileName)                    
            sortedByQuality 
        logger.Trace "Performing Alignments: finished"
        // if diagCharts then 
        //     logger.Trace "Plotting file distances"
        //     let chart = 
        //         alignmentFilesOrdered
        //         |> Array.map (fun (target,sources) ->
        //                 Chart.Point(sources |> Array.mapi (fun i x -> (snd x).FileName, fst x))
        //                 |> Chart.withTraceName target.FileName
        //                 |> Chart.withX_AxisStyle("FileNames")
        //                 |> Chart.withY_AxisStyle("Median absolute difference of peptide ion scan times")
        //                 |> Chart.withSize(1000.,1000.)
        //                 |> Chart.SaveHtmlAs(getPlotFilePathFilePath "differences" target.FileName)
        //             )
        //     logger.Trace "Plotting file distances: finished"
        logger.Trace "Transfer identifications"
        let result = 
            alignmentsSortedByQuality
            |> Array.fold (fun target alignment -> 
                    let peptideIonsToTransfer,peptideIonsStillMissing = 
                        target.MissingPeptides
                        |> Array.fold (fun (pepsToTransfer,stillMissingPeps) missingPep -> 
                            match Map.tryFind missingPep alignment.SourceFile.QuantifiedPeptides with 
                            | Some pepToTransfer -> (pepToTransfer::pepsToTransfer,stillMissingPeps)
                            | None -> (pepsToTransfer,missingPep::stillMissingPeps)
                            ) ([],[])
                    let alignmentResults : AlignmentResult [] = 
                        peptideIonsToTransfer 
                        |> List.toArray 
                        |> Array.map alignment.AlignFunc 
                    {target with MissingPeptides = peptideIonsStillMissing |> Array.ofList; GainedPeptides = Array.append target.GainedPeptides alignmentResults}
                ) targetAlignmentFile
        logger.Trace "Writing Results"
        let outFilePath =
            let fileName = (result.FileName) + ".align"
            Path.Combine [|outputDir;fileName|]
        logger.Trace (sprintf "outFilePath:%s" outFilePath)
        result.GainedPeptides
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        logger.Trace "Writing Results:finished"
                
        

        




































