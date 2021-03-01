namespace ProteomIQon

open System.IO
open System
open FSharpAux.Colors.Table.StatisticalGraphics24
open FSharpAux.IO
open FSharp.Stats
open FSharpAux.IO.SchemaReader
open FSharp.Plotly
open BioFSharp
open Microsoft
open Microsoft.ML
open Microsoft.ML.Data   
open Dto
open Dto.QuantificationResult
open Newtonsoft.Json
open MzIO.IO
open MzIO.Processing
open MzIO.Processing.MzIOLinq
open MzIO.MetaData.PSIMSExtension
open MzIO.MetaData.ParamEditExtension
open MzIO.MetaData.UO.UO
open MzIO.Binary
open MzIO.Model
open MzIO.Model.CvParam
open MzIO.Commons.Arrays
open System.Linq

module ConsensusSpectralLibrary =
    
    /////
    //type SparsePeakArray = { 
    //    Length: int
    //    Data: System.Collections.Generic.IDictionary<int,float>
    //    }
   
    /////
    //let dot (x:SparsePeakArray) (y:SparsePeakArray) =
    //    x.Data
    //    |> Seq.fold (fun (acc:float) xi -> 
    //        let present,yi = y.Data.TryGetValue xi.Key
    //        if present then acc + (yi * xi.Value) else acc
    //        ) 0.
           
    /////
    //let getBinIdx' width offset x = int ((x / width) + offset)

    /////
    //let peaksToNearestBinVector binWidth offset (pkarr:BioFSharp.Mz.PeakArray<_>) (minMassBoarder:float) (maxMassBoarder:float) = 
    //    let minMassBoarder = getBinIdx' binWidth offset minMassBoarder
    //    let maxMassBoarder = getBinIdx' binWidth offset maxMassBoarder
    //    let maxIndex = maxMassBoarder - minMassBoarder + 1        
    //    let keyValues = 
    //        pkarr 
    //        |> Array.choose (fun p ->  
    //            let index = (getBinIdx' binWidth offset p.Mz) - minMassBoarder
    //            if index < maxIndex-1 && index > -1 then Some (index, p.Intensity) else None
    //            )
    //        |> Array.groupBy fst 
    //        |> Array.map (fun (idx,data) -> idx, data |> Array.sumBy snd)
    //    { 
    //        Length = maxIndex
    //        Data = keyValues |> dict
    //    } 

    ///
    [<CLIMutable>]
    type PeptideForLearning = 
        {
            [<ColumnName("Sequence")>]
            Sequence                     : string
            [<ColumnName("GlobalMod")>]
            GlobalMod                    : int
            [<ColumnName("Charge")>]
            Charge                       : int
            [<ColumnName("PepSequenceID")>]
            PepSequenceID                : int
            [<ColumnName("ModSequenceID")>]
            ModSequenceID                : int
            [<ColumnName("SourceScanTime")>]
            SourceScanTime               : float32
            [<ColumnName("TargetScanTime")>]
            TargetScanTime               : float32
        }

    let toPeptideForLearning targetScanTime pepIon = 
        match targetScanTime with 
        | Some t -> 
            {
                Sequence                     = pepIon.StringSequence     
                GlobalMod                    = pepIon.GlobalMod    
                Charge                       = pepIon.Charge       
                PepSequenceID                = pepIon.PepSequenceID
                ModSequenceID                = pepIon.ModSequenceID
                SourceScanTime               = float32 pepIon.ScanTime
                TargetScanTime               = float32 t
            }
        | _ -> 
            {
                Sequence                     = pepIon.StringSequence     
                GlobalMod                    = pepIon.GlobalMod    
                Charge                       = pepIon.Charge       
                PepSequenceID                = pepIon.PepSequenceID
                ModSequenceID                = pepIon.ModSequenceID
                SourceScanTime               = float32 pepIon.ScanTime
                TargetScanTime               = float32 nan
            }

    ///
    type Library =
        {
        FileName: string
        Data: ProteomIQon.Dto.PeptideIon[]
        }

    ///
    let readLibraryFrom path= 
        let tmp: ProteomIQon.Dto.PeptideIon[] = ProteomIQon.Json.ReadAndDeserialize path
        {
            FileName = Path.GetFileNameWithoutExtension path
            Data = tmp
        }

    ///
    type LibToAlign = 
        {
            FileName                    : string
            PeptideIons                 : ProteomIQon.Dto.PeptideIon[]
            IdentifiedPeptideFeatures   : PeptideForLearning []
        }

    ///
    let createLibToAlign (lib:Library) (identifiedFeatures: PeptideForLearning [])= 
        {
            FileName                    = lib.FileName
            PeptideIons                 = lib.Data
            IdentifiedPeptideFeatures   = identifiedFeatures
        }

    type RegressionMetrics = {
        RSquared: float
        }
        
    ///
    type AlignedLib = 
        {
            FileName                    : string
            PeptideIons                 : ProteomIQon.Dto.PeptideIon[]
            IdentifiedPeptideFeatures   : PeptideForLearning []
            RegressionMetric            : RegressionMetrics
            AlignedPeptideIons          : ProteomIQon.Dto.PeptideIon[]
        }

    ///
    let createAlignedLib (lib:LibToAlign) regressionMetrics alignedIons = 
        {
            FileName                    = lib.FileName
            PeptideIons                 = lib.PeptideIons
            IdentifiedPeptideFeatures   = lib.IdentifiedPeptideFeatures
            RegressionMetric            = regressionMetrics 
            AlignedPeptideIons          = alignedIons
        }

    type IonSpecies = {
        ModSequenceID                               : int
        Charge                                      : int
        GlobalMod                                   : int
        }
    
    let createIonSpecies modSequenceID charge globalMod = {
        ModSequenceID  = modSequenceID
        Charge         = charge
        GlobalMod      = globalMod
        }

    ///
    [<CLIMutable>]
    type ScanTimePrediction = 
        {
            [<ColumnName("Score")>]
            TargetScanTime : float32
        }

    ///
    type SwathIndexer.SwathIndexer with
            
        member this.GetRTProfilesFirstWnd(dataReader:IMzIODataReader, query: SwathQuery, getLockMz: bool, spectrumSelector: seq<SwathIndexer.MSSwath> -> seq<SwathIndexer.MSSwath> list ,?mzRangeSelector: Peak1DArray * RangeQuery -> Peak1D) =

            let getClosestMz (peaks: Peak1DArray, mzRange: RangeQuery) =
                peaks.Peaks
                    .DefaultIfEmpty(new Peak1D(0., mzRange.LockValue))
                    .ItemAtMin(fun x -> Math.Abs(x.Mz - mzRange.LockValue))

            let mzRangeSelector = defaultArg mzRangeSelector getClosestMz

            let swathSpectra = 
                this.SwathList.SearchAllTargetMz(query.TargetMz)
                |> spectrumSelector

            swathSpectra
            |> List.map (fun spec ->
                let swathSpectrum = spec.SelectMany(fun x -> x.SearchAllRt(query)).ToArray()

                if swathSpectrum.Length > 0 then

                    let profile = Array2D.create query.CountMS2Masses swathSpectrum.Length (new Peak2D())

                    for specIdx = 0 to swathSpectrum.Length - 1 do

                        let swathSpec = swathSpectrum.[specIdx]
                        let pa = dataReader.ReadSpectrumPeaks(swathSpec.SpectrumID)

                        for ms2MassIndex = 0 to query.CountMS2Masses - 1 do
                
                            let mzRange = query.Ms2Masses.[ms2MassIndex]
                            let p = mzRangeSelector(pa, mzRange)

                            if getLockMz then
                                profile.[ms2MassIndex, specIdx] <- new Peak2D(p.Intensity, mzRange.LockValue, swathSpec.Rt)
                            else
                                profile.[ms2MassIndex, specIdx] <- new Peak2D(p.Intensity, p.Mz, swathSpec.Rt)

                    Some profile
                else
                    None
            )

    ///
    let getClosestMz (peaks: Peak1DArray, mzRange: RangeQuery) =
        let p1d = 
            peaks.Peaks
                .DefaultIfEmpty(new Peak1D(0., mzRange.LockValue))
                .ItemAtMin(fun x -> Math.Abs(x.Mz - mzRange.LockValue))
        if p1d.Mz > mzRange.HighValue || p1d.Mz < mzRange.LowValue then
            new Peak1D(0., mzRange.LockValue)
        else
            p1d

    ///
    let getRTProfiles (swathIndexer: SwathIndexer.SwathIndexer) (reader: IMzIODataReader) (*(spectrumSelector: seq<SwathIndexer.MSSwath> -> seq<SwathIndexer.MSSwath> list)*) (swathQuery: SwathQuery) =
        swathIndexer.GetRTProfilesFirstWnd(reader, swathQuery, false, (fun x -> [x.Take(1)]) ,getClosestMz)
    
    ///
    let getPlotFilePathFilePath outputDir (swathFileName:string) (plotName:string) (libraryFileName:string) =
        let swathFileName = swathFileName + "_" + (libraryFileName) + "_" + plotName 
        Path.Combine [|outputDir;swathFileName|]
    
    ///
    let getBinIdx width scantime = int ((scantime / width))    

    ///
    let identifyFeatures getPlotFilePathFilePath getRTProfiles (processParams:Domain.ConsensusSpectralLibraryParams) (libraryToAlign:Library) = 
        let library = libraryToAlign.Data
        let libBinned = 
            library
            |> Array.map (fun pepIons -> {pepIons with Fragments = pepIons.Fragments |> List.filter (fun x -> x.Number > processParams.MinFragmentLadderIdx) })
            |> Array.filter (fun pepIons -> pepIons.Fragments.Length > processParams.MinFragmentCount)
            |> Array.filter (fun pepIons -> pepIons.StringSequence |> String.filter Char.IsUpper |> String.length >= processParams.MinPeptideLength)
            |> Array.groupBy (fun ion -> getBinIdx processParams.BinningWindowWidth ion.ScanTime)
            |> Array.sortBy fst
        let filtered =
            libBinned 
            |> Array.map (fun (binIdx,pepIons) -> 
                let itemsSorted = pepIons |> Array.sortByDescending (fun x -> x.MeasuredApex)
                let frac = ((itemsSorted.Length |> float) * processParams.FractionOfMostAbundandIonsPerBin) |> int
                binIdx,
                itemsSorted
                |> Array.take (Math.Min(itemsSorted.Length-1,frac))
                |> Array.sortBy (fun x -> x.ScanTime)
                ) 
        
        [
        Chart.Histogram(libBinned |> Array.map snd  |> Array.concat |> Array.map (fun x -> x.Fragments.Length))
        |> Chart.withTraceName "NumberOfFragments: after Filter"
        Chart.Histogram(filtered |> Array.map snd |> Array.concat |> Array.map (fun x -> x.Fragments.Length))
        |> Chart.withTraceName (sprintf "NumberOfFragments: top %f" processParams.FractionOfMostAbundandIonsPerBin)
        ]
        |> Chart.Combine
        |> Chart.Show

        [
        Chart.Column(libBinned |> Array.map (fun (x,y) -> x,y.Length))
        |> Chart.withTraceName "libBinned"
        Chart.Column(filtered |> Array.map (fun (x,y) -> x,y.Length))
        |> Chart.withTraceName (sprintf "libBinned_filtered: top %f" processParams.FractionOfMostAbundandIonsPerBin)
        ]
        |> Chart.Combine
        |> Chart.Show
        
        let estimatedScanTimesInSwathFile = 
            filtered
            |> Array.map snd
            |> Array.concat
            |> Array.choose (fun pepIon -> 
                    let fragmentMzs = 
                        pepIon.Fragments
                        |> List.map (fun f -> f.MeanFragMz,f.MeanRelativeIntensity_Frags)
                        |> Array.ofSeq
                        |> BioFSharp.Mz.PeakArray.zipMzInt
                    let fragmentVector = 
                        fragmentMzs 
                        //|> fun x -> BioFSharp.Mz.PeakArray.peaksToNearestUnitDaltonBinVector x 100. 2000.
                        |> SparsePeakArray'.peaksToNearestBinVector processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset (processParams.MS2ScanRange |> fst) (processParams.MS2ScanRange |> snd)
                    let rtQuery   = MzIO.Processing.Query.createRangeQuery pepIon.ScanTime processParams.RtWindowWidth
                    let mzQueries = 
                        fragmentMzs 
                        |> Array.map (fun pk -> MzIO.Processing.Query.createRangeQuery pk.Mz 0.05)
                    
                    let swathQuery = MzIO.Processing.Query.createSwathQuery pepIon.PrecursorMZ rtQuery mzQueries
                    let profiles :Peak2D [,] option list  = getRTProfiles swathQuery
                    let correlationTrace = 
                        match profiles with 
                        | (Some queryRes)::t -> 
                            let fragVecs = 
                                queryRes
                                //Out: mz inner: intensitäten zu rts 
                                |> JaggedArray.ofArray2D
                                |> JaggedArray.transpose
                                |> Array.map (fun rtGroup ->
                                    rtGroup.[0].Rt, 
                                    rtGroup
                                    |> Array.map (fun p -> p.Mz,p.Intensity)
                                    |> Array.ofSeq
                                    |> BioFSharp.Mz.PeakArray.zipMzInt
                                    |> SparsePeakArray'.peaksToNearestBinVector processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset (processParams.MS2ScanRange |> fst) (processParams.MS2ScanRange |> snd)
                    
                                    )
                                |> Array.sortBy fst
                            let correlationTrace = 
                                fragVecs
                                |> Array.map (fun (rt, measuredFrags) -> 
                                    rt, SparsePeakArray'.dot fragmentVector measuredFrags
                                    )
                            [
                            Chart.Point correlationTrace
                            Chart.Point [pepIon.ScanTime,correlationTrace |> Array.maxBy snd |> snd]
                            ]
                            |> Chart.Combine
                            |> Chart.SaveHtmlAs (getPlotFilePathFilePath ((pepIon.StringSequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + pepIon.GlobalMod.ToString() + "Ch" + pepIon.Charge.ToString()) (libraryToAlign.FileName) ) 
                            
                            let xData,yData = correlationTrace |> Array.unzip
                            // TODO: refactor parameters to parameter type
                            try
                            let s = FSharpStats'.Wavelet.identify FSharpStats'.Wavelet.p xData yData
                            if Array.isEmpty s then 
                                None
                            elif s.Length = 1 then
                                toPeptideForLearning (s.[0].Apex.XVal |> float32 |> Some) pepIon
                                |> Some
                                
                            else
                                let sorted = s |> Array.sortByDescending (fun x -> x.Apex.YVal)
                                let mainPeak,sndPeak = sorted.[0],sorted.[1]
                                let ratio = sndPeak.Apex.YVal / mainPeak.Apex.YVal
                                if ratio > processParams.MaxRatioMostAbundandVsSecondAbundandPeak then None else
                                toPeptideForLearning (mainPeak.Apex.XVal |> float32 |> Some) pepIon 
                                |> Some
                            with 
                            | _ -> None
                        | _ -> None 
                    //let peaks = 
                    correlationTrace
                            
                            
                )
        estimatedScanTimesInSwathFile

    let plotIdentifiedFeatures getPlotFilePathFilePath libraryFiles swathFileFeatures =         
        Array.map2 (fun (lib:Library) (d:PeptideForLearning []) ->
            d
            |> Array.map (fun x -> x.SourceScanTime ,x.TargetScanTime)
            |> Chart.Point
            |> Chart.withX_AxisStyle "Source ScanTime"
            |> Chart.withY_AxisStyle "Target ScanTime"
            |> Chart.withTraceName lib.FileName
            ) libraryFiles swathFileFeatures 
        |> Chart.Combine
        |> Chart.withSize(1200.,1200.)
        |> Chart.SaveHtmlAs (getPlotFilePathFilePath "identifiedPepFeatures" "" )

    let plotDifferencesBeforeAlignment getPlotFilePathFilePath libraryFiles swathFileFeatures =         
        Array.map2 (fun (lib:Library) (d:PeptideForLearning []) ->
            d
            |> Array.map (fun x -> x.SourceScanTime - x.TargetScanTime)
            |> Chart.BoxPlot
            |> Chart.withX_AxisStyle "Differences before Alignment"
            |> Chart.withTraceName lib.FileName
            ) libraryFiles swathFileFeatures 
        |> Chart.Combine
        |> Chart.withSize(1200.,1200.)
        |> Chart.SaveHtmlAs (getPlotFilePathFilePath "differencesBeforeAlignment" "" )


    ///
    let downcastPipeline (x : IEstimator<_>) = 
        match x with 
        | :? IEstimator<ITransformer> as y -> y
        | _ -> failwith "downcastPipeline: expecting a IEstimator<ITransformer>"
       
    ///
    let createAlignmentResult (peptideIon:PeptideIon) (scanTimePrediction:ScanTimePrediction) :PeptideIon = 
        {
        StringSequence                              = peptideIon.StringSequence       
        GlobalMod                                   = peptideIon.GlobalMod            
        Charge                                      = peptideIon.Charge               
        PepSequenceID                               = peptideIon.PepSequenceID        
        ModSequenceID                               = peptideIon.ModSequenceID        
        PrecursorMZ                                 = peptideIon.PrecursorMZ          
        MeasuredMass                                = peptideIon.MeasuredMass         
        TheoMass                                    = peptideIon.TheoMass             
        AbsDeltaMass                                = peptideIon.AbsDeltaMass         
        MeanPercolatorScore                         = peptideIon.MeanPercolatorScore  
        QValue                                      = peptideIon.QValue               
        PEPValue                                    = peptideIon.PEPValue             
        ProteinNames                                = peptideIon.ProteinNames         
        Quant                                       = peptideIon.Quant                
        RelativeQuant                               = peptideIon.RelativeQuant        
        MeasuredApex                                = peptideIon.MeasuredApex         
        RelativeMeasuredApex                        = peptideIon.RelativeMeasuredApex 
        Seo                                         = peptideIon.Seo                  
        ScanTime                                    = float scanTimePrediction.TargetScanTime             
        ElutionWidth                                = peptideIon.ElutionWidth         
        Fragments                                   = peptideIon.Fragments            
        }    
        

    ///
    let initAlignFastTree (ctx:MLContext) (pepsForLearning: PeptideForLearning []) = 
        let data = ctx.Data.LoadFromEnumerable(pepsForLearning)
        let split = ctx.Data.TrainTestSplit(data, testFraction= 0.1)
        let trainer = ctx.Regression.Trainers.Gam(featureColumnName="Features",labelColumnName="TargetScanTime")
        let pipeline =    
            (ctx.Transforms.Concatenate("Features","SourceScanTime")|> downcastPipeline)
              .Append(trainer)
        let model = pipeline.Fit(split.TrainSet)    
        let metrics = 
            let tmp = ctx.Regression.Evaluate(model.Transform(split.TestSet),labelColumnName="TargetScanTime")
            {RSquared = tmp.RSquared}
        let predF = ctx.Model.CreatePredictionEngine<PeptideForLearning,ScanTimePrediction>(model)
        let predict pepIon = 
            pepIon
            |> toPeptideForLearning None
            |> predF.Predict
            |> createAlignmentResult pepIon
        metrics, predict

    ///
    let initAlignSpline (processParams:Domain.ConsensusSpectralLibraryParams) (pepsForLearning: PeptideForLearning []) = 
        let knotX,train,test = 
            let knotX,train,test = 
                pepsForLearning
                |> Array.groupBy (fun ion -> getBinIdx processParams.BinningWindowWidth (float ion.SourceScanTime))
                |> Array.map (fun (binIdx,ions) -> 
                    let dataS = ions |> Array.shuffleFisherYates
                    let knotX = 
                        dataS 
                        |> Array.map (fun x -> x.SourceScanTime)
                        |> Array.max
                    let train,test =
                        let nTest = 
                            (float dataS.Length) * 0.9
                            |> int
                        dataS.[.. nTest], dataS.[nTest+1 ..] 
                    float knotX, train, test
                    )
                |> Array.unzip3
            knotX |> Array.sort, train |> Array.concat |> Array.sortBy (fun x -> x.SourceScanTime), test |> Array.concat |> Array.sortBy (fun x -> x.SourceScanTime)
        let trainer lambda = 
            let train' = 
                train 
                |> Array.map (fun x -> float x.SourceScanTime, float x.TargetScanTime)
            let test' = 
                test 
                |> Array.map (fun x -> float x.SourceScanTime, float x.TargetScanTime)
            let fit = FSharp.Stats.Fitting.Spline.smoothingSpline train' (knotX) lambda 
            let rSquared = 
                let x,y,yHat = 
                    test'
                    |> Array.map (fun (x,y) -> 
                        x, y, fit x
                        )
                    |> Array.unzip3
                let rs = FSharp.Stats.Fitting.GoodnessOfFit.calculateDeterminationFromValue yHat y
                [
                Chart.Point(x,y)
                Chart.Line(x,yHat)
                ]
                |> Chart.Combine
                |> Chart.withTraceName (sprintf "rs: %f, l:%f" rs lambda)
                |> Chart.Show
                rs
            rSquared, fit     
        let rSquared,model = 
            [|0.01.. 0.05 .. 0.5|]
            |> Array.map trainer
            |> Array.maxBy fst
        let metrics = {RSquared = rSquared}          
        let predict pepIon = 
            (pepIon|> toPeptideForLearning None).SourceScanTime
            |> float 
            |> model
            |> fun x -> {TargetScanTime = float32 x}
            |> createAlignmentResult pepIon
        metrics, predict


    let plotAlignments getPlotFilePathFilePath alignedLibs =         
        alignedLibs
        |> Array.map (fun (alignedLib:AlignedLib) ->
            let source = 
                alignedLib.IdentifiedPeptideFeatures
                |> Array.map (fun x -> (x.ModSequenceID,x.GlobalMod,x.Charge),x)
                |> Map.ofSeq
            [
            alignedLib.IdentifiedPeptideFeatures
            |> Array.map (fun x -> x.SourceScanTime ,x.TargetScanTime)
            |> Chart.Point
            alignedLib.AlignedPeptideIons
            |> Array.choose (fun x -> 
                    match Map.tryFind (x.ModSequenceID,x.GlobalMod,x.Charge) source with 
                    | Some identPep -> 
                           (identPep.SourceScanTime,x.ScanTime)
                           |> Some
                    | None -> None
                )
            |> Array.sortBy fst
            |> Chart.Line
            ]
            |> Chart.Combine
            |> Chart.withTraceName (sprintf "%s, Rsquared:%f" alignedLib.FileName alignedLib.RegressionMetric.RSquared)
            ) 
        |> Chart.Combine
        |> Chart.withX_AxisStyle "Differences after Alignment"
        |> Chart.withSize(1200.,1200.)
        |> Chart.SaveHtmlAs (getPlotFilePathFilePath "Alignments" "" )

    let plotDifferencesAfterAlignment getPlotFilePathFilePath alignedLibs =         
        alignedLibs
        |> Array.map (fun (alignedLib:AlignedLib) ->
            let source = 
                alignedLib.IdentifiedPeptideFeatures
                |> Array.map (fun x -> (x.ModSequenceID,x.GlobalMod,x.Charge),x)
                |> Map.ofSeq
            alignedLib.AlignedPeptideIons
            |> Array.choose (fun x -> 
                    match Map.tryFind (x.ModSequenceID,x.GlobalMod,x.Charge) source with 
                    | Some identPep -> 
                            (float identPep.TargetScanTime - x.ScanTime)
                            |> Some
                    | None -> None
                )
            |> Chart.BoxPlot
            |> Chart.withTraceName (sprintf "%s, Rsquared:%f" alignedLib.FileName alignedLib.RegressionMetric.RSquared)
            |> Chart.withX_AxisStyle "Source ScanTime"
            ) 
        |> Chart.Combine
        |> Chart.withSize(1200.,1200.)
        |> Chart.SaveHtmlAs (getPlotFilePathFilePath "differencesAfterAlignment" "" )
    /////
    //let performAlignment outDir align (source: LibToAlign) =
    
    //    ///
    //    let peptidesForLearning = source.IdentifiedPeptideFeatures
 
    //    ///
    //    let metrics,model: RegressionMetrics*(PeptideIon->AlignmentResult) = 
    //        align peptidesForLearning
        
    //    metrics,model

    ///
    let createSpectralLibrary (logger:NLog.Logger) (processParams:Domain.ConsensusSpectralLibraryParams) (outputDir:string) (libraryFiles:string) (targetSwathFile:string) = 
        logger.Trace (sprintf "Input directory containing library files: %s" libraryFiles)
        logger.Trace (sprintf "Input Swath file to align to: %s" targetSwathFile)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension targetSwathFile) + ".csl"
            Path.Combine [|outputDir;fileName|]
        
        logger.Trace "Init align function"
        let align = 
            match processParams.ConsensusAlignmentAlgorithm with 
            | Domain.ConsensusAlignmentAlgorithm.FastTree ->
                let ctx = new ML.MLContext()
                let align = initAlignFastTree ctx
                align
            | Domain.ConsensusAlignmentAlgorithm.Spline -> 
                let align = initAlignSpline processParams
                align
        logger.Trace "Init align function: finished"
         
        logger.Trace "Reading and preparing .sl files for alignment"
        let libraryFilePaths = System.IO.Directory.GetFiles (libraryFiles, "*.sl")
        let libraryFiles = libraryFilePaths |> Array.map readLibraryFrom 
        logger.Trace "Reading and preparing .sl files for alignment: finished"
        
        logger.Trace "Init connection to swath mass spectrum data."
        let inReader = Core.MzIO.Reader.getReader targetSwathFile :?> MzIO.MzSQL.MzSQL
        inReader.Connection.Open()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader
        let inTr = inReader.BeginTransaction()
        logger.Trace "Init connection to swath mass spectrum data.:finished"

        logger.Trace "Create Swath RetentionTime index"
        let retTimeIdxed = Query.getSwathIdx inReader inRunID
        let getRTProfiles = getRTProfiles retTimeIdxed inReader
        logger.Trace "Create Swath RetentionTime index:finished"

        logger.Trace "Identify features in Swath file"
        let swathFileName = Path.GetFileNameWithoutExtension targetSwathFile
        let getPlotFilePathFilePath = getPlotFilePathFilePath outputDir swathFileName                    
        let swathFileFeatures = 
            libraryFiles 
            |> Array.map (fun libraryToAlign -> identifyFeatures getPlotFilePathFilePath getRTProfiles processParams libraryToAlign)
        plotIdentifiedFeatures getPlotFilePathFilePath libraryFiles swathFileFeatures
        plotDifferencesBeforeAlignment getPlotFilePathFilePath libraryFiles swathFileFeatures
        logger.Trace "Identify features in Swath file: finished"
        
        logger.Trace "Pruning peptide features across libaries"
        let prunedSwathfileFeatures = swathFileFeatures
        logger.Trace "Pruning peptide features across libaries: finished"
        logger.Trace "Performing alignments"
        let libsToAlign = Array.map2 createLibToAlign libraryFiles swathFileFeatures
        let alignedLibs = 
            libsToAlign
            |> Array.map (fun libtoAlign ->
                    let regMetrics,alignFunc = align libtoAlign.IdentifiedPeptideFeatures
                    let alignmentResult = libtoAlign.PeptideIons |> Array.map alignFunc
                    createAlignedLib libtoAlign regMetrics alignmentResult 
                )
            |> Array.sortByDescending (fun (alib) -> alib.RegressionMetric.RSquared)
        plotAlignments getPlotFilePathFilePath alignedLibs
        plotDifferencesAfterAlignment getPlotFilePathFilePath alignedLibs
        logger.Trace "Performing alignments:finished"
        logger.Trace "Transfering PeptideIons to Consensus library"
        let missingPeptideIons = 
            alignedLibs
            |> Array.map (fun x -> 
                x.AlignedPeptideIons 
                |> Array.map (fun x -> createIonSpecies x.ModSequenceID x.Charge x.GlobalMod)
                )
            |> Array.concat
            |> Set.ofArray
        let consensusLibrary =
            let stillMissingPeptides,transferedIons = 
                alignedLibs
                |> Array.fold (fun ((missingPeps:Set<IonSpecies>),transferedPeps) alignedLib ->
                        logger.Trace (sprintf "%i peptide Ions still missing" missingPeps.Count)
                        let ionsToTransfer = 
                            alignedLib.AlignedPeptideIons
                            |> Array.filter (fun x -> missingPeps |> Set.contains (createIonSpecies x.ModSequenceID x.Charge x.GlobalMod))
                        missingPeps-(ionsToTransfer |> Array.map (fun x -> createIonSpecies x.ModSequenceID x.Charge x.GlobalMod)|> Set.ofArray),
                        ionsToTransfer::transferedPeps
                        ) (missingPeptideIons,[])
            logger.Trace (sprintf "%i peptide Ions still missing" stillMissingPeptides.Count)
            transferedIons
            |> Array.concat
            |> Array.sortBy (fun x -> x.ScanTime)
        logger.Trace "Transfering PeptideIons to Consensus library:finished"
        Json.serializeAndWrite outFilePath consensusLibrary
        logger.Trace "Creating Consensus library: finished"
        
            


        
        
            
        
