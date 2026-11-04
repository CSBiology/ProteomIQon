namespace ProteomIQon

open System.IO
open System.Data.SQLite
open ProteomIQon.Core
open Core.MzIO
open Dto
open FSharp.Stats
open BioFSharp.Mz.Quantification
open BioFSharp.Mz
open FSharpAux.IO.SchemaReader
open Plotly.NET
open BioFSharp
open MzIO.Processing
open BioFSharp.Mz.SearchDB
open PeakQuantification
open PSMQuantification

module PSMBasedQuantification =
    module Query = 
        open System
        open System.Collections.Generic
        open System.Linq
        open MzIO.Commons.Arrays
        open MzIO.Processing.MzIOLinq
        open MzIO.Binary
         
        /// Extract a rt profile for specified target mass and rt range.
        /// Mz range peak aggregation is closest lock mz.
        /// Profile array with index corresponding to continous mass spectra over rt range and mz range given.
        let initRTProfile (readspecPeaks:string -> Peak1DArray)  (rtIndex: IMzIOArray<RtIndexEntry>) (rtRange: RangeQuery) (mzRange: RangeQuery) =
            let entries = RtIndexEntry.Search(rtIndex, rtRange).ToArray()
            //printfn "RtProfile %i" entries.Length
            let profile = Array.zeroCreate<Peak2D> entries.Length
            for rtIdx = 0 to entries.Length-1 do
                let entry = entries.[rtIdx]
                let peaks = (readspecPeaks entry.SpectrumID).Peaks
                let p = (RtIndexEntry.MzSearch (peaks, mzRange)).DefaultIfEmpty(Peak1D(0., mzRange.LockValue))
                        |> fun x -> RtIndexEntry.ClosestMz (x, mzRange.LockValue)
                        |> fun x -> RtIndexEntry.AsPeak2D (x, entry.Rt)
                profile.[rtIdx] <- p
            profile

        /// The retention time index and the scan times of the MS1 spectra, both sorted by scan time, from one pass
        /// over the spectrum descriptions with the streaming parser (SpectrumIndex.readMs1Index).
        let readMs1Index (reader: MzIO.MzSQL.MzSQL) (runID: string) (warn: string -> unit) =
            let ms1SortedByScanTime = SpectrumIndex.readMs1Index reader runID warn ignore
            let rtIndex =
                ms1SortedByScanTime
                |> Array.map (fun (scanTime, id) -> RtIndexEntry(scanTime, id))
                |> MzIOArray.ToMzIOArray
            rtIndex, ms1SortedByScanTime


    ///
    let initSpline binningWindowWidth (scanTimeVsDiff: (float*float) []) = 
        let getBinIdx width scantime = int ((scantime / width))    
        let knotX,train,test = 
            let knotX,train,test = 
                scanTimeVsDiff
                |> Array.groupBy (fun (s,d) -> getBinIdx binningWindowWidth (s))
                |> Array.map (fun (binIdx,ions) -> 
                    let dataS = ions |> Array.shuffleFisherYates
                    let knotX = 
                        dataS 
                        |> Array.map fst
                        |> Array.max
                    let train,test =
                        let nTest = 
                            (float dataS.Length) * 0.9
                            |> int
                        dataS.[.. nTest], dataS.[nTest+1 ..] 
                    float knotX, train, test
                    )
                |> Array.unzip3
            knotX |> Array.sort, train |> Array.concat |> Array.sortBy fst, test |> Array.concat |> Array.sortBy fst
        let trainer lambda = 
            let train' = train 
            let test' = test 
            let fit = FSharp.Stats.Fitting.Spline.smoothingSpline train' (knotX) lambda 
            let rSquared = 
                let x,y,yHat = 
                    test'
                    |> Array.map (fun (x,y) -> 
                        x, y, fit x
                        )
                    |> Array.unzip3
                let rs = FSharp.Stats.Fitting.GoodnessOfFit.calculateDeterminationFromValue yHat y
                rs
            rSquared, fit     
        let rSquared,model = 
            [|0.01.. 0.05 .. 0.5|]
            |> Array.map trainer
            |> Array.maxBy fst          
        rSquared, model

    ///
    let getSpec (reader:MzIO.IO.IMzIODataReader) (ms1ID: string)  =
        Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(ms1ID).Peaks)
        |> fun (mzData,intensityData) -> PeakArray.zip mzData intensityData
    
    ///
    let lightQualityFilter lowerBorder upperBorder (quantResults:QuantificationResult[]) =
        let medianApexIntensities = 
            quantResults
            |> Array.map (fun x -> x.MeasuredApex_Light)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        let medianQuantIntensities = 
            quantResults
            |> Array.map (fun x -> x.Quant_Light)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        quantResults
        |> Array.filter (fun x -> 
            let qualR = 
                let apexNorm = x.MeasuredApex_Light / medianApexIntensities
                let quantNorm =  x.Quant_Light / medianQuantIntensities
                quantNorm / apexNorm
                |> FSharp.Stats.Ops.log2
            (qualR > lowerBorder && qualR < upperBorder) || (nan.Equals(qualR) )
            ) 

    ///
    let heavyQualityFilter lowerBorder upperBorder (quantResults:QuantificationResult[]) =
        let medianApexIntensities = 
            quantResults
            |> Array.map (fun x -> x.MeasuredApex_Heavy)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        let medianQuantIntensities = 
            quantResults
            |> Array.map (fun x -> x.Quant_Heavy)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        quantResults
        |> Array.filter (fun x -> 
            let qualR = 
                let apexNorm = x.MeasuredApex_Heavy / medianApexIntensities
                let quantNorm =  x.Quant_Heavy / medianQuantIntensities
                quantNorm / apexNorm
                |> FSharp.Stats.Ops.log2
            (qualR > lowerBorder && qualR < upperBorder) || (nan.Equals(qualR) )
            )
    
    ///
    let substractBaseLine (logger: NLog.Logger) (baseLineParams:Domain.BaseLineCorrection) (yData:float []) =
        if yData.Length > 500 then
            yData
        else
            let baseLine = BaselineAls.baselineAls baseLineParams.MaxIterations baseLineParams.Lambda baseLineParams.P yData
            Array.map2 (fun y b ->
                           let c = y - b
                           if c < 0. then 0. else c
                       ) yData baseLine

    ///
    let initGetProcessedXIC logger (baseLineCorrection:Domain.BaseLineCorrection option) getPeaks idx scanTimeWindow mzWindow_Da meanScanTime meanPrecMz =
        let rtQuery = Query.createRangeQuery meanScanTime scanTimeWindow
        let mzQuery = Query.createRangeQuery meanPrecMz mzWindow_Da
        let retData',itzData' =
            let tmp =
                getPeaks idx rtQuery mzQuery
                |> Array.map (fun (p:MzIO.Binary.Peak2D) -> p.Rt , p.Intensity)
            tmp
            |> Array.mapi (fun i (rt,intensity) ->
                            if i = 0 || i = tmp.Length-1 || intensity > 0. then
                                Some (rt,intensity)
                            else
                                let rt',intensity' = tmp.[i-1]
                                if intensity' = 0. then
                                    Some (rt,intensity)
                                elif intensity' > (100. * (intensity+1.)) then
                                    None
                                else
                                    Some (rt,intensity)
                          )
            |> Array.choose id
            |> Array.unzip
        match baseLineCorrection with
        | Some baseLineParams ->
            retData', substractBaseLine logger baseLineParams itzData', itzData'
        | None ->
            retData',itzData',itzData'

    type InferredQuantification = {
        Model                       :HULQ.PeakModel option
        Area                        :float
        StandardErrorOfPrediction   :float
        MeasuredApexIntensity       :float
        Correlation_Light_Heavy     :float
        SearchRTMinusFittedRT       :float
        ClusterComparison           :ClusterComparison               
        EstimatedParams             :float[]
        X_Xic                       :float[]
        Y_Xic                       :float[]
        Y_Xic_uncorrected           :float[]
        xPeak                       :float[]
        yFitted                     :float[]
        }

    ///
    let calcCorrelation (xValues:float []) (quantifiedPeak:HULQ.QuantifiedPeak) (inferredPeak:HULQ.QuantifiedPeak) = 
        let getValue (model:HULQ.PeakModel) estParams =
            match model with
            | HULQ.PeakModel.Gaussian m -> 
                m.GetFunctionValue (vector estParams)
            | HULQ.PeakModel.EMG m -> 
                m.GetFunctionValue (vector estParams)
        match quantifiedPeak.Model, inferredPeak.Model with 
        | Some q , Some i -> 
            let xValuesBW = 
                [|for i = 1 to xValues.Length-1 do abs(xValues.[i] - xValues.[i-1])|]
                |> Array.median
            let xValues = [|xValues.[0] .. xValuesBW/2. .. xValues.[xValues.Length-1]|]
            let fQ = getValue q quantifiedPeak.EstimatedParams
            let yQ = xValues |> Array.map fQ
            let fI = getValue i inferredPeak.EstimatedParams
            let yI = xValues |> Array.map fI
            FSharp.Stats.Correlation.Seq.pearson yQ yI
        | _ -> nan

    ///
    let quantifyPeptides diagCharts zipCharts (processParams:Domain.QuantificationParams) (outputDir:string) (cn:SQLiteConnection) (instrumentOutput:string) (scoredPSMs:string)  =
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension scoredPSMs)
        logger.Trace (sprintf "Input file: %s" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        logger.Trace (sprintf "Now performing Quantification using: %s and %s, Results will be written to: %s" instrumentOutput scoredPSMs outputDir)
        // initialize Reader and Transaction
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".quant"
            Path.Combine [|outputDir;fileName|]
        logger.Trace (sprintf "outFilePath:%s" outFilePath)
        logger.Trace "Read scored PSMs"
        ///
        let qpsms =
            Csv.CsvReader<PSMStatisticsResult>(SchemaMode=Csv.Fill).ReadFile(scoredPSMs,'\t',false,1)
            |> Array.ofSeq
        logger.Trace "Read scored PSMs:finished"
        if Array.isEmpty qpsms then
            logger.Warn (sprintf "%s contains no PSMs, %s gets no quantifications." scoredPSMs outFilePath)
            Seq.empty<QuantificationResult>
            |> SeqIO'.csv "\t" true false
            |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        else
        //
        let plotDirectory =
            let fileName = sprintf "%s_plots" (Path.GetFileNameWithoutExtension instrumentOutput)
            let path = Path.Combine [|outputDir;fileName|]
            if System.IO.Directory.Exists path then
                path
            else
                System.IO.Directory.CreateDirectory path |> ignore
                path
        logger.Trace (sprintf "plotDirectory:%s" plotDirectory)
        logger.Trace "Copy peptide DB into Memory"
        let memoryDB = SearchDB.copyDBIntoMemory cn
        logger.Trace "Copy peptide DB into Memory: finished"
        logger.Trace "Get peptide lookUp function"
        let dBParams     = getSDBParamsByCn memoryDB
        //let massLookUp = prepareSelectMassByModSequenceAndGlobalMod memoryDB
        let peptideLookUp = getThreadSafePeptideLookUpFromFileBySequenceAndGMod memoryDB dBParams
        let isotopePatterns = new IsotopePatterns.IsotopePatternLookUp(logger, memoryDB)
        let calcIonSeries aal = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal
        logger.Trace "Get peptide lookUp function: finished"
        // initialize Reader and Transaction
        logger.Trace "Init connection to mass spectrum data."
        let inReader = Core.MzIO.Reader.getReader instrumentOutput :?> MzIO.MzSQL.MzSQL
        inReader.Connection.Open()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader       
        let inTr = inReader.BeginTransaction()
        let countMatchedMasses (peptide: LookUpResult<AminoAcids.AminoAcid>)(psms: PSMStatisticsResult []) =
            let frag = 
                let ionSeries = (calcIonSeries peptide.BioSequence).TargetMasses
                [1. .. 2.]
                |> List.collect (fun ch -> 
                    ionSeries 
                    |> List.map (fun x -> Mass.toMZ x.MainPeak.Mass ch)
                )                  
            psms
            |> Array.map (fun psm -> 
                let spec = inReader.ReadSpectrumPeaks psm.PSMId
                let sum = 
                    spec.Peaks 
                    |> Seq.filter (fun peak -> 
                        frag
                        |> List.exists (fun ion -> abs (ion - peak.Mz) <= (Mass.deltaMassByPpm 100. peak.Mz))
                    )
                    |> Seq.sumBy (fun x -> x.Intensity)
                psm,sum
            )
        logger.Trace "Create RetentionTime index and read the ms1s"
        /// retTimeIdxed and ms1SortedByScanTime: the MS1s sorted by their scanTime
        let retTimeIdxed, ms1SortedByScanTime = Query.readMs1Index inReader inRunID (fun message -> logger.Warn message)
        logger.Trace (sprintf "Create RetentionTime index and read the ms1s:finished, %i MS1s" ms1SortedByScanTime.Length)
        
        logger.Trace "Estimate precursor mz standard deviation and mz correction."
        ///
        let ms1AccuracyEstimate,scanTimeToMzCorrection = 
            let scanTimeVsDelta = 
                qpsms
                |> Array.map (fun x -> 
                    let precMz = x.PrecursorMZ
                    let theMz  = Mass.toMZ x.TheoMass (float x.Charge)
                    let diff = precMz - theMz 
                    x.ScanTime,diff
                    )
            let borders =
                scanTimeVsDelta
                |> Seq.map snd
                |> Array.ofSeq
                |> FSharp.Stats.Testing.Outliers.tukey 3.
            let filteredValues =
                scanTimeVsDelta
                |> Array.ofSeq
                |> Array.filter (fun (s,d) -> (d <= borders.Upper && d >= borders.Lower) && d < 0.1)
            let scanTimeToMzCorrection =
                let runTime = scanTimeVsDelta |> Array.maxBy fst |> fst
                if runTime > 20. && qpsms.Length > 500 then
                    let binWidth =
                        System.Math.Min(runTime / 2., 20.)
                    let stabw = filteredValues |> Seq.stDevBy snd
                    try
                        let r,f = initSpline binWidth filteredValues
                        f
                    with ex ->
                        // the smoothing spline needs enough filtered points per bin; sparse samples fall back to the median correction below
                        logger.Trace (sprintf "Mz correction spline failed (%s), falling back to median correction." ex.Message)
                        let m =
                            filteredValues
                            |> Seq.map snd
                            |> Seq.median
                        fun scanTime -> m
                else
                    let m = 
                        filteredValues 
                        |> Seq.map snd 
                        |> Seq.median 
                    fun scanTime -> m
            let stDev = 
                filteredValues 
                |> Seq.map snd 
                |> Seq.stDev
            if diagCharts then 
                [
                Chart.Point(scanTimeVsDelta)
                |> Chart.withTraceInfo "Raw"
                Chart.Line(scanTimeVsDelta |> Array.sortBy fst |> Array.map (fun (st,d) -> st, scanTimeToMzCorrection st))
                ]
                |> Chart.combine
                |> Chart.saveHtml(Path.Combine[|plotDirectory; "mzErrorAndCorrection"|])
            stDev, scanTimeToMzCorrection 
        logger.Trace (sprintf "Estimate precursor mz standard deviation and mz correction.:finished, standard deviation: %f" ms1AccuracyEstimate) 

        logger.Trace "init lookup functions"        
        ///
        let comparePredictedAndMeasuredIsotopicCluster = initComparePredictedAndMeasuredIsotopicCluster isotopePatterns (fun ms1ID lowMz highMz -> getSpec inReader ms1ID |> Array.filter (fun x -> x.Mz < highMz && x.Mz > lowMz)) (SpectrumIndex.initClosestMs1 ms1SortedByScanTime) ms1AccuracyEstimate     
        
        let mzWindow = 
            match processParams.XicExtraction.MzWindow_Da with 
            | Domain.Window.Fixed v -> v 
            | Domain.Window.Estimate -> 
                let mzW = ms1AccuracyEstimate*4.
                logger.Trace (sprintf "optimal mz Window for XIC look up found by estimation :%f Da" mzW)  
                mzW

        let readSpecPeaksWithMem = FSharpAux.Memoization.memoize inReader.ReadSpectrumPeaks
        ///
        let getPeaks = Query.initRTProfile readSpecPeaksWithMem
        let getXIC = initGetProcessedXIC logger processParams.BaseLineCorrection getPeaks retTimeIdxed processParams.XicExtraction.ScanTimeWindow mzWindow    
        // PSMBasedQuantification has no ion mobility.
        let getXICWithIM scanTime mz (_: float) = getXIC scanTime mz
        
        match processParams.XicExtraction.XicProcessing with
        | Domain.XicProcessing.Gabor3D _ -> failwith "Gabor3D peak detection is only available in PSMBasedQuantificationTIMs."
        | _ -> ()
        ///
        let identifyPeaks = initIdentifyPeaks processParams.XicExtraction.XicProcessing
        logger.Trace "init lookup functions:finished"
        
        logger.Trace "init quantification functions"
        ///
        let labledQuantification (pepIon:PeptideIon) (psms:PSMStatisticsResult []) = 
            try
            let bestQValue,bestPepValue,prots = psms |> Array.minBy (fun x -> x.QValue) |> fun x -> x.QValue, x.PEPValue,x.ProteinNames
            let unlabledPeptide = peptideLookUp pepIon.Sequence 0
            let labeledPeptide  = peptideLookUp pepIon.Sequence 1
            let targetPeptide = if pepIon.GlobalMod = 0 then unlabledPeptide else labeledPeptide            
            let psmsWithMatchedSums = countMatchedMasses targetPeptide psms 
            let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
            let theoMz = Mass.toMZ targetPeptide.Mass (float pepIon.Charge)
            let averagePSM = average (fun _ -> nan) getXICWithIM scanTimeToMzCorrection theoMz psmsWithMatchedSums
            let avgMass = Mass.ofMZ (averagePSM.MeanPrecMz) (pepIon.Charge |> float)
            let peaks = identifyPeaks averagePSM.X_Xic averagePSM.Y_Xic
            if Array.isEmpty peaks then 
                if diagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "noPeaks" plotDirectory
                logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                None
            else
            let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks averagePSM.WeightedAvgScanTime
            let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify                        
            let searchRTMinusFittedRT = searchRTMinusFittedRtTarget averagePSM.WeightedAvgScanTime quantP
            if quantP.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantP.EstimatedParams then 
                if diagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "fittingFailed" plotDirectory    
                logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                None
            else
            let inferredScanTime = chooseScanTime processParams.XicExtraction.ScanTimeWindow searchRTMinusFittedRT averagePSM.WeightedAvgScanTime quantP 
            let clusterComparisonTarget = comparePredictedAndMeasuredIsotopicCluster averagePSM.X_Xic averagePSM.Y_Xic averagePSM.Y_Xic_uncorrected pepIon.Charge targetPeptide quantP.EstimatedParams.[1] averagePSM.MeanPrecMz            
            if pepIon.GlobalMod = 0 then
                let mzHeavy = 
                    let mz = Mass.toMZ (labeledPeptide.Mass) (pepIon.Charge|> float)
                    let correctedMz = scanTimeToMzCorrection inferredScanTime + mz
                    correctedMz 
                let inferredQuant = 
                    let inferredXicHeavy = getInferredXic getXICWithIM averagePSM.WeightedAvgScanTime mzHeavy nan    
                    let inferredPeaksHeavy = identifyPeaks inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic
                    if Array.isEmpty inferredPeaksHeavy then 
                        if diagCharts then saveErrorChart inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic pepIon "noInferredPeaks" plotDirectory
                        logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    else
                    let peakToQuantifyHeavy = BioFSharp.Mz.Quantification.HULQ.getPeakBy inferredPeaksHeavy inferredScanTime
                    let quantPHeavy         = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantifyHeavy          
                    if quantPHeavy.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantPHeavy.EstimatedParams then 
                        if diagCharts then saveErrorChart inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic pepIon "fittingInferredFailed" plotDirectory    
                        logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    else
                    // let inferred_Heavy = quantifyInferredPeak getXIC identifyPeaks mz_Heavy averagePSM.WeightedAvgScanTime inferredScanTime
                    let searchRTMinusFittedRTHeavy = searchRTMinusFittedRtTarget inferredScanTime quantPHeavy
                    let clusterComparisonHeavy = comparePredictedAndMeasuredIsotopicCluster inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic inferredXicHeavy.Y_Xic_uncorrected pepIon.Charge labeledPeptide quantP.EstimatedParams.[1] mzHeavy
                    let corrLightHeavy  = calcCorrelation averagePSM.X_Xic quantP quantPHeavy
                    {
                    Model                       = quantPHeavy.Model
                    Area                        = quantPHeavy.Area
                    StandardErrorOfPrediction   = quantPHeavy.StandardErrorOfPrediction
                    MeasuredApexIntensity       = quantPHeavy.MeasuredApexIntensity
                    Correlation_Light_Heavy     = corrLightHeavy
                    SearchRTMinusFittedRT       = searchRTMinusFittedRTHeavy 
                    ClusterComparison           = clusterComparisonHeavy     
                    EstimatedParams             = quantPHeavy.EstimatedParams
                    X_Xic                       = inferredXicHeavy.X_Xic
                    Y_Xic                       = inferredXicHeavy.Y_Xic 
                    Y_Xic_uncorrected           = inferredXicHeavy.Y_Xic_uncorrected 
                    xPeak                       = peakToQuantifyHeavy.XData 
                    yFitted                     = quantPHeavy.YPredicted
                    }
                    |> Some                      
                match inferredQuant with 
                | Some successfulQuant -> 
                    if diagCharts then 
                        saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                        peakToQuantify.XData peakToQuantify.YData quantP.YPredicted successfulQuant.X_Xic successfulQuant.Y_Xic successfulQuant.xPeak successfulQuant.yFitted peaks clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = pepIon.Sequence
                    GlobalMod                                   = pepIon.GlobalMod
                    Charge                                      = pepIon.Charge
                    PepSequenceID                               = pepIon.PepSequenceID
                    ModSequenceID                               = pepIon.ModSequenceID
                    PrecursorMZ                                 = averagePSM.MeanPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = unlabledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-unlabledPeptide.Mass)
                    MeanPercolatorScore                         = averagePSM.MeanScore
                    QValue                                      = bestQValue
                    PEPValue                                    = bestPepValue
                    ProteinNames                                = prots
                    QuantMz_Light                               = averagePSM.MeanPrecMz
                    Quant_Light                                 = quantP.Area
                    MeasuredApex_Light                          = quantP.MeasuredApexIntensity
                    Seo_Light                                   = quantP.StandardErrorOfPrediction
                    Params_Light                                = quantP.EstimatedParams            
                    Difference_SearchRT_FittedRT_Light          = searchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Light            = clusterComparisonTarget.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Light   = clusterComparisonTarget.KLDiv_Corrected
                    QuantMz_Heavy                               = mzHeavy
                    Quant_Heavy                                 = successfulQuant.Area
                    MeasuredApex_Heavy                          = successfulQuant.MeasuredApexIntensity
                    Seo_Heavy                                   = successfulQuant.StandardErrorOfPrediction
                    Params_Heavy                                = successfulQuant.EstimatedParams 
                    Difference_SearchRT_FittedRT_Heavy          = successfulQuant.SearchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Heavy            = successfulQuant.ClusterComparison.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Heavy   = successfulQuant.ClusterComparison.KLDiv_Corrected
                    Correlation_Light_Heavy                     = successfulQuant.Correlation_Light_Heavy
                    QuantificationSource                        = QuantificationSource.PSM
                    IsotopicPatternMz_Light                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Light     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Light    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Light                               = averagePSM.X_Xic 
                    IntensityTrace_Observed_Light               = averagePSM.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Light              = averagePSM.Y_Xic
                    IsotopicPatternMz_Heavy                     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Heavy     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Heavy    = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Heavy                               = successfulQuant.X_Xic 
                    IntensityTrace_Observed_Heavy               = successfulQuant.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Heavy              = successfulQuant.Y_Xic
                    AlignmentScore                              = nan
                    AlignmentQValue                             = 0.
                    IonMobility                                 = nan
                    }
                    |> Option.Some
                | None -> 
                    if diagCharts then 
                        saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                        peakToQuantify.XData peakToQuantify.YData quantP.YPredicted [||] [||] [||] [||] peaks clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = pepIon.Sequence
                    GlobalMod                                   = pepIon.GlobalMod
                    Charge                                      = pepIon.Charge
                    PepSequenceID                               = pepIon.PepSequenceID
                    ModSequenceID                               = pepIon.ModSequenceID
                    PrecursorMZ                                 = averagePSM.MeanPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = unlabledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-unlabledPeptide.Mass)
                    MeanPercolatorScore                         = averagePSM.MeanScore
                    QValue                                      = bestQValue
                    PEPValue                                    = bestPepValue
                    ProteinNames                                = prots
                    QuantMz_Light                               = averagePSM.MeanPrecMz
                    Quant_Light                                 = quantP.Area
                    MeasuredApex_Light                          = quantP.MeasuredApexIntensity
                    Seo_Light                                   = quantP.StandardErrorOfPrediction
                    Params_Light                                = quantP.EstimatedParams            
                    Difference_SearchRT_FittedRT_Light          = searchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Light            = clusterComparisonTarget.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Light   = clusterComparisonTarget.KLDiv_Corrected
                    QuantMz_Heavy                               = mzHeavy
                    Quant_Heavy                                 = nan
                    MeasuredApex_Heavy                          = nan
                    Seo_Heavy                                   = nan
                    Params_Heavy                                = [||]
                    Difference_SearchRT_FittedRT_Heavy          = nan
                    KLDiv_Observed_Theoretical_Heavy            = nan
                    KLDiv_CorrectedObserved_Theoretical_Heavy   = nan
                    Correlation_Light_Heavy                     = nan
                    QuantificationSource                        = QuantificationSource.PSM
                    IsotopicPatternMz_Light                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Light     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Light    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Light                               = averagePSM.X_Xic 
                    IntensityTrace_Observed_Light               = averagePSM.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Light              = averagePSM.Y_Xic
                    IsotopicPatternMz_Heavy                     = [||]
                    IsotopicPatternIntensity_Observed_Heavy     = [||]
                    IsotopicPatternIntensity_Corrected_Heavy    = [||]
                    RtTrace_Heavy                               = [||]
                    IntensityTrace_Observed_Heavy               = [||]
                    IntensityTrace_Corrected_Heavy              = [||]
                    AlignmentScore                              = nan
                    AlignmentQValue                             = 0.
                    IonMobility                                 = nan
                    }
                    |> Option.Some
            else
                let mzLight = 
                    let mz = Mass.toMZ (unlabledPeptide.Mass) (pepIon.Charge|> float)
                    let correctedMz = scanTimeToMzCorrection inferredScanTime + mz
                    correctedMz   
                let inferredQuant = 
                    let inferredXicLight = getInferredXic getXICWithIM averagePSM.WeightedAvgScanTime mzLight nan    
                    let inferredPeaksLight = identifyPeaks inferredXicLight.X_Xic inferredXicLight.Y_Xic
                    if Array.isEmpty inferredPeaksLight then 
                        if diagCharts then saveErrorChart inferredXicLight.X_Xic inferredXicLight.Y_Xic pepIon "noInferredPeaks" plotDirectory
                        logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    else
                    let peakToQuantifyLight = BioFSharp.Mz.Quantification.HULQ.getPeakBy inferredPeaksLight inferredScanTime
                    let quantPLight         = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantifyLight          
                    if quantPLight.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantPLight.EstimatedParams then 
                        if diagCharts then saveErrorChart inferredXicLight.X_Xic inferredXicLight.Y_Xic pepIon "fittingInferredFailed" plotDirectory    
                        logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    else
                    // let inferred_Heavy = quantifyInferredPeak getXIC identifyPeaks mz_Heavy averagePSM.WeightedAvgScanTime inferredScanTime
                    let searchRTMinusFittedRTLight = searchRTMinusFittedRtTarget inferredScanTime quantPLight
                    let clusterComparisonLight = comparePredictedAndMeasuredIsotopicCluster inferredXicLight.X_Xic inferredXicLight.Y_Xic inferredXicLight.Y_Xic_uncorrected pepIon.Charge unlabledPeptide quantP.EstimatedParams.[1] mzLight
                    let corrLightHeavy  = calcCorrelation averagePSM.X_Xic quantP quantPLight
                    {
                    Model                       = quantPLight.Model
                    Area                        = quantPLight.Area
                    StandardErrorOfPrediction   = quantPLight.StandardErrorOfPrediction
                    MeasuredApexIntensity       = quantPLight.MeasuredApexIntensity
                    Correlation_Light_Heavy     = corrLightHeavy
                    SearchRTMinusFittedRT       = searchRTMinusFittedRTLight 
                    ClusterComparison           = clusterComparisonLight     
                    EstimatedParams             = quantPLight.EstimatedParams
                    X_Xic                       = inferredXicLight.X_Xic
                    Y_Xic                       = inferredXicLight.Y_Xic 
                    Y_Xic_uncorrected           = inferredXicLight.Y_Xic_uncorrected 
                    xPeak                       = peakToQuantifyLight.XData 
                    yFitted                     = quantPLight.YPredicted
                    }
                    |> Some                      
                match inferredQuant with 
                | Some successfulQuant -> 
                    if diagCharts then 
                        saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                        peakToQuantify.XData peakToQuantify.YData quantP.YPredicted successfulQuant.X_Xic successfulQuant.Y_Xic successfulQuant.xPeak successfulQuant.yFitted peaks clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = pepIon.Sequence
                    GlobalMod                                   = pepIon.GlobalMod
                    Charge                                      = pepIon.Charge
                    PepSequenceID                               = pepIon.PepSequenceID
                    ModSequenceID                               = pepIon.ModSequenceID
                    PrecursorMZ                                 = averagePSM.MeanPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = labeledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-labeledPeptide.Mass)
                    MeanPercolatorScore                         = averagePSM.MeanScore
                    QValue                                      = bestQValue
                    PEPValue                                    = bestPepValue
                    ProteinNames                                = prots
                    QuantMz_Light                               = mzLight
                    Quant_Light                                 = successfulQuant.Area
                    MeasuredApex_Light                          = successfulQuant.MeasuredApexIntensity
                    Seo_Light                                   = successfulQuant.StandardErrorOfPrediction
                    Params_Light                                = successfulQuant.EstimatedParams 
                    Difference_SearchRT_FittedRT_Light          = successfulQuant.SearchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Light            = successfulQuant.ClusterComparison.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Light   = successfulQuant.ClusterComparison.KLDiv_Corrected
                    QuantMz_Heavy                               = averagePSM.MeanPrecMz
                    Quant_Heavy                                 = quantP.Area
                    MeasuredApex_Heavy                          = quantP.MeasuredApexIntensity
                    Seo_Heavy                                   = quantP.StandardErrorOfPrediction
                    Params_Heavy                                = quantP.EstimatedParams 
                    Difference_SearchRT_FittedRT_Heavy          = searchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Heavy            = clusterComparisonTarget.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterComparisonTarget.KLDiv_Corrected
                    Correlation_Light_Heavy                     = successfulQuant.Correlation_Light_Heavy
                    QuantificationSource                        = QuantificationSource.PSM
                    IsotopicPatternMz_Light                     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Light     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Light    = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Light                               = successfulQuant.X_Xic 
                    IntensityTrace_Observed_Light               = successfulQuant.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Light              = successfulQuant.Y_Xic
                    IsotopicPatternMz_Heavy                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Heavy     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Heavy    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Heavy                               = averagePSM.X_Xic 
                    IntensityTrace_Observed_Heavy               = averagePSM.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Heavy              = averagePSM.Y_Xic
                    AlignmentScore                              = nan
                    AlignmentQValue                             = 0.
                    IonMobility                                 = nan
                    }
                    |> Option.Some
                | None ->
                    if diagCharts then 
                        saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                            peakToQuantify.XData peakToQuantify.YData quantP.YPredicted [||] [||] [||] [||] peaks clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = pepIon.Sequence
                    GlobalMod                                   = pepIon.GlobalMod
                    Charge                                      = pepIon.Charge
                    PepSequenceID                               = pepIon.PepSequenceID
                    ModSequenceID                               = pepIon.ModSequenceID
                    PrecursorMZ                                 = averagePSM.MeanPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = labeledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-labeledPeptide.Mass)
                    MeanPercolatorScore                         = averagePSM.MeanScore
                    QValue                                      = bestQValue
                    PEPValue                                    = bestPepValue
                    ProteinNames                                = prots
                    QuantMz_Light                               = mzLight
                    Quant_Light                                 = nan
                    MeasuredApex_Light                          = nan
                    Seo_Light                                   = nan
                    Params_Light                                = [||] 
                    Difference_SearchRT_FittedRT_Light          = nan
                    KLDiv_Observed_Theoretical_Light            = nan
                    KLDiv_CorrectedObserved_Theoretical_Light   = nan
                    QuantMz_Heavy                               = averagePSM.MeanPrecMz
                    Quant_Heavy                                 = quantP.Area
                    MeasuredApex_Heavy                          = quantP.MeasuredApexIntensity
                    Seo_Heavy                                   = quantP.StandardErrorOfPrediction
                    Params_Heavy                                = quantP.EstimatedParams 
                    Difference_SearchRT_FittedRT_Heavy          = searchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Heavy            = clusterComparisonTarget.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterComparisonTarget.KLDiv_Corrected
                    Correlation_Light_Heavy                     = nan
                    QuantificationSource                        = QuantificationSource.PSM
                    IsotopicPatternMz_Light                     = [||] 
                    IsotopicPatternIntensity_Observed_Light     = [||] 
                    IsotopicPatternIntensity_Corrected_Light    = [||] 
                    RtTrace_Light                               = [||] 
                    IntensityTrace_Observed_Light               = [||] 
                    IntensityTrace_Corrected_Light              = [||] 
                    IsotopicPatternMz_Heavy                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Heavy     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Heavy    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Heavy                               = averagePSM.X_Xic 
                    IntensityTrace_Observed_Heavy               = averagePSM.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Heavy              = averagePSM.Y_Xic
                    AlignmentScore                              = nan
                    AlignmentQValue                             = 0.
                    IonMobility                                 = nan
                    }
                    |> Option.Some
            with
            | ex ->             
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None

        let lableFreeQuantification (pepIon:PeptideIon)  (psms:PSMStatisticsResult []) = 
            try
            if pepIon.GlobalMod <> 0 then None 
            else
            let bestQValue,bestPepValue,prots = psms |> Array.minBy (fun x -> x.QValue) |> fun x -> x.QValue, x.PEPValue,x.ProteinNames
            let unlabledPeptide = peptideLookUp pepIon.Sequence 0
            let psmsWithMatchedSums = countMatchedMasses unlabledPeptide psms 
            let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
            let theoMz = Mass.toMZ unlabledPeptide.Mass (float pepIon.Charge)
            let averagePSM = average (fun _ -> nan) getXICWithIM scanTimeToMzCorrection theoMz psmsWithMatchedSums
            let avgMass = Mass.ofMZ (averagePSM.MeanPrecMz) (pepIon.Charge |> float)      
            let peaks = identifyPeaks averagePSM.X_Xic averagePSM.Y_Xic               
            if Array.isEmpty peaks then 
                if diagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "noPeaks" plotDirectory
                logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                None
            else
            let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks averagePSM.WeightedAvgScanTime
            let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify
            let searchRTMinusFittedRT = searchRTMinusFittedRtTarget averagePSM.WeightedAvgScanTime quantP 
            if quantP.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantP.EstimatedParams then 
                if diagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "fittingFailed" plotDirectory    
                logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                None
            else
            let clusterComparisonTarget = comparePredictedAndMeasuredIsotopicCluster averagePSM.X_Xic averagePSM.Y_Xic averagePSM.Y_Xic_uncorrected pepIon.Charge unlabledPeptide quantP.EstimatedParams.[1] averagePSM.MeanPrecMz
            if diagCharts then 
                saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                    peakToQuantify.XData peakToQuantify.YData quantP.YPredicted  [||] [||] [||] [||] peaks clusterComparisonTarget.PeakComparisons plotDirectory
            {
            StringSequence                              = pepIon.Sequence
            GlobalMod                                   = pepIon.GlobalMod
            Charge                                      = pepIon.Charge
            PepSequenceID                               = pepIon.PepSequenceID
            ModSequenceID                               = pepIon.ModSequenceID
            PrecursorMZ                                 = averagePSM.MeanPrecMz
            MeasuredMass                                = avgMass
            TheoMass                                    = unlabledPeptide.Mass
            AbsDeltaMass                                = abs(avgMass-unlabledPeptide.Mass)
            MeanPercolatorScore                         = averagePSM.MeanScore
            QValue                                      = bestQValue
            PEPValue                                    = bestPepValue
            ProteinNames                                = prots
            QuantMz_Light                               = averagePSM.MeanPrecMz
            Quant_Light                                 = quantP.Area
            MeasuredApex_Light                          = quantP.MeasuredApexIntensity
            Seo_Light                                   = quantP.StandardErrorOfPrediction
            Params_Light                                = quantP.EstimatedParams 
            Difference_SearchRT_FittedRT_Light          = searchRTMinusFittedRT
            KLDiv_Observed_Theoretical_Light            = clusterComparisonTarget.KLDiv_UnCorrected
            KLDiv_CorrectedObserved_Theoretical_Light   = clusterComparisonTarget.KLDiv_Corrected
            QuantMz_Heavy                               = nan
            Quant_Heavy                                 = nan
            MeasuredApex_Heavy                          = nan
            Seo_Heavy                                   = nan
            Params_Heavy                                = [||]
            Difference_SearchRT_FittedRT_Heavy          = nan
            KLDiv_Observed_Theoretical_Heavy            = nan
            KLDiv_CorrectedObserved_Theoretical_Heavy   = nan
            Correlation_Light_Heavy                     = nan
            QuantificationSource                        = QuantificationSource.PSM
            IsotopicPatternMz_Light                     = clusterComparisonTarget.PeakComparisons|> Array.map (fun x -> x.Mz)
            IsotopicPatternIntensity_Observed_Light     = clusterComparisonTarget.PeakComparisons|> Array.map (fun x -> x.MeasuredIntensity)
            IsotopicPatternIntensity_Corrected_Light    = clusterComparisonTarget.PeakComparisons|> Array.map (fun x -> x.MeasuredIntensityCorrected)
            RtTrace_Light                               = averagePSM.X_Xic 
            IntensityTrace_Observed_Light               = averagePSM.Y_Xic_uncorrected
            IntensityTrace_Corrected_Light              = averagePSM.Y_Xic
            IsotopicPatternMz_Heavy                     = [||]
            IsotopicPatternIntensity_Observed_Heavy     = [||]
            IsotopicPatternIntensity_Corrected_Heavy    = [||]
            RtTrace_Heavy                               = [||]
            IntensityTrace_Observed_Heavy               = [||]
            IntensityTrace_Corrected_Heavy              = [||]
            AlignmentScore                              = nan
            AlignmentQValue                             = 0.
            IonMobility                                 = nan
            }
            |> Option.Some
            with
            | ex ->
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None
        logger.Trace "init quantification functions:finished"
        
        logger.Trace "executing quantification"
        let quantResults = 
            qpsms 
            |> Array.groupBy (fun x -> 
                {
                    Sequence             = x.StringSequence     
                    GlobalMod            = x.GlobalMod    
                    Charge               = x.Charge       
                    ModSequenceID        = x.ModSequenceID
                    PepSequenceID        = x.PepSequenceID
                }
                )
            |> Array.map (fun (pepIon,psms) -> 
                    match processParams.XicExtraction.TopKPSMs with 
                    | Some x when psms.Length > x -> 
                        pepIon,
                        psms
                        |> Array.sortByDescending (fun x -> x.SequestScore)
                        |> Array.take x
                    | _ -> 
                        pepIon,
                        psms                  
                )
            |> Array.mapi (fun i (pepIon,psms) -> 
                if i % 100 = 0 then logger.Trace (sprintf "%i peptides quantified" i)
                
                match processParams.PerformLabeledQuantification with 
                |Domain.Labeling.N15Labeling | Domain.Labeling.N15LabelingOnly | Domain.Labeling.Labelshift-> labledQuantification pepIon psms
                |Domain.Labeling.Unlabeled -> lableFreeQuantification pepIon psms
                )
            |> Array.choose id
                    
        let filteredResults = 
            match processParams.PerformLabeledQuantification with 
            |Domain.Labeling.N15Labeling ->
                quantResults
                |> heavyQualityFilter -2. 2.
                |> lightQualityFilter -2. 2.
            |Domain.Labeling.N15LabelingOnly ->
                quantResults
                |> heavyQualityFilter -2. 2.
            |Domain.Labeling.Unlabeled ->
                quantResults
                |> lightQualityFilter -2. 2.
            |Domain.Labeling.Labelshift ->
                quantResults
        filteredResults 
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        if zipCharts then
            plotDirectory
            |> Zipping.zipDirectory "*.html" logger
            |> fun zipped ->
                match zipped with
                | Ok byteArr ->
                    byteArr
                    |> Zipping.saveZippedDirectory outputDir logger ((Path.GetFileNameWithoutExtension instrumentOutput) + "_plots")
                    |> fun saved ->
                        match saved with
                        | Ok save -> 
                            save
                            Directory.Delete (plotDirectory, true)
                        | Error ex -> logger.Trace (sprintf "Error saving zipped directory: %A" ex)
                | Error ex -> logger.Trace (sprintf "Error zipping directory: %A" ex)
        logger.Trace (sprintf "Isotope patterns: %i read from the peptide data base, %i computed" isotopePatterns.Read isotopePatterns.Computed)
        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()
        (isotopePatterns :> System.IDisposable).Dispose()
        memoryDB.Dispose()
        logger.Trace "executing quantification:finished"
        