namespace ProteomIQon

open System.IO
open Argu
open System
open ProteomIQon.Core
open Core.MzIO
open Dto
open FSharp.Stats
open BioFSharp.Mz.Quantification
open BioFSharp.Mz
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute
open Plotly.NET
open BioFSharp
open MzIO.Processing
open BioFSharp.Mz.SearchDB
open System.Data

module AlignmentBasedQuantification =
    
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

    ///
    type AlignmentQuantMetrics = 
       {
           [<FieldAttribute(0)>]
           Sequence                             : string
           [<FieldAttribute(1)>]
           GlobalMod                            : int
           [<FieldAttribute(2)>]
           Charge                               : int
           [<FieldAttribute(3)>]
           PepSequenceID                        : int
           [<FieldAttribute(4)>]
           ModSequenceID                        : int
           [<FieldAttribute(5)>]
           X_FileName                           : string 
           [<FieldAttribute(6)>]
           X_Intensities                        : float 
           [<FieldAttribute(7)>]
           X_Stabw                              : float
           [<FieldAttribute(8)>]
           X_Test                               : float 
           [<FieldAttribute(9)>] [<TraceConverter>]
           X_IsotopicPatternMz                  : float []
           [<FieldAttribute(10)>] [<TraceConverter>]
           X_IsotopicPatternIntensity_Observed  : float []
           [<FieldAttribute(11)>] [<TraceConverter>]
           X_RtTrace                            : float []
           [<FieldAttribute(12)>] [<TraceConverter>]
           X_IntensityTrace                     : float []   
           [<FieldAttribute(13)>]
           Y_Test                               : float 
           [<FieldAttribute(14)>]
           YHat_Test                            : float 
           [<FieldAttribute(15)>]
           YHat_Refined_Test                    : float
           [<FieldAttribute(16)>]
           Y_Intensity                          : float
           [<FieldAttribute(17)>] [<TraceConverter>]
           Y_IsotopicPatternMz                  : float []
           [<FieldAttribute(18)>] [<TraceConverter>]
           Y_IsotopicPatternIntensity_Observed  : float [] 
           [<FieldAttribute(19)>] [<TraceConverter>]
           Y_RtTrace                            : float []
           [<FieldAttribute(20)>] [<TraceConverter>]
           Y_IntensityTrace                     : float []
           [<FieldAttribute(21)>]
           DtwDistanceBefore                    : float 
           [<FieldAttribute(22)>]
           Y_ReQuant                            : float 
       }

    type PeakComparison = {
        Mz: float
        MeasuredIntensity: float
        MeasuredIntensityCorrected:float
        PredictedRelFrequency:float
        }

    type ClusterComparison = {
        PeakComparisons     : PeakComparison []
        KLDiv_UnCorrected   : float
        KLDiv_Corrected     : float
        }


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
    let getBaseLineCorrectionOffsetAt tarRT x_Xic y_Xic y_Xic_uncorrected =
        let (rt,y,yUncorr) = 
            Array.zip3 x_Xic y_Xic y_Xic_uncorrected 
            |> Array.minBy (fun (rt,y,yUncorr) -> abs (rt - tarRT))
        yUncorr - y
       
    ///
    let getClosestMs1 (ms1s: (float*MzIO.Model.MassSpectrum) []) scanTime = 
        ms1s
        |> Array.minBy (fun ms -> abs (fst ms - scanTime))
        |> snd

    ///
    let getSpec (reader:MzIO.IO.IMzIODataReader) (ms1: MzIO.Model.MassSpectrum)  =
        Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(ms1.ID).Peaks)
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
                |> log2
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
                |> log2
            (qualR > lowerBorder && qualR < upperBorder) || (nan.Equals(qualR) )
            )
    
    /// Can be used to filter out all QuantResults whos absolute difference between fitted scan time and initial scantime
    /// exceed a given (positive) standard deviation factor (stDevFactor).
    /// Example: if the search scan time was 30. and fitted Scantime 30.5 the difference between
    /// fitted scan time and initial scantime is 0.5. If the standard deviation of the peak
    /// was 0.25, factor is 2. If stDevFactor is 1.5 this quantification is discarded, if it is 3. the quantification is kept.
    let lightScanTimeDifferenceFilter stDevFactor (quantResults:QuantificationResult[]) =
        if stDevFactor <= 0. then failwith "stabwfactor has to be positive." else
        quantResults
        |> Array.filter (fun qP -> 
            match QuantificationResult.tryLightGetStabw qP with 
            | Some stdev -> 
                let differenceStDevRatio = (abs qP.Difference_SearchRT_FittedRT_Light) / stdev
                differenceStDevRatio < stDevFactor
            | None -> false
            )    

    /// Can be used to filter out all QuantResults whos absolute difference between fitted scan time and initial scantime
    /// exceed a given (positive) standard deviation factor (stDevFactor).
    /// Example: if the search scan time was 30. and fitted Scantime 30.5 the difference between
    /// fitted scan time and initial scantime is 0.5. If the standard deviation of the peak
    /// was 0.25, factor is 2. If stDevFactor is 1.5 this quantification is discarded, if it is 3. the quantification is kept.
    let heavyScanTimeDifferenceFilter stDevFactor (quantResults:QuantificationResult[]) =
        if stDevFactor <= 0. then failwith "stabwfactor has to be positive." else
        quantResults
        |> Array.filter (fun qP -> 
            match QuantificationResult.tryHeavyGetStabw qP with 
            | Some stdev -> 
                let differenceStDevRatio = (abs qP.Difference_SearchRT_FittedRT_Heavy) / stdev
                differenceStDevRatio < stDevFactor
            | None -> false
            )
            
    /// Calculates the Kullback-Leibler divergence Dkl(p||q) from q (theory, model, description, or approximation of p) 
    /// to p (the "true" distribution of data, observations, or a precisely calculated theoretical distribution).
    let klDiv (p:float []) (q:float []) = 
        Array.fold2 (fun acc p q -> (System.Math.Log(p/q)*p) + acc ) 0. p q
     
    ///
    let substractBaseLine (logger: NLog.Logger) (baseLineParams:Domain.BaseLineCorrection) (yData:float []) =
        if yData.Length > 500 then
            yData
        else
            let baseLine = FSharp.Stats.Signal.Baseline.baselineAls' baseLineParams.MaxIterations baseLineParams.Lambda baseLineParams.P yData |> Array.ofSeq
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

    ///
    let initGetIsotopicEnvelope reader idx scanTimeWindow mzWindow_Da ch meanScanTime meanPrecMz =
        let rtQuery   = Query.createRangeQuery meanScanTime scanTimeWindow
        let mzQueries = 
            [|
                for i = 0 to 2 do 
                    let mz = meanPrecMz + (float i) * (Mass.Table.PMassInU / (ch|> float)) 
                    Query.createRangeQuery mz mzWindow_Da 
            |]
            |> Array.filter (fun x ->  abs (meanPrecMz-x.LockValue) < 1.)                 
        let retData',itzData' =
            let query = Query.getXICs reader idx rtQuery mzQueries 
            [|
                for i = 0 to query.[*,0].Length-1 do 
                let tmp = query.[i,*]
                yield (tmp.[0].Rt,tmp |> Array.sumBy (fun p -> p.Intensity))

            |]
            |> Array.unzip 
        retData',itzData'

    type InferredPeak = {
        Model                       :HULQ.PeakModel option
        Area                        :float
        StandardErrorOfPrediction   :float
        MeasuredApexIntensity       :float
        EstimatedParams             :float[]
        X_Xic                       :float[]
        Y_Xic                       :float[]
        Y_Xic_uncorrected           :float[]
        xPeak                       :float[]
        yFitted                     :float[]
        InferredScanTimeRefined     : float
        AlignedSource               :(float*float) list
        }

    type XIC = {
        FinalPrecMz         : float
        FinalTargetScanTime : float
        X_Xic               : float []
        Y_Xic               : float []
        Y_Xic_uncorrected   : float []
        AlignedSource       : (float*float) []
        }

    ///
    let getRefinedXic refineByDTW getXic xSource ySource sourceScanTime scanTimeToMzCorrection theoMz scanTimePrediciton =
            let correctedMz = scanTimeToMzCorrection scanTimePrediciton + theoMz
            let (retData,itzDataCorrected,ItzDataUncorrected) = getXic scanTimePrediciton correctedMz
            let maxItz = itzDataCorrected |> Array.max
            let inferredScanTime, alignedSource = 
                if not refineByDTW then 
                    scanTimePrediciton,[||]
                else 
                    let target = Array.zip retData (DTW'.zNorm itzDataCorrected)
                    let source = Array.zip xSource (DTW'.zNorm ySource)
                    let alignedSource = DTW'.align target source 
                    let inferredScanTimeRefined = DTW'.align' target source sourceScanTime |> snd
                    let normFactor = alignedSource |> List.maxBy snd |> snd
                    inferredScanTimeRefined, alignedSource |> Array.ofList |> Array.map (fun (x,y) -> x,(y / normFactor) * maxItz)
            {
                FinalPrecMz         = correctedMz
                FinalTargetScanTime = inferredScanTime
                X_Xic               = retData
                Y_Xic               = itzDataCorrected
                Y_Xic_uncorrected   = ItzDataUncorrected
                AlignedSource       = alignedSource
            }

    type InferredXic = {
        X_Xic                       :float[]
        Y_Xic                       :float[]
        Y_Xic_uncorrected           :float[]
        }

    let getInferredXic getXic targetScanTime targetMz =
        let (retData,itzData,uncorrectedItzData)   =
                getXic targetScanTime targetMz 
        {
        X_Xic               = retData
        Y_Xic               = itzData
        Y_Xic_uncorrected   = uncorrectedItzData
        }


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

    /// Calculates the difference between the scan time used to retreave the peak and the fitted peak midpoint.
    let searchRTMinusFittedRtTarget searchRT (fit:HULQ.QuantifiedPeak) = 
        try
            searchRT - fit.EstimatedParams.[1] 
        with
        | _ -> nan

    ///
    let saveMetricChart sequence globalMod ch (xXic:float[]) (yXic:float[]) avgScanTime (xScantime:float) 
            (xToQuantify:float[]) (ypToQuantify:float[]) (fitY:float[]) 
            inferredScanTimeRefined alignedSource
            (xXicInferred:float[]) (yXicinferred:float[]) (xInferred:float[]) (inferredFit:float[]) 
                (pattern:PeakComparison []) plotDirectory =
        let xic = 
            [
            Chart.Point(xXic, yXic)                     |> Chart.withTraceName "Target XIC"
            Chart.Point(alignedSource)                  |> Chart.withTraceName "Aligned Source XIC"
            Chart.Point([avgScanTime],[1.])             |> Chart.withTraceName "Inferred Scan Time"
            Chart.Point([inferredScanTimeRefined],[1.]) |> Chart.withTraceName "Inferred Scan Time Refined"
            Chart.Point([xScantime],[1.])               |> Chart.withTraceName "ScanTime Source"
            Chart.Point((xToQuantify), (ypToQuantify))  |> Chart.withTraceName "Identified Target Peak"
            Chart.Line(xToQuantify,fitY)                |> Chart.withTraceName "Fit of target Peak"
            Chart.Point(xXicInferred, yXicinferred)     |> Chart.withTraceName "Inferred XIC"
            Chart.Line(xInferred,inferredFit)           |> Chart.withTraceName "Fit of inferred Peak"

            ]
            |> Chart.Combine
        let pattern = 
            [
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensity))          |> Chart.withTraceName "Measured"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensityCorrected)) |> Chart.withTraceName "Measured Corrected"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.PredictedRelFrequency))      |> Chart.withTraceName "Predicted Relative Frequency"
            ]
            |> Chart.Combine
        [xic;pattern]
        |> Chart.Stack(2, 0.1)
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i" sequence globalMod)
        |> Chart.withSize(2500.,800.)
        |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + globalMod.ToString() + "Ch" + ch.ToString())|])
        
    ///
    let saveChart sequence globalMod ch (xXic:float[]) (yXic:float[]) avgScanTime 
            (xToQuantify:float[]) (ypToQuantify:float[]) (fitY:float[]) 
            inferredScanTimeRefined alignedSource
            (xXicInferred:float[]) (yXicinferred:float[]) (xInferred:float[]) (inferredFit:float[]) 
                (pattern:PeakComparison []) plotDirectory =
        let xic = 
            [
            Chart.Point(xXic, yXic)                     |> Chart.withTraceName "Target XIC"
            Chart.Point(alignedSource)                  |> Chart.withTraceName "Aligned Source XIC"
            Chart.Point([avgScanTime],[1.])             |> Chart.withTraceName "Inferred Scan Time"
            Chart.Point([inferredScanTimeRefined],[1.]) |> Chart.withTraceName "Inferred Scan Time Refined"
            Chart.Point((xToQuantify), (ypToQuantify))  |> Chart.withTraceName "Identified Target Peak"
            Chart.Line(xToQuantify,fitY)                |> Chart.withTraceName "Fit of target Peak"
            Chart.Point(xXicInferred, yXicinferred)     |> Chart.withTraceName "Inferred XIC"
            Chart.Line(xInferred,inferredFit)           |> Chart.withTraceName "Fit of inferred Peak"

            ]
            |> Chart.Combine
        let pattern = 
            [
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensity))          |> Chart.withTraceName "Measured"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensityCorrected)) |> Chart.withTraceName "Measured Corrected"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.PredictedRelFrequency))      |> Chart.withTraceName "Predicted Relative Frequency"
            ]
            |> Chart.Combine
        [xic;pattern]
        |> Chart.Stack(2, 0.1)
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i" sequence globalMod)
        |> Chart.withSize(2500.,800.)
        |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + globalMod.ToString() + "Ch" + ch.ToString())|])
            
    let saveErrorChart (xXic:float[]) (yXic:float[]) pepSeq gMod ch desc plotDirectory =      
        Chart.Point(xXic, yXic)
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i_%s" pepSeq gMod desc)
        |> Chart.withSize(1500.,800.)
        |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((pepSeq |> String.filter (fun x -> x <> '*')) + "_GMod_" + gMod.ToString() + "Ch" + ch.ToString() + "_notQuantified")|])
               
    // Method is based on: https://doi.org/10.1021/ac0600196
    /// Estimates the autocorrelation at lag 1 of a blank signal (containing only noise). Subsequently, the signal of interest is smoothed
    /// several times by a savitzky golay filter using constant polynomial order and variing windowWidth. For each iteration, the deviation
    /// of the smoothed to the original signal is computed and the autocorrelation at lag 1 of this residual noise is computed. The function returns the optimized
    /// window width yielding a autocorrelation at lag 1 closest to the value computed for the blank signal.
    let optimizeWindowWidth polOrder (windowWidthToTest:int[]) noiseAutoCorr (signalOfInterest:float[]) =
        let signalOfInterest' = signalOfInterest |> vector
        //let noiseAutoCorr = Correlation.Vector.autoCorrelation 1 (blankSignal |> vector)
        let filterF w yData = FSharp.Stats.Signal.Filtering.savitzkyGolay w polOrder 0 0 yData |> vector
        let windowWidthToTest' = windowWidthToTest |> Array.filter (fun x -> x%2 <> 0)
        let optimizedWindowWidth =
            windowWidthToTest'
            |> Array.map (fun w ->
                let smoothedY = filterF w signalOfInterest
                let noise = smoothedY - (signalOfInterest')
                w, Correlation.Vector.autoCorrelation 1 noise
                )
            |> Array.minBy (fun (w,ac) -> (ac - noiseAutoCorr) |> abs )
            |> fst
        optimizedWindowWidth
    
    ///
    let initGetWindowWidth (windowEst:Domain.WindowSize) polynomOrder (windowWidthToTest:int[]) =
        match windowEst with
        | Domain.WindowSize.Fixed w  -> fun yData -> w
        | Domain.WindowSize.EstimateUsingAutoCorrelation noiseAutoCorr -> fun yData -> optimizeWindowWidth polynomOrder windowWidthToTest noiseAutoCorr yData
    
    ///
    let initIdentifyPeaks (peakDetectionParams:Domain.XicProcessing) =
        match peakDetectionParams with 
        | Domain.XicProcessing.SecondDerivative parameters ->
            let getWindowWidth = initGetWindowWidth parameters.WindowSize parameters.PolynomOrder [|5 .. 2 .. 60|] 
            (fun xData yData -> 
                let  windowSize = getWindowWidth yData
                FSharp.Stats.Signal.PeakDetection.SecondDerivative.getPeaks parameters.MinSNR parameters.PolynomOrder windowSize xData yData
                )
        | Domain.XicProcessing.Wavelet parameters ->
            (fun xData yData -> 
                FSharpStats'.Wavelet.identify parameters xData yData
                )
        
    ///
    let searchRTMinusFittedRt searchRT fit = 
        try
            searchRT - fit.EstimatedParams.[1] 
        with
        | _ -> nan

    /// Aims to provide the best scan time estimate when quanitfying a inferred peak.
    let chooseScanTime maxDiff searchRTMinusFittedRT initialScanTime (quantP:HULQ.QuantifiedPeak) = 
        if Array.isEmpty quantP.EstimatedParams then
            initialScanTime 
        elif abs searchRTMinusFittedRT > maxDiff then
            initialScanTime 
        else
            quantP.EstimatedParams.[1]

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

    ///Predicts an isotopic distribution of the given formula at the given charge, normalized by the sum of probabilities, using the MIDAs algorithm
    let generateIsotopicDistributionOfFormulaBySum (charge:int) (seq:AminoAcids.AminoAcid list) =
        seq
        |> BioList.toFormula
        |> Formula.add Formula.Table.H2O
        |> IsotopicDistribution.MIDA.ofFormula IsotopicDistribution.MIDA.normalizeByProbSum 0.01 0.001 charge 
        |> Array.ofList

    ///
    let initComparePredictedAndMeasuredIsotopicCluster inReader ms1s ms1AccuracyEstimate (x_Xic:float[]) (y_Xic:float[]) y_Xic_uncorrected ch peptideSequence tarRt tarMz =    
        /// IsotopicCluster
        let targetIsotopicPattern_predicted = 
            generateIsotopicDistributionOfFormulaBySum ch peptideSequence
        let baseLineCorrectionF = getBaseLineCorrectionOffsetAt tarRt x_Xic y_Xic y_Xic_uncorrected
        let closestMS1 = getClosestMs1 ms1s tarRt
        let peaks' = 
            getSpec inReader closestMS1
            |> Array.filter (fun x -> x.Mz < tarMz + 1. && x.Mz > tarMz - 0.6)
        let recordedVsPredictedPattern = 
            targetIsotopicPattern_predicted
            |> Array.choose (fun (mz,relFreq) ->
                if peaks' |> Array.isEmpty then None
                else
                    let closestRealPeak = peaks' |> Array.minBy (fun peak -> abs(peak.Mz - mz)) 
                    if (abs(closestRealPeak.Mz - mz) < 4. * ms1AccuracyEstimate) then  
                        Some (closestRealPeak,relFreq)
                    else None
                )
            |> Array.groupBy fst
            |> Array.map (fun ((peak),list) -> 
                let (mz,measuredIntensity,predictedRelFrequency) = peak.Mz,peak.Intensity,list |> Array.sumBy snd 
                {Mz=mz;MeasuredIntensity=measuredIntensity;MeasuredIntensityCorrected=measuredIntensity - baseLineCorrectionF;PredictedRelFrequency= predictedRelFrequency}
                )
            |> Array.filter (fun (isoP:PeakComparison) -> isoP.MeasuredIntensityCorrected > 0.)
        let recordedVsPredictedPatternNorm = 
            let sumMeasured = recordedVsPredictedPattern |> Array.sumBy (fun x -> x.MeasuredIntensity)
            let sumMeasuredCorr = recordedVsPredictedPattern |> Array.sumBy (fun x -> x.MeasuredIntensityCorrected)
            let sumMeasuredPred = recordedVsPredictedPattern |> Array.sumBy (fun x -> x.PredictedRelFrequency)
            recordedVsPredictedPattern
            |>Array.map (fun isoP -> 
                {isoP with 
                    MeasuredIntensity = isoP.MeasuredIntensity / sumMeasured; 
                    MeasuredIntensityCorrected = isoP.MeasuredIntensityCorrected / sumMeasuredCorr; 
                    PredictedRelFrequency = isoP.PredictedRelFrequency / sumMeasuredPred }
                )
        let klUnCorr = klDiv (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.MeasuredIntensity)) (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.PredictedRelFrequency))
        let klCorr   = klDiv (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.MeasuredIntensityCorrected)) (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.PredictedRelFrequency))
        {
            PeakComparisons     = recordedVsPredictedPatternNorm
            KLDiv_UnCorrected   = klUnCorr
            KLDiv_Corrected     = klCorr
        }
    
    
    ///
    let quantifyPeptides diagCharts (processParams:Domain.AlignmentBasedQuantificationParams) (outputDir:string) (d:string) (instrumentOutput:string) (psmbasedQuant:string) (alignmentMetrics:string) (alignmentResults:string)  =

        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension alignmentResults)
        logger.Trace (sprintf "Input file: %s" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        logger.Trace (sprintf "Now performing Quantification using: %s and %s, Results will be written to: %s" instrumentOutput alignmentResults outputDir)
        // initialize Reader and Transaction
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".quant"
            Path.Combine [|outputDir;fileName|]
        logger.Trace (sprintf "outFilePath:%s" outFilePath)
        let outFilePathMetrics =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".alignquantMetrics"
            Path.Combine [|outputDir;fileName|]
        logger.Trace (sprintf "outFilePath:%s" outFilePath)
        //
        let plotDirectory =
            let fileName = sprintf "%s_plots" (Path.GetFileNameWithoutExtension instrumentOutput)
            let path = Path.Combine [|outputDir;fileName|]
            if System.IO.Directory.Exists path then
                path
            else
                System.IO.Directory.CreateDirectory path |> ignore
                path
        let metricPlotDirectory =
            let fileName = sprintf "%s_metricPlots" (Path.GetFileNameWithoutExtension instrumentOutput)
            let path = Path.Combine [|outputDir;fileName|]
            if System.IO.Directory.Exists path then
                path
            else
                System.IO.Directory.CreateDirectory path |> ignore
                path
        logger.Trace (sprintf "plotDirectory:%s" plotDirectory)
        
        logger.Trace "Copy peptide DB into Memory"
        let cn =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        let memoryDB = SearchDB.copyDBIntoMemory cn
        cn.Dispose()
        logger.Trace "Copy peptide DB into Memory: finished"
        
        logger.Trace "Get peptide lookUp function"
        let dBParams = SearchDB'.getSDBParams memoryDB
        let peptideLookUp = SearchDB'.getThreadSafePeptideLookUpFromFileBySequenceAndGMod memoryDB dBParams
        let calcIonSeries aal = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal
        logger.Trace "Get peptide lookUp function: finished"
        // initialize Reader and Transaction
        logger.Trace "Init connection to mass spectrum data."
        let inReader = Core.MzIO.Reader.getReader instrumentOutput :?> MzIO.MzSQL.MzSQL
        inReader.Connection.Open()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader
        let inTr = inReader.BeginTransaction()

        logger.Trace "Create RetentionTime index"
        let retTimeIdxed = Query.getMS1RTIdx inReader inRunID
        logger.Trace "Create RetentionTime index:finished"
        
        logger.Trace "Read and sort ms1s"
        /// Returns a Sequence containing all MassSpectra of a single MS-run
        let massSpectra = inReader.ReadMassSpectra(inRunID)
        /// Returns a Array that contains all MS1s sorted by their scanTime
        let ms1SortedByScanTime =
            massSpectra
            |> Seq.filter (fun ms -> MassSpectrum.getMsLevel ms = 1)
            |> Seq.map (fun ms -> MassSpectrum.getScanTime ms, ms)
            |> Seq.sortBy fst
            |> Array.ofSeq
        logger.Trace "Read and sort ms1s:finished"

        logger.Trace "Read alignmentResults"
        ///
        let psmbasedQuant =
            Csv.CsvReader<QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(psmbasedQuant,'\t',false,1)
            |> Array.ofSeq
            |> Array.filter (fun x -> QuantificationResult.tryTargetGetScanTime x |> Option.isSome)
        logger.Trace "Read scored alignmentResult:finished"        

        logger.Trace "Read alignmentResults"
        ///
        let alignmentMetrics =
            Csv.CsvReader<AlignmentMetricsDTO>(SchemaMode=Csv.Fill).ReadFile(alignmentMetrics,'\t',false,1)
            |> Array.ofSeq
        logger.Trace "Read scored alignmentResult:finished"        

        logger.Trace "Read alignmentResults"
        ///
        let alignmentResults =
            Csv.CsvReader<AlignmentResult>(SchemaMode=Csv.Fill).ReadFile(alignmentResults,'\t',false,1)
            |> Array.ofSeq
        logger.Trace "Read scored alignmentResult:finished"

        logger.Trace "Estimate precursor mz standard deviation and mz correction."
        /// auf quantresult umstellen.
        let ms1AccuracyEstimate,scanTimeToMzCorrection = 
            let scanTimeVsDelta = 
                psmbasedQuant
                //|> Array.filter (Dto.QuantificationResult.tryTargetGetScanTime >> Option.isSome)
                |> Array.choose (fun x -> 
                    match Dto.QuantificationResult.tryTargetGetScanTime x with 
                    | None -> None
                    | Some scanTime -> 
                        let precMz = x.PrecursorMZ
                        let theMz  = Mass.toMZ x.TheoMass (float x.Charge)
                        let diff = precMz - theMz 
                        (scanTime, diff)
                        |> Option.Some
                    )
                |> Array.filter (fun (x,y) -> (nan.Equals(x) |> not) && (nan.Equals(y) |> not))
                |> Array.sortBy fst
            //scanTimeVsDelta
            //|> Array.map (fun x -> sprintf "%f;%f" (fst x) (snd x))
            //|> FSharpAux.IO.SeqIO.Seq.writeOrAppend (@"C:\Users\Public\David\DataAnalysis\1_Pipeline\AlignQuantN14N15\booo.csv")
            let borders =
                scanTimeVsDelta
                |> Seq.map snd
                |> Array.ofSeq
                |> FSharp.Stats.Testing.Outliers.tukey 3.
            let filteredValues =
                if nan.Equals(borders.Upper) || nan.Equals(borders.Lower) then 
                    scanTimeVsDelta
                else
                    scanTimeVsDelta
                    |> Array.filter (fun (s,d) -> d <= borders.Upper && d >= borders.Lower)
            [
            Chart.Point(scanTimeVsDelta)
            |> Chart.withTraceName "Raw"
            Chart.Point(filteredValues)
            |> Chart.withTraceName "Filtered"
            ]
            |> Chart.Combine
            |> Chart.withTitle (sprintf "upper border:%f, lower border:%f" borders.Upper borders.Lower)
            |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; "mzCorrectionData"|])
            let scanTimeToMzCorrection =
                let runTime = scanTimeVsDelta |> Array.maxBy fst |> fst
                if runTime > 20. && psmbasedQuant.Length > 500 then
                    let binWidth = 
                        System.Math.Min(runTime / 2., 20.)
                    let stabw = filteredValues |> Seq.stDevBy snd
                    let r,f = initSpline binWidth filteredValues
                    f
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
                |> Chart.withTraceName "Raw"
                Chart.Line(scanTimeVsDelta |> Array.sortBy fst |> Array.map (fun (st,d) -> st, scanTimeToMzCorrection st))
                ]
                |> Chart.Combine
                |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; "mzErrorAndCorrection"|])
            //stDev, scanTimeToMzCorrection 
            0.05, scanTimeToMzCorrection 
            
        logger.Trace (sprintf "Estimate precursor mz standard deviation and mz correction.:finished, standard deviation: %f" ms1AccuracyEstimate) 

        
        logger.Trace "init lookup functions"   
        ///
        let comparePredictedAndMeasuredIsotopicCluster = initComparePredictedAndMeasuredIsotopicCluster inReader ms1SortedByScanTime ms1AccuracyEstimate     
        
        let readSpecPeaksWithMem = FSharpAux.Memoization.memoize inReader.ReadSpectrumPeaks
        let getPeaks = Query.initRTProfile readSpecPeaksWithMem
        ///
        let getXIC = initGetProcessedXIC logger processParams.BaseLineCorrection getPeaks retTimeIdxed processParams.XicExtraction.ScanTimeWindow ms1AccuracyEstimate     
        
        ///
        let identifyPeaks = initIdentifyPeaks processParams.XicExtraction.XicProcessing
        
        logger.Trace "init lookup functions:finished"
        
        logger.Trace "init quantification functions"  
        let quantifyTestDataSet (metrics:AlignmentMetricsDTO []) = 
            metrics
            |> Array.choose (fun testPep -> 
                try
                let unlabledPeptide = peptideLookUp testPep.Sequence 0
                let labeledPeptide  = peptideLookUp testPep.Sequence 1
                let targetPeptide = if testPep.GlobalMod = 0 then unlabledPeptide else labeledPeptide            
                let targetMz = Mass.toMZ (targetPeptide.Mass) (testPep.Charge|> float)
                let refinedXIC = 
                    getRefinedXic processParams.PerformLocalWarp getXIC 
                        testPep.X_RtTrace testPep.X_IntensityTrace testPep.X_Test 
                            scanTimeToMzCorrection targetMz testPep.YHat_Test
                let avgMass = Mass.ofMZ (targetMz) (testPep.Charge |> float)
                let peaks = identifyPeaks refinedXIC.X_Xic refinedXIC.Y_Xic
                if Array.isEmpty peaks then 
                    if diagCharts then saveErrorChart refinedXIC.X_Xic refinedXIC.Y_Xic testPep.Sequence testPep.GlobalMod testPep.Charge "noPeaks" plotDirectory
                    logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (testPep.Sequence |> String.filter (fun x -> x <> '*')) (testPep.GlobalMod.ToString()) (testPep.Charge.ToString()))
                    None
                else
                let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks refinedXIC.FinalTargetScanTime
                let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify                        
                let searchRTMinusFittedRT = searchRTMinusFittedRtTarget refinedXIC.FinalTargetScanTime quantP
                if quantP.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantP.EstimatedParams then 
                    if diagCharts then saveErrorChart refinedXIC.X_Xic refinedXIC.Y_Xic testPep.Sequence testPep.GlobalMod testPep.Charge "fittingFailed" plotDirectory    
                    logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (testPep.Sequence |> String.filter (fun x -> x <> '*')) (testPep.GlobalMod.ToString()) (testPep.Charge.ToString()))
                    None
                else
                let inferredScanTime = chooseScanTime processParams.XicExtraction.ScanTimeWindow searchRTMinusFittedRT refinedXIC.FinalTargetScanTime quantP 
                let clusterComparisonTarget = comparePredictedAndMeasuredIsotopicCluster refinedXIC.X_Xic refinedXIC.Y_Xic refinedXIC.Y_Xic_uncorrected testPep.Charge targetPeptide.BioSequence quantP.EstimatedParams.[1] refinedXIC.FinalPrecMz 
                if diagCharts then 
                    saveMetricChart testPep.Sequence testPep.GlobalMod testPep.Charge refinedXIC.X_Xic refinedXIC.Y_Xic testPep.Y_Test testPep.X_Test
                                peakToQuantify.XData quantP.YPredicted quantP.YPredicted 
                                refinedXIC.FinalTargetScanTime refinedXIC.AlignedSource
                                [||] [||] [||] [||]
                                clusterComparisonTarget.PeakComparisons metricPlotDirectory
                {
                    Sequence                             = testPep.Sequence                            
                    GlobalMod                            = testPep.GlobalMod                           
                    Charge                               = testPep.Charge                              
                    PepSequenceID                        = testPep.PepSequenceID                       
                    ModSequenceID                        = testPep.ModSequenceID                       
                    X_FileName                           = testPep.X_FileName                          
                    X_Intensities                        = testPep.X_Intensities                       
                    X_Stabw                              = testPep.X_Stabw                             
                    X_Test                               = testPep.X_Test                              
                    X_IsotopicPatternMz                  = testPep.X_IsotopicPatternMz                 
                    X_IsotopicPatternIntensity_Observed  = testPep.X_IsotopicPatternIntensity_Observed 
                    X_RtTrace                            = testPep.X_RtTrace                           
                    X_IntensityTrace                     = testPep.X_IntensityTrace                      
                    Y_Test                               = testPep.Y_Test                              
                    YHat_Test                            = testPep.YHat_Test                           
                    YHat_Refined_Test                    = testPep.YHat_Refined_Test                   
                    Y_Intensity                          = testPep.Y_Intensity                         
                    Y_IsotopicPatternMz                  = testPep.Y_IsotopicPatternMz                 
                    Y_IsotopicPatternIntensity_Observed  = testPep.Y_IsotopicPatternIntensity_Observed 
                    Y_RtTrace                            = testPep.Y_RtTrace                           
                    Y_IntensityTrace                     = testPep.Y_IntensityTrace                    
                    DtwDistanceBefore                    = testPep.DtwDistanceBefore                   
                    Y_ReQuant                            = quantP.Area
                }
                |> Some
                with
                | _ -> None
                )
                
        ///
        let labledQuantification (alignmentResult:AlignmentResult) = 
            try          
            let unlabledPeptide = peptideLookUp alignmentResult.StringSequence 0
            let labeledPeptide  = peptideLookUp alignmentResult.StringSequence 1
            let targetPeptide = if alignmentResult.GlobalMod = 0 then unlabledPeptide else labeledPeptide            
            let targetMz = Mass.toMZ (targetPeptide.Mass) (alignmentResult.Charge|> float)
            let refinedXIC = 
                getRefinedXic processParams.PerformLocalWarp getXIC 
                    alignmentResult.RtTrace_SourceFile alignmentResult.IntensityTrace_SourceFile alignmentResult.ScanTime_SourceFile 
                        scanTimeToMzCorrection targetMz alignmentResult.PredictedScanTime
            let avgMass = Mass.ofMZ (targetMz) (alignmentResult.Charge |> float)
            let peaks = identifyPeaks refinedXIC.X_Xic refinedXIC.Y_Xic
            if Array.isEmpty peaks then 
                if diagCharts then saveErrorChart refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge "noPeaks" plotDirectory
                logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (alignmentResult.StringSequence |> String.filter (fun x -> x <> '*')) (alignmentResult.GlobalMod.ToString()) (alignmentResult.Charge.ToString()))
                None
            else
            let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks refinedXIC.FinalTargetScanTime
            let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify                        
            let searchRTMinusFittedRT = searchRTMinusFittedRtTarget refinedXIC.FinalTargetScanTime quantP
            if quantP.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantP.EstimatedParams then 
                if diagCharts then saveErrorChart refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge "fittingFailed" plotDirectory    
                logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (alignmentResult.StringSequence |> String.filter (fun x -> x <> '*')) (alignmentResult.GlobalMod.ToString()) (alignmentResult.Charge.ToString()))
                None
            else
            let inferredScanTime = chooseScanTime processParams.XicExtraction.ScanTimeWindow searchRTMinusFittedRT refinedXIC.FinalTargetScanTime quantP 
            let clusterComparisonTarget = comparePredictedAndMeasuredIsotopicCluster refinedXIC.X_Xic refinedXIC.Y_Xic refinedXIC.Y_Xic_uncorrected alignmentResult.Charge targetPeptide.BioSequence quantP.EstimatedParams.[1] refinedXIC.FinalPrecMz            
            if alignmentResult.GlobalMod = 0 then
                let mzHeavy = 
                    let mz = Mass.toMZ (labeledPeptide.Mass) (alignmentResult.Charge|> float)
                    let correctedMz = scanTimeToMzCorrection inferredScanTime + mz
                    correctedMz 
                let inferredQuant = 
                    let inferredXicHeavy = getInferredXic getXIC refinedXIC.FinalTargetScanTime mzHeavy    
                    let inferredPeaksHeavy = identifyPeaks inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic
                    if Array.isEmpty inferredPeaksHeavy then 
                        if diagCharts then saveErrorChart refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge "noInferredPeaks" plotDirectory
                        logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (alignmentResult.StringSequence |> String.filter (fun x -> x <> '*')) (alignmentResult.GlobalMod.ToString()) (alignmentResult.Charge.ToString()))
                        None
                    else
                    let peakToQuantifyHeavy = BioFSharp.Mz.Quantification.HULQ.getPeakBy inferredPeaksHeavy inferredScanTime
                    let quantPHeavy         = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantifyHeavy          
                    if quantPHeavy.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantPHeavy.EstimatedParams then 
                        if diagCharts then saveErrorChart refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge "fittingInferredFailed" plotDirectory    
                        logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (alignmentResult.StringSequence |> String.filter (fun x -> x <> '*')) (alignmentResult.GlobalMod.ToString()) (alignmentResult.Charge.ToString()))
                        None
                    else
                    // let inferred_Heavy = quantifyInferredPeak getXIC identifyPeaks mz_Heavy averagePSM.WeightedAvgScanTime inferredScanTime
                    let searchRTMinusFittedRTHeavy = searchRTMinusFittedRtTarget inferredScanTime quantPHeavy
                    let clusterComparisonHeavy = comparePredictedAndMeasuredIsotopicCluster inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic inferredXicHeavy.Y_Xic_uncorrected alignmentResult.Charge labeledPeptide.BioSequence quantP.EstimatedParams.[1] mzHeavy
                    let corrLightHeavy  = calcCorrelation refinedXIC.X_Xic quantP quantPHeavy
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
                        saveChart alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.PredictedScanTime
                                    peakToQuantify.XData quantP.YPredicted quantP.YPredicted 
                                    refinedXIC.FinalTargetScanTime refinedXIC.AlignedSource
                                    successfulQuant.X_Xic successfulQuant.Y_Xic successfulQuant.xPeak successfulQuant.yFitted
                                    clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = alignmentResult.StringSequence
                    GlobalMod                                   = alignmentResult.GlobalMod
                    Charge                                      = alignmentResult.Charge
                    PepSequenceID                               = alignmentResult.PepSequenceID
                    ModSequenceID                               = alignmentResult.ModSequenceID
                    PrecursorMZ                                 = refinedXIC.FinalPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = unlabledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-unlabledPeptide.Mass)
                    MeanPercolatorScore                         = nan
                    QValue                                      = nan
                    PEPValue                                    = nan
                    ProteinNames                                = alignmentResult.ProteinNames
                    QuantMz_Light                               = refinedXIC.FinalPrecMz
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
                    QuantificationSource                        = QuantificationSource.Alignment
                    IsotopicPatternMz_Light                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Light     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Light    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Light                               = refinedXIC.X_Xic 
                    IntensityTrace_Observed_Light               = refinedXIC.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Light              = refinedXIC.Y_Xic
                    IsotopicPatternMz_Heavy                     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Heavy     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Heavy    = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Heavy                               = successfulQuant.X_Xic 
                    IntensityTrace_Observed_Heavy               = successfulQuant.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Heavy              = successfulQuant.Y_Xic
                    AlignmentScore                              = nan
                    AlignmentQValue                             = nan
                    }
                    |> Option.Some
                | None -> 
                    if diagCharts then 
                        saveChart alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.PredictedScanTime
                                    peakToQuantify.XData quantP.YPredicted quantP.YPredicted 
                                    refinedXIC.FinalTargetScanTime refinedXIC.AlignedSource
                                    [||] [||] [||] [||]
                                    clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = alignmentResult.StringSequence
                    GlobalMod                                   = alignmentResult.GlobalMod
                    Charge                                      = alignmentResult.Charge
                    PepSequenceID                               = alignmentResult.PepSequenceID
                    ModSequenceID                               = alignmentResult.ModSequenceID
                    PrecursorMZ                                 = refinedXIC.FinalPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = unlabledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-unlabledPeptide.Mass)
                    MeanPercolatorScore                         = nan
                    QValue                                      = nan
                    PEPValue                                    = nan
                    ProteinNames                                = alignmentResult.ProteinNames
                    QuantMz_Light                               = refinedXIC.FinalPrecMz
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
                    QuantificationSource                        = QuantificationSource.Alignment
                    IsotopicPatternMz_Light                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Light     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Light    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Light                               = refinedXIC.X_Xic 
                    IntensityTrace_Observed_Light               = refinedXIC.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Light              = refinedXIC.Y_Xic
                    IsotopicPatternMz_Heavy                     = [||]
                    IsotopicPatternIntensity_Observed_Heavy     = [||]
                    IsotopicPatternIntensity_Corrected_Heavy    = [||]
                    RtTrace_Heavy                               = [||]
                    IntensityTrace_Observed_Heavy               = [||]
                    IntensityTrace_Corrected_Heavy              = [||]
                    AlignmentScore                              = nan
                    AlignmentQValue                             = nan
                    }
                    |> Option.Some
            else
                let mzLight = 
                    let mz = Mass.toMZ (unlabledPeptide.Mass) (alignmentResult.Charge|> float)
                    let correctedMz = scanTimeToMzCorrection inferredScanTime + mz
                    correctedMz   
                let inferredQuant = 
                    let inferredXicLight = getInferredXic getXIC refinedXIC.FinalTargetScanTime mzLight    
                    let inferredPeaksLight = identifyPeaks inferredXicLight.X_Xic inferredXicLight.Y_Xic
                    if Array.isEmpty inferredPeaksLight then 
                        if diagCharts then saveErrorChart refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge "noInferredPeaks" plotDirectory
                        logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (alignmentResult.StringSequence |> String.filter (fun x -> x <> '*')) (alignmentResult.GlobalMod.ToString()) (alignmentResult.Charge.ToString()))
                        None
                    else
                    let peakToQuantifyLight = BioFSharp.Mz.Quantification.HULQ.getPeakBy inferredPeaksLight inferredScanTime
                    let quantPLight         = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantifyLight          
                    if quantPLight.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantPLight.EstimatedParams then 
                        if diagCharts then saveErrorChart refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge "fittingInferredFailed" plotDirectory    
                        logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (alignmentResult.StringSequence |> String.filter (fun x -> x <> '*')) (alignmentResult.GlobalMod.ToString()) (alignmentResult.Charge.ToString()))
                        None
                    else
                    // let inferred_Heavy = quantifyInferredPeak getXIC identifyPeaks mz_Heavy averagePSM.WeightedAvgScanTime inferredScanTime
                    let searchRTMinusFittedRTLight = searchRTMinusFittedRtTarget inferredScanTime quantPLight
                    let clusterComparisonLight = comparePredictedAndMeasuredIsotopicCluster inferredXicLight.X_Xic inferredXicLight.Y_Xic inferredXicLight.Y_Xic_uncorrected alignmentResult.Charge unlabledPeptide.BioSequence quantP.EstimatedParams.[1] mzLight
                    let corrLightHeavy  = calcCorrelation refinedXIC.X_Xic quantP quantPLight
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
                        saveChart alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.PredictedScanTime
                                    peakToQuantify.XData quantP.YPredicted quantP.YPredicted 
                                    refinedXIC.FinalTargetScanTime refinedXIC.AlignedSource
                                    successfulQuant.X_Xic successfulQuant.Y_Xic successfulQuant.xPeak successfulQuant.yFitted
                                    clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = alignmentResult.StringSequence
                    GlobalMod                                   = alignmentResult.GlobalMod
                    Charge                                      = alignmentResult.Charge
                    PepSequenceID                               = alignmentResult.PepSequenceID
                    ModSequenceID                               = alignmentResult.ModSequenceID
                    PrecursorMZ                                 = refinedXIC.FinalPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = labeledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-labeledPeptide.Mass)
                    MeanPercolatorScore                         = nan
                    QValue                                      = nan
                    PEPValue                                    = nan
                    ProteinNames                                = alignmentResult.ProteinNames
                    QuantMz_Light                               = mzLight
                    Quant_Light                                 = successfulQuant.Area
                    MeasuredApex_Light                          = successfulQuant.MeasuredApexIntensity
                    Seo_Light                                   = successfulQuant.StandardErrorOfPrediction
                    Params_Light                                = successfulQuant.EstimatedParams 
                    Difference_SearchRT_FittedRT_Light          = successfulQuant.SearchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Light            = successfulQuant.ClusterComparison.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Light   = successfulQuant.ClusterComparison.KLDiv_Corrected
                    QuantMz_Heavy                               = refinedXIC.FinalPrecMz
                    Quant_Heavy                                 = quantP.Area
                    MeasuredApex_Heavy                          = quantP.MeasuredApexIntensity
                    Seo_Heavy                                   = quantP.StandardErrorOfPrediction
                    Params_Heavy                                = quantP.EstimatedParams 
                    Difference_SearchRT_FittedRT_Heavy          = searchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Heavy            = clusterComparisonTarget.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterComparisonTarget.KLDiv_Corrected
                    Correlation_Light_Heavy                     = successfulQuant.Correlation_Light_Heavy
                    QuantificationSource                        = QuantificationSource.Alignment
                    IsotopicPatternMz_Light                     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Light     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Light    = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Light                               = successfulQuant.X_Xic 
                    IntensityTrace_Observed_Light               = successfulQuant.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Light              = successfulQuant.Y_Xic
                    IsotopicPatternMz_Heavy                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Heavy     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Heavy    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Heavy                               = refinedXIC.X_Xic 
                    IntensityTrace_Observed_Heavy               = refinedXIC.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Heavy              = refinedXIC.Y_Xic
                    AlignmentScore                              = nan
                    AlignmentQValue                             = nan
                    }
                    |> Option.Some
                | None -> 
                    saveChart alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.PredictedScanTime
                                peakToQuantify.XData quantP.YPredicted quantP.YPredicted 
                                refinedXIC.FinalTargetScanTime refinedXIC.AlignedSource
                                [||] [||] [||] [||]
                                clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = alignmentResult.StringSequence
                    GlobalMod                                   = alignmentResult.GlobalMod
                    Charge                                      = alignmentResult.Charge
                    PepSequenceID                               = alignmentResult.PepSequenceID
                    ModSequenceID                               = alignmentResult.ModSequenceID
                    PrecursorMZ                                 = refinedXIC.FinalPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = labeledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-labeledPeptide.Mass)
                    MeanPercolatorScore                         = nan
                    QValue                                      = nan
                    PEPValue                                    = nan
                    ProteinNames                                = alignmentResult.ProteinNames
                    QuantMz_Light                               = mzLight
                    Quant_Light                                 = nan
                    MeasuredApex_Light                          = nan
                    Seo_Light                                   = nan
                    Params_Light                                = [||] 
                    Difference_SearchRT_FittedRT_Light          = nan
                    KLDiv_Observed_Theoretical_Light            = nan
                    KLDiv_CorrectedObserved_Theoretical_Light   = nan
                    QuantMz_Heavy                               = refinedXIC.FinalPrecMz
                    Quant_Heavy                                 = quantP.Area
                    MeasuredApex_Heavy                          = quantP.MeasuredApexIntensity
                    Seo_Heavy                                   = quantP.StandardErrorOfPrediction
                    Params_Heavy                                = quantP.EstimatedParams 
                    Difference_SearchRT_FittedRT_Heavy          = searchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Heavy            = clusterComparisonTarget.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterComparisonTarget.KLDiv_Corrected
                    Correlation_Light_Heavy                     = nan
                    QuantificationSource                        = QuantificationSource.Alignment
                    IsotopicPatternMz_Light                     = [||] 
                    IsotopicPatternIntensity_Observed_Light     = [||] 
                    IsotopicPatternIntensity_Corrected_Light    = [||] 
                    RtTrace_Light                               = [||] 
                    IntensityTrace_Observed_Light               = [||] 
                    IntensityTrace_Corrected_Light              = [||] 
                    IsotopicPatternMz_Heavy                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Heavy     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Heavy    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Heavy                               = refinedXIC.X_Xic 
                    IntensityTrace_Observed_Heavy               = refinedXIC.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Heavy              = refinedXIC.Y_Xic
                    AlignmentScore                              = nan
                    AlignmentQValue                             = nan
                    }
                    |> Option.Some
            with
            | ex ->             
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None

        let lableFreeQuantification (alignmentResult:AlignmentResult) = 
            try
            if alignmentResult.GlobalMod <> 0 then 
                None
            else
                let targetPeptide  = peptideLookUp alignmentResult.StringSequence 0          
                let targetMz = Mass.toMZ (targetPeptide.Mass) (alignmentResult.Charge|> float)
                let refinedXIC = 
                    getRefinedXic processParams.PerformLocalWarp getXIC 
                        alignmentResult.RtTrace_SourceFile alignmentResult.IntensityTrace_SourceFile alignmentResult.ScanTime_SourceFile 
                            scanTimeToMzCorrection targetMz alignmentResult.PredictedScanTime
                let avgMass = Mass.ofMZ (targetMz) (alignmentResult.Charge |> float)
                let peaks = identifyPeaks refinedXIC.X_Xic refinedXIC.Y_Xic
                if Array.isEmpty peaks then 
                    if diagCharts then saveErrorChart refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge "noPeaks" plotDirectory
                    logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (alignmentResult.StringSequence |> String.filter (fun x -> x <> '*')) (alignmentResult.GlobalMod.ToString()) (alignmentResult.Charge.ToString()))
                    None
                else
                let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks refinedXIC.FinalTargetScanTime
                let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify                        
                let searchRTMinusFittedRT = searchRTMinusFittedRtTarget refinedXIC.FinalTargetScanTime quantP
                if quantP.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantP.EstimatedParams then 
                    if diagCharts then saveErrorChart refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge "fittingFailed" plotDirectory    
                    logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (alignmentResult.StringSequence |> String.filter (fun x -> x <> '*')) (alignmentResult.GlobalMod.ToString()) (alignmentResult.Charge.ToString()))
                    None
                else
                let inferredScanTime = chooseScanTime processParams.XicExtraction.ScanTimeWindow searchRTMinusFittedRT refinedXIC.FinalTargetScanTime quantP 
                let clusterComparisonTarget = comparePredictedAndMeasuredIsotopicCluster refinedXIC.X_Xic refinedXIC.Y_Xic refinedXIC.Y_Xic_uncorrected alignmentResult.Charge targetPeptide.BioSequence quantP.EstimatedParams.[1] refinedXIC.FinalPrecMz
                if diagCharts then 
                    saveChart alignmentResult.StringSequence alignmentResult.GlobalMod alignmentResult.Charge refinedXIC.X_Xic refinedXIC.Y_Xic alignmentResult.PredictedScanTime
                                peakToQuantify.XData quantP.YPredicted quantP.YPredicted 
                                refinedXIC.FinalTargetScanTime refinedXIC.AlignedSource
                                [||] [||] [||] [||]
                                clusterComparisonTarget.PeakComparisons plotDirectory
                {
                StringSequence                              = alignmentResult.StringSequence
                GlobalMod                                   = alignmentResult.GlobalMod
                Charge                                      = alignmentResult.Charge
                PepSequenceID                               = alignmentResult.PepSequenceID
                ModSequenceID                               = alignmentResult.ModSequenceID
                PrecursorMZ                                 = refinedXIC.FinalPrecMz
                MeasuredMass                                = avgMass
                TheoMass                                    = targetPeptide.Mass
                AbsDeltaMass                                = abs(avgMass-targetPeptide.Mass)
                MeanPercolatorScore                         = nan
                QValue                                      = nan
                PEPValue                                    = nan
                ProteinNames                                = alignmentResult.ProteinNames
                QuantMz_Light                               = refinedXIC.FinalPrecMz
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
                QuantificationSource                        = QuantificationSource.Alignment
                IsotopicPatternMz_Light                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                IsotopicPatternIntensity_Observed_Light     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                IsotopicPatternIntensity_Corrected_Light    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                RtTrace_Light                               = refinedXIC.X_Xic 
                IntensityTrace_Observed_Light               = refinedXIC.Y_Xic_uncorrected
                IntensityTrace_Corrected_Light              = refinedXIC.Y_Xic
                IsotopicPatternMz_Heavy                     = [||]
                IsotopicPatternIntensity_Observed_Heavy     = [||]
                IsotopicPatternIntensity_Corrected_Heavy    = [||]
                RtTrace_Heavy                               = [||]
                IntensityTrace_Observed_Heavy               = [||]
                IntensityTrace_Corrected_Heavy              = [||]
                AlignmentScore                              = nan
                AlignmentQValue                             = nan
                }
                |> Option.Some
            with
            | ex ->
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None
        logger.Trace "init quantification functions:finished"
        
        logger.Trace (sprintf "executing quantification of %i metric peptides" alignmentMetrics.Length)
        quantifyTestDataSet alignmentMetrics
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePathMetrics)
        logger.Trace "executing quantification metric peptides: finished"      
        
        logger.Trace "executing quantification"
        let quantResults = 
            alignmentResults        
            |> Array.choose (fun alignmentResult -> 
                if processParams.PerformLabeledQuantification then 
                    labledQuantification alignmentResult
                else
                    lableFreeQuantification alignmentResult
                )      
        let filteredResults = 
            if processParams.PerformLabeledQuantification then 
                quantResults
                |> lightScanTimeDifferenceFilter 4. 
                |> heavyScanTimeDifferenceFilter 4. 
                |> heavyQualityFilter -2. 2.
                |> lightQualityFilter -2. 2.          
            else
                quantResults
                |> lightScanTimeDifferenceFilter 4. 
                |> lightQualityFilter -2. 2.
        let appendedResults = 
            Array.append psmbasedQuant filteredResults       
            |> Array.groupBy (fun x -> x.StringSequence,x.ModSequenceID,x.Charge,x.GlobalMod, x.QuantificationSource)
            |> Array.choose (fun (pepIon,ions) -> 
                if processParams.PerformLabeledQuantification then
                    let min =
                        ions
                        |> Array.minBy (fun x -> 
                            if x.GlobalMod = 0 then
                                x.Difference_SearchRT_FittedRT_Light
                            else
                                x.Difference_SearchRT_FittedRT_Heavy
                            )
                    if nan.Equals min then None else Some min
                else
                    let min =
                        ions
                        |> Array.minBy (fun x -> 
                                x.Difference_SearchRT_FittedRT_Light
                            )

                    if nan.Equals min then None else Some min
                )
        appendedResults
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()
        memoryDB.Dispose()
        logger.Trace "executing quantification:finished"

        