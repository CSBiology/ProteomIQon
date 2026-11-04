namespace ProteomIQon

open FSharp.Stats
open BioFSharp
open BioFSharp.Mz
open BioFSharp.Mz.Quantification
open BioFSharp.Mz.SearchDB

/// Peak detection, peak choice and isotope pattern comparison of the XIC based quantification.
module PeakQuantification =

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
    let getBaseLineCorrectionOffsetAt tarRT x_Xic y_Xic y_Xic_uncorrected =
        let (rt,y,yUncorr) = 
            Array.zip3 x_Xic y_Xic y_Xic_uncorrected 
            |> Array.minBy (fun (rt,y,yUncorr) -> abs (rt - tarRT))
        yUncorr - y
       
    /// Calculates the Kullback-Leibler divergence Dkl(p||q) from q (theory, model, description, or approximation of p) 
    /// to p (the "true" distribution of data, observations, or a precisely calculated theoretical distribution).
    let klDiv (p:float []) (q:float []) = 
        Array.fold2 (fun acc p q -> (System.Math.Log(p/q)*p) + acc ) 0. p q
     

    /// Calculates the difference between the scan time used to retreave the peak and the fitted peak midpoint.
    let searchRTMinusFittedRtTarget searchRT (fit:HULQ.QuantifiedPeak) = 
        try
            searchRT - fit.EstimatedParams.[1] 
        with
        | _ -> nan


    /// Aims to provide the best scan time estimate when quanitfying a inferred peak.
    let chooseScanTime (maxDiff: float) (searchRTMinusFittedRT: float) initialScanTime (quantP:HULQ.QuantifiedPeak) = 
        if Array.isEmpty quantP.EstimatedParams then
            initialScanTime 
        elif abs searchRTMinusFittedRT > maxDiff then
            initialScanTime 
        else
            quantP.EstimatedParams.[1]

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
                try
                let  windowSize = getWindowWidth yData
                FSharp.Stats.Signal.PeakDetection.SecondDerivative.getPeaks parameters.MinSNR parameters.PolynomOrder windowSize xData yData
                with 
                | ex ->
                    // logger.Trace (sprintf "Quant failed: Peak detection failed with: %A" ex)
                    [||]
                )
        | Domain.XicProcessing.Wavelet parameters ->
            (fun xData yData -> 
                try
                FSharpStats'.Wavelet.identify parameters xData yData
                with 
                | ex ->
                    // logger.Trace (sprintf "Quant failed: Peak detection failed with: %A" ex)
                    [||]
                )
        | Domain.XicProcessing.Gabor3D _ ->
            fun _ _ -> [||]
        
    ///
    /// getMs1Peaks gives the peaks of a MS1 spectrum (by ID) with lowMz < m/z < highMz in file order.
    let initComparePredictedAndMeasuredIsotopicCluster (isotopePatterns: IsotopePatterns.IsotopePatternLookUp) (getMs1Peaks: string -> float -> float -> Peak[]) (getClosestMs1: float -> string) ms1AccuracyEstimate (x_Xic:float[]) (y_Xic:float[]) y_Xic_uncorrected ch (peptide: LookUpResult<AminoAcids.AminoAcid>) tarRt tarMz =    
        /// IsotopicCluster
        let targetIsotopicPattern_predicted = 
            isotopePatterns.Get (peptide, ch)
        let baseLineCorrectionF = getBaseLineCorrectionOffsetAt tarRt x_Xic y_Xic y_Xic_uncorrected
        let closestMS1 = getClosestMs1 tarRt
        let peaks' = 
            getMs1Peaks closestMS1 (tarMz - 0.6) (tarMz + 1.)
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
       
