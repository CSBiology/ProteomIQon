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

    /// The peak detection of the XIC processing as a function of the scan times and intensities. It does not catch exceptions of
    /// the peak detection. Gabor3D is not an XIC based peak detection, its function fails when it is called.
    let initIdentifyPeaks (peakDetectionParams:Domain.XicProcessing) : float [] -> float [] -> FSharp.Stats.Signal.PeakDetection.IdentifiedPeak [] =
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
        | Domain.XicProcessing.Gabor3D _ ->
            fun _ _ -> failwith "Gabor3D peak detection works on the ion mobility dimension, not on an XIC."

    /// The peak that contains the scan time, None when no peak does.
    let tryPeakContaining (peaks: FSharp.Stats.Signal.PeakDetection.IdentifiedPeak []) (scanTime: float) =
        peaks |> Array.tryFind (fun p -> p.LeftEnd.XVal < scanTime && p.RightEnd.XVal > scanTime)

    /// The peak to quantify at the scan time. With requirePeakAtScanTime only a peak that contains the scan time, none otherwise.
    /// Without it the peak that contains the scan time, else the peak with the closest apex (HULQ.getPeakBy), none only for no peaks.
    let selectPeak (requirePeakAtScanTime: bool) (peaks: FSharp.Stats.Signal.PeakDetection.IdentifiedPeak []) (scanTime: float) =
        if requirePeakAtScanTime || Array.isEmpty peaks then tryPeakContaining peaks scanTime
        else Some (HULQ.getPeakBy peaks scanTime)

    /// A fit is valid when it has at least 3 parameters, all parameters and the area are finite.
    let isValidFit (q: HULQ.QuantifiedPeak) =
        q.EstimatedParams.Length >= 3 && (q.EstimatedParams |> Array.forall Numerics.isFinite) && Numerics.isFinite q.Area

    /// Correlation of the fitted light and heavy peak over the scan times of the XIC; nan if a peak has no model, the scan time grid
    /// has fewer than 2 points or no positive spacing.
    let calcCorrelation (xValues:float []) (quantifiedPeak:HULQ.QuantifiedPeak) (inferredPeak:HULQ.QuantifiedPeak) =
        let getValue (model:HULQ.PeakModel) estParams =
            match model with
            | HULQ.PeakModel.Gaussian m -> m.GetFunctionValue (vector estParams)
            | HULQ.PeakModel.EMG m -> m.GetFunctionValue (vector estParams)
        match quantifiedPeak.Model, inferredPeak.Model with
        | Some q , Some i when xValues.Length > 1 ->
            let step = (Numerics.medianSpacing (Array.sort xValues)) / 2.
            if not (step > 0.) then nan else
            let xValues = [|xValues.[0] .. step .. xValues.[xValues.Length-1]|]
            let fQ = getValue q quantifiedPeak.EstimatedParams
            let fI = getValue i inferredPeak.EstimatedParams
            FSharp.Stats.Correlation.Seq.pearson (xValues |> Array.map fQ) (xValues |> Array.map fI)
        | _ -> nan

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
       
