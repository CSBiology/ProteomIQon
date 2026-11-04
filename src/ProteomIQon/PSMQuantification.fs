namespace ProteomIQon

open System
open System.IO
open FSharp.Stats
open Plotly.NET
open BioFSharp
open BioFSharp.Mz
open BioFSharp.Mz.Quantification
open BioFSharp.Mz.SearchDB
open Dto
open PeakQuantification

/// Quantification of the PSMs of one peptide ion: the shared types and helpers of PSMBasedQuantification and PSMBasedQuantificationTIMs.
module PSMQuantification =

    type PeptideIon =
        {
            Sequence             : string
            GlobalMod            : int
            Charge               : int
            ModSequenceID        : int
            PepSequenceID        : int
        }
            
    type AveragePSM = {
        MeanPrecMz   : float
        MeanScanTime : float
        WeightedAvgScanTime:float
        WeightedAvgIM:float
        MeanScore   : float
        X_Xic         : float []
        Y_Xic         : float []
        Y_Xic_uncorrected: float []
        }

    ///
    let createAveragePSM meanPrecMz meanScanTime weightedAvgScanTime weightedAvgIM meanScore xXic yXic yXic_uncorrected = {
        MeanPrecMz    = meanPrecMz
        MeanScanTime  = meanScanTime
        WeightedAvgScanTime= weightedAvgScanTime
        WeightedAvgIM = weightedAvgIM
        MeanScore = meanScore
        X_Xic         = xXic
        Y_Xic         = yXic
        Y_Xic_uncorrected= yXic_uncorrected
        }

    ///
    let weightedMean (weights:seq<float>) (items:seq<float>) =
        let sum,n = Seq.fold2 (fun (sum,n) w i -> w*i+sum,n + w ) (0.,0.) weights items
        sum / n
        
    /// The PSM with the lowest QValue. A PSM without QValue (NaN) is taken only when no PSM has one.
    let bestPsmByQValue (psms: PSMStatisticsResult []) =
        match psms |> Array.filter (fun x -> not (Double.IsNaN x.QValue)) with
        | [||] -> psms.[0]
        | withQValue -> withQValue |> Array.minBy (fun x -> x.QValue)

    ///
    /// ionMobilityOf gives the ion mobility of the MS2 spectrum of a PSM by its PSMId. The top 3 PSMs are taken
    /// among the PSMs whose spectrum has an ion mobility, among all when none has one.
    let average (ionMobilityOf: string -> float) getXic scanTimeToMzCorrection theoMz (psms:(PSMStatisticsResult*float) []) =
            //let meanPrecMz   = psms |> Seq.meanBy (fun (psm,m) -> psm.PrecursorMZ)

            let meanScanTime = psms |> Seq.meanBy (fun (psm,m) -> psm.ScanTime)
            let meanScore = psms |> Seq.averageBy (fun (psm,m) -> psm.ModelScore)
            let psms' =
                let withIonMobility = psms |> Array.filter (fun (psm,m) -> Numerics.isFinite (ionMobilityOf psm.PSMId))
                let tmp = Array.sortByDescending (fun (psm,m) -> m) (if Array.isEmpty withIonMobility then psms else withIonMobility)
                if tmp.Length > 3 then tmp.[..2] else tmp 
            let weightedAvgScanTime =
                let scanTimes = psms' |> Array.map (fun (psm,m) -> psm.ScanTime)
                let weights = psms' |> Array.map snd
                weightedMean weights scanTimes
            let weightedAvgIM =
                let ims = psms' |> Array.map (fun (psm,m) -> ionMobilityOf psm.PSMId)
                let weights = psms' |> Array.map snd
                weightedMean weights ims
            let correctedMz = scanTimeToMzCorrection weightedAvgScanTime + theoMz
            let (retData,itzDataCorrected,ItzDataUncorrected) = getXic weightedAvgScanTime correctedMz weightedAvgIM
            createAveragePSM correctedMz meanScanTime weightedAvgScanTime weightedAvgIM meanScore retData itzDataCorrected ItzDataUncorrected


    type InferredXic = {
        X_Xic                       :float[]
        Y_Xic                       :float[]
        Y_Xic_uncorrected           :float[]
        }
    
    let getInferredXic getXic targetScanTime targetMz targetIM =
        let (retData,itzData,uncorrectedItzData)   =
                getXic targetScanTime targetMz targetIM
        {
        X_Xic               = retData
        Y_Xic               = itzData
        Y_Xic_uncorrected   = uncorrectedItzData
        }

    ///
    let saveChart sequence globalMod ch (xXic:float[]) (yXic:float[]) ms2s avgScanTime (xToQuantify:float[]) (ypToQuantify:float[]) (fitY:float[])
            (xXicInferred:float[]) (yXicinferred:float[]) (xInferred:float[]) (inferredFit:float[]) (*(xEnvelopeSum:float[]) (yEnvelopeSum:float[])*) (peaks:FSharp.Stats.Signal.PeakDetection.IdentifiedPeak []) (pattern:PeakComparison []) plotDirectory =
        let xic = 
            [
            Chart.Point(xXic, yXic)                     |> Chart.withTraceInfo "Target XIC"
            peaks
            |> Array.map (fun x -> Chart.Point(x.XData,x.YData)) 
            |> Chart.combine
            Chart.Point(ms2s)                           |> Chart.withTraceInfo "MS2s with scores"
            Chart.Point([avgScanTime],[1.])             |> Chart.withTraceInfo "Weighted Mean of Ms2 scan times"
            Chart.Point((xToQuantify), (ypToQuantify))  |> Chart.withTraceInfo "Identified Target Peak"
            Chart.Line(xToQuantify,fitY)                |> Chart.withTraceInfo "Fit of target Peak"
            Chart.Point(xXicInferred, yXicinferred)     |> Chart.withTraceInfo "Inferred XIC"
            Chart.Line(xInferred,inferredFit)           |> Chart.withTraceInfo "Fit of inferred Peak"
            //Chart.Point(xEnvelopeSum, yEnvelopeSum)     |> Chart.withTraceInfo "Target Envelope Sum"
            ]
            |> Chart.combine
        let pattern = 
            [
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensity))          |> Chart.withTraceInfo "Measured"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensityCorrected)) |> Chart.withTraceInfo "Measured Corrected"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.PredictedRelFrequency))      |> Chart.withTraceInfo "Predicted Relative Frequency"
            //Chart.Point(xEnvelopeSum, yEnvelopeSum)     |> Chart.withTraceInfo "Target Envelope Sum"
            ]
            |> Chart.combine
        [xic;pattern]
        |> Chart.Grid(1, 2, XGap = 0.1, YGap = 0.1)
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i" sequence globalMod)
        |> Chart.withSize(2500,800)
        |> Chart.saveHtml(Path.Combine[|plotDirectory; ((sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + globalMod.ToString() + "Ch" + ch.ToString())|])

    let saveErrorChart (xXic:float[]) (yXic:float[]) pepIon desc plotDirectory =      
        Chart.Point(xXic, yXic)
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i_%s" pepIon.Sequence pepIon.GlobalMod desc)
        |> Chart.withSize(1500,800)
        |> Chart.saveHtml(Path.Combine[|plotDirectory; ((pepIon.Sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + pepIon.GlobalMod.ToString() + "Ch" + pepIon.Charge.ToString() + "_notQuantified")|])
            

    type private InferredQuantification = {
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

    /// The functions the quantification of a peptide ion uses; they depend on the run (reader, calibration, parameters).
    type QuantificationFunctions =
        {
            Logger                                      : NLog.Logger
            DiagCharts                                  : bool
            PlotDirectory                               : string
            /// Half width of the scan time window of the XIC (min).
            ScanTimeWindow                              : float
            PeptideLookUp                               : string -> int -> LookUpResult<AminoAcids.AminoAcid>
            /// The PSMs with the summed intensity of their fragment ions that match the peptide.
            CountMatchedMasses                          : LookUpResult<AminoAcids.AminoAcid> -> PSMStatisticsResult [] -> (PSMStatisticsResult * float) []
            ScanTimeToMzCorrection                      : float -> float
            /// The ion mobility of the MS2 spectrum of a PSM by its PSMId, nan without ion mobility.
            IonMobilityOf                               : string -> float
            /// Scan time, m/z and ion mobility to the scan times, corrected and uncorrected intensities of the XIC.
            GetXic                                      : float -> float -> float -> float [] * float [] * float []
            IdentifyPeaks                               : float [] -> float [] -> FSharp.Stats.Signal.PeakDetection.IdentifiedPeak []
            ComparePredictedAndMeasuredIsotopicCluster  : float [] -> float [] -> float [] -> int -> LookUpResult<AminoAcids.AminoAcid> -> float -> float -> ClusterComparison
        }

    /// Quantifies the light and heavy ion of a peptide ion from its PSMs; None when it fails.
    let quantifyLabeled (f: QuantificationFunctions) (pepIon:PeptideIon) (psms:PSMStatisticsResult []) = 
        // a failing peak detection of the partner leaves the partner without quantification, the quantification of the ion stays
        let identifyPartnerPeaks xData yData =
            try
                f.IdentifyPeaks xData yData
            with ex ->
                f.Logger.Trace (sprintf "Quant failed: Peak detection of the labeled partner failed with: %A" ex)
                [||]
        try
        let bestQValue,bestPepValue,prots = bestPsmByQValue psms |> fun x -> x.QValue, x.PEPValue,x.ProteinNames
        let unlabledPeptide = f.PeptideLookUp pepIon.Sequence 0
        let labeledPeptide  = f.PeptideLookUp pepIon.Sequence 1
        let targetPeptide = if pepIon.GlobalMod = 0 then unlabledPeptide else labeledPeptide            
        let psmsWithMatchedSums = f.CountMatchedMasses targetPeptide psms 
        let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
        let theoMz = Mass.toMZ targetPeptide.Mass (float pepIon.Charge)
        let averagePSM = average f.IonMobilityOf f.GetXic f.ScanTimeToMzCorrection theoMz psmsWithMatchedSums
        let avgMass = Mass.ofMZ (averagePSM.MeanPrecMz) (pepIon.Charge |> float)
        let peaks = f.IdentifyPeaks averagePSM.X_Xic averagePSM.Y_Xic
        if Array.isEmpty peaks then 
            if f.DiagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "noPeaks" f.PlotDirectory
            f.Logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
            None
        else
        let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks averagePSM.WeightedAvgScanTime
        let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify                        
        let searchRTMinusFittedRT = searchRTMinusFittedRtTarget averagePSM.WeightedAvgScanTime quantP
        if not (isValidFit quantP) then 
            if f.DiagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "fittingFailed" f.PlotDirectory    
            f.Logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
            None
        else
        let inferredScanTime = chooseScanTime f.ScanTimeWindow searchRTMinusFittedRT averagePSM.WeightedAvgScanTime quantP 
        let clusterComparisonTarget = f.ComparePredictedAndMeasuredIsotopicCluster averagePSM.X_Xic averagePSM.Y_Xic averagePSM.Y_Xic_uncorrected pepIon.Charge targetPeptide quantP.EstimatedParams.[1] averagePSM.MeanPrecMz            
        if pepIon.GlobalMod = 0 then
            let mzHeavy = 
                let mz = Mass.toMZ (labeledPeptide.Mass) (pepIon.Charge|> float)
                let correctedMz = f.ScanTimeToMzCorrection inferredScanTime + mz
                correctedMz 
            let inferredQuant = 
                let inferredXicHeavy = getInferredXic f.GetXic averagePSM.WeightedAvgScanTime mzHeavy averagePSM.WeightedAvgIM
                let inferredPeaksHeavy = identifyPartnerPeaks inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic
                if Array.isEmpty inferredPeaksHeavy then 
                    if f.DiagCharts then saveErrorChart inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic pepIon "noInferredPeaks" f.PlotDirectory
                    f.Logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                    None
                else
                let peakToQuantifyHeavy = BioFSharp.Mz.Quantification.HULQ.getPeakBy inferredPeaksHeavy inferredScanTime
                let quantPHeavy         = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantifyHeavy          
                if not (isValidFit quantPHeavy) then 
                    if f.DiagCharts then saveErrorChart inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic pepIon "fittingInferredFailed" f.PlotDirectory    
                    f.Logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                    None
                else
                let searchRTMinusFittedRTHeavy = searchRTMinusFittedRtTarget inferredScanTime quantPHeavy
                let clusterComparisonHeavy = f.ComparePredictedAndMeasuredIsotopicCluster inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic inferredXicHeavy.Y_Xic_uncorrected pepIon.Charge labeledPeptide quantP.EstimatedParams.[1] mzHeavy
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
                if f.DiagCharts then 
                    saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                    peakToQuantify.XData peakToQuantify.YData quantP.YPredicted successfulQuant.X_Xic successfulQuant.Y_Xic successfulQuant.xPeak successfulQuant.yFitted peaks clusterComparisonTarget.PeakComparisons f.PlotDirectory
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
                IonMobility                                 = averagePSM.WeightedAvgIM
                }
                |> Option.Some
            | None -> 
                if f.DiagCharts then 
                    saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                    peakToQuantify.XData peakToQuantify.YData quantP.YPredicted [||] [||] [||] [||] peaks clusterComparisonTarget.PeakComparisons f.PlotDirectory
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
                IonMobility                                 = averagePSM.WeightedAvgIM
                }
                |> Option.Some
        else
            let mzLight = 
                let mz = Mass.toMZ (unlabledPeptide.Mass) (pepIon.Charge|> float)
                let correctedMz = f.ScanTimeToMzCorrection inferredScanTime + mz
                correctedMz   
            let inferredQuant = 
                let inferredXicLight = getInferredXic f.GetXic averagePSM.WeightedAvgScanTime mzLight averagePSM.WeightedAvgIM
                let inferredPeaksLight = identifyPartnerPeaks inferredXicLight.X_Xic inferredXicLight.Y_Xic
                if Array.isEmpty inferredPeaksLight then 
                    if f.DiagCharts then saveErrorChart inferredXicLight.X_Xic inferredXicLight.Y_Xic pepIon "noInferredPeaks" f.PlotDirectory
                    f.Logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                    None
                else
                let peakToQuantifyLight = BioFSharp.Mz.Quantification.HULQ.getPeakBy inferredPeaksLight inferredScanTime
                let quantPLight         = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantifyLight          
                if not (isValidFit quantPLight) then 
                    if f.DiagCharts then saveErrorChart inferredXicLight.X_Xic inferredXicLight.Y_Xic pepIon "fittingInferredFailed" f.PlotDirectory    
                    f.Logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                    None
                else
                let searchRTMinusFittedRTLight = searchRTMinusFittedRtTarget inferredScanTime quantPLight
                let clusterComparisonLight = f.ComparePredictedAndMeasuredIsotopicCluster inferredXicLight.X_Xic inferredXicLight.Y_Xic inferredXicLight.Y_Xic_uncorrected pepIon.Charge unlabledPeptide quantP.EstimatedParams.[1] mzLight
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
                if f.DiagCharts then 
                    saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                    peakToQuantify.XData peakToQuantify.YData quantP.YPredicted successfulQuant.X_Xic successfulQuant.Y_Xic successfulQuant.xPeak successfulQuant.yFitted peaks clusterComparisonTarget.PeakComparisons f.PlotDirectory
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
                IonMobility                                 = averagePSM.WeightedAvgIM
                }
                |> Option.Some
            | None ->
                if f.DiagCharts then 
                    saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                        peakToQuantify.XData peakToQuantify.YData quantP.YPredicted [||] [||] [||] [||] peaks clusterComparisonTarget.PeakComparisons f.PlotDirectory
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
                IonMobility                                 = averagePSM.WeightedAvgIM
                }
                |> Option.Some
        with
        | ex ->             
            f.Logger.Trace (sprintf "Quantfailed: %A" ex)
            Option.None

    /// Quantifies a peptide ion without a labeled partner from its PSMs; None when it fails.
    let quantifyLabelFree (f: QuantificationFunctions) (pepIon:PeptideIon)  (psms:PSMStatisticsResult []) = 
        try
        if pepIon.GlobalMod <> 0 then None 
        else
        let bestQValue,bestPepValue,prots = bestPsmByQValue psms |> fun x -> x.QValue, x.PEPValue,x.ProteinNames
        let unlabledPeptide = f.PeptideLookUp pepIon.Sequence 0
        let psmsWithMatchedSums = f.CountMatchedMasses unlabledPeptide psms 
        let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
        let theoMz = Mass.toMZ unlabledPeptide.Mass (float pepIon.Charge)
        let averagePSM = average f.IonMobilityOf f.GetXic f.ScanTimeToMzCorrection theoMz psmsWithMatchedSums
        let avgMass = Mass.ofMZ (averagePSM.MeanPrecMz) (pepIon.Charge |> float)      
        let peaks = f.IdentifyPeaks averagePSM.X_Xic averagePSM.Y_Xic               
        if Array.isEmpty peaks then 
            if f.DiagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "noPeaks" f.PlotDirectory
            f.Logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
            None
        else
        let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks averagePSM.WeightedAvgScanTime
        let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify
        let searchRTMinusFittedRT = searchRTMinusFittedRtTarget averagePSM.WeightedAvgScanTime quantP 
        if not (isValidFit quantP) then 
            if f.DiagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "fittingFailed" f.PlotDirectory    
            f.Logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
            None
        else
        let clusterComparisonTarget = f.ComparePredictedAndMeasuredIsotopicCluster averagePSM.X_Xic averagePSM.Y_Xic averagePSM.Y_Xic_uncorrected pepIon.Charge unlabledPeptide quantP.EstimatedParams.[1] averagePSM.MeanPrecMz
        if f.DiagCharts then 
            saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                peakToQuantify.XData peakToQuantify.YData quantP.YPredicted  [||] [||] [||] [||] peaks clusterComparisonTarget.PeakComparisons f.PlotDirectory
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
        IonMobility                                 = averagePSM.WeightedAvgIM
        }
        |> Option.Some
        with
        | ex ->
            f.Logger.Trace (sprintf "Quantfailed: %A" ex)
            Option.None
