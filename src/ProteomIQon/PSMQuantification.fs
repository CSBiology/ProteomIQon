namespace ProteomIQon

open System
open System.IO
open FSharp.Stats
open Plotly.NET
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
            
