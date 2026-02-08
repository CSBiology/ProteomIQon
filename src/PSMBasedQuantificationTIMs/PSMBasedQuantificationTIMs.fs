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
open SearchDB'
open Newtonsoft.Json.Linq

module PSMBasedQuantificationTIMs =
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
        let initRTProfile (readspecPeaks:string -> Peak1DArray)  (rtIndex: IMzIOArray<RtIndexEntry>) (rtRange: RangeQuery) (mzRange: RangeQuery) (ionMobilityRange: RangeQuery)=
            let entries = RtIndexEntry.Search(rtIndex, rtRange).ToArray()
            let profile = Array.zeroCreate<Peak2D> entries.Length
            let mzLow = mzRange.LowValue
            let mzHigh = mzRange.HighValue
            let mzLock = mzRange.LockValue
            let imLow = ionMobilityRange.LowValue
            let imHigh = ionMobilityRange.HighValue
            let imLock = ionMobilityRange.LockValue
            for rtIdx = 0 to entries.Length - 1 do
                let entry = entries.[rtIdx]
                let peaks = (readspecPeaks entry.SpectrumID).Peaks
                let mutable bestMz = 0.
                let mutable bestDelta = System.Double.PositiveInfinity
                let mutable bestIntensity = 0.
                for peak in RtIndexEntry.MzSearch(peaks, mzRange) do
                    if peak.Mz < mzHigh && peak.Mz > mzLow then
                        let ionMobility = peak.IonMobility.Value
                        if ionMobility < imHigh && ionMobility > imLow then
                            let delta = abs (peak.Mz - mzLock)
                            if delta < bestDelta then
                                bestDelta <- delta
                                bestMz <- peak.Mz
                                bestIntensity <- peak.Intensity
                            elif peak.Mz = bestMz then
                                bestIntensity <- bestIntensity + peak.Intensity
                if System.Double.IsPositiveInfinity bestDelta then
                    profile.[rtIdx] <- new Peak2D(0., mzLock, entry.Rt, imLock)
                else
                    profile.[rtIdx] <- new Peak2D(bestIntensity, bestMz, entry.Rt)
            profile

        /// Like `initRTProfile`, but emits retention time and intensity arrays directly.
        /// Avoids `Peak2D` object materialization in the quantification hot path.
        let initRTIntensityProfile (readspecPeaks:string -> Peak1DArray) (rtIndex: IMzIOArray<RtIndexEntry>) (rtRange: RangeQuery) (mzRange: RangeQuery) (ionMobilityRange: RangeQuery) =
            let entries = RtIndexEntry.Search(rtIndex, rtRange).ToArray()
            let rtData = Array.zeroCreate<float> entries.Length
            let intensityData = Array.zeroCreate<float> entries.Length
            let mzLow = mzRange.LowValue
            let mzHigh = mzRange.HighValue
            let mzLock = mzRange.LockValue
            let imLow = ionMobilityRange.LowValue
            let imHigh = ionMobilityRange.HighValue
            for rtIdx = 0 to entries.Length - 1 do
                let entry = entries.[rtIdx]
                rtData.[rtIdx] <- entry.Rt
                let peaks = (readspecPeaks entry.SpectrumID).Peaks
                let mutable bestMz = 0.
                let mutable bestDelta = System.Double.PositiveInfinity
                let mutable bestIntensity = 0.
                for peak in RtIndexEntry.MzSearch(peaks, mzRange) do
                    if peak.Mz < mzHigh && peak.Mz > mzLow then
                        let ionMobility = peak.IonMobility.Value
                        if ionMobility < imHigh && ionMobility > imLow then
                            let delta = abs (peak.Mz - mzLock)
                            if delta < bestDelta then
                                bestDelta <- delta
                                bestMz <- peak.Mz
                                bestIntensity <- peak.Intensity
                            elif peak.Mz = bestMz then
                                bestIntensity <- bestIntensity + peak.Intensity
                if System.Double.IsPositiveInfinity bestDelta then
                    intensityData.[rtIdx] <- 0.
                else
                    intensityData.[rtIdx] <- bestIntensity
            rtData, intensityData

    type RuntimeTuning =
        {
            PeptideDbMode            : string option
            PeakCacheMode            : string option
            PeakCacheMax             : int option
            PeakCacheMaxPeaks        : int64 option
            FragmentMatchMode        : string option
            IsotopicPatternCacheMode : string option
            RtIndexMode              : string option
        }

    module RuntimeTuning =
        let defaultValue =
            {
                PeptideDbMode            = None
                PeakCacheMode            = None
                PeakCacheMax             = None
                PeakCacheMaxPeaks        = None
                FragmentMatchMode        = None
                IsotopicPatternCacheMode = None
                RtIndexMode              = None
            }

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
    let getClosestMs1 (ms1s: (float*string) []) scanTime = 
        if Array.isEmpty ms1s then
            invalidArg "ms1s" "ms1 index must not be empty"
        elif scanTime <= fst ms1s.[0] then
            snd ms1s.[0]
        elif scanTime >= fst ms1s.[ms1s.Length - 1] then
            snd ms1s.[ms1s.Length - 1]
        else
            let mutable lo = 0
            let mutable hi = ms1s.Length - 1
            while hi - lo > 1 do
                let mid = lo + ((hi - lo) >>> 1)
                if fst ms1s.[mid] <= scanTime then
                    lo <- mid
                else
                    hi <- mid
            let rtLo,idLo = ms1s.[lo]
            let rtHi,idHi = ms1s.[hi]
            if abs (scanTime - rtLo) <= abs (rtHi - scanTime) then idLo else idHi

    let private tryParseFirstValueAsInt (token: JToken) =
        if obj.ReferenceEquals(token, null) then
            None
        else
            let valuesToken = token.["Values"]
            if obj.ReferenceEquals(valuesToken, null) || valuesToken.Type <> JTokenType.Array then
                None
            else
                let values = valuesToken :?> JArray
                if values.Count = 0 then
                    None
                else
                    let raw = values.[0].ToString()
                    match System.Int32.TryParse(raw, System.Globalization.NumberStyles.Any, System.Globalization.CultureInfo.InvariantCulture) with
                    | true, parsed -> Some parsed
                    | _ -> None

    let private tryParseFirstValueAsFloat (token: JToken) =
        if obj.ReferenceEquals(token, null) then
            None
        else
            let valuesToken = token.["Values"]
            if obj.ReferenceEquals(valuesToken, null) || valuesToken.Type <> JTokenType.Array then
                None
            else
                let values = valuesToken :?> JArray
                if values.Count = 0 then
                    None
                else
                    let raw = values.[0].ToString()
                    match System.Double.TryParse(raw, System.Globalization.NumberStyles.Any, System.Globalization.CultureInfo.InvariantCulture) with
                    | true, parsed -> Some parsed
                    | _ -> None

    let private tryGetMsLevelAndScanTime (description: string) =
        try
            let descriptionJson = JObject.Parse(description)
            let msLevel =
                descriptionJson.SelectToken("properties['MS:1000511']")
                |> tryParseFirstValueAsInt
            let scansPropertiesToken = descriptionJson.SelectToken("Scans.properties")
            let scanTime =
                if obj.ReferenceEquals(scansPropertiesToken, null) || scansPropertiesToken.Type <> JTokenType.Object then
                    None
                else
                    let scansProperties = scansPropertiesToken :?> JObject
                    scansProperties.Properties()
                    |> Seq.tryPick (fun scanProp ->
                        let scanStartToken = scanProp.Value.SelectToken("properties['MS:1000016']")
                        tryParseFirstValueAsFloat scanStartToken
                    )
            match msLevel, scanTime with
            | Some level, Some rt -> Some (level, rt)
            | _ -> None
        with
        | _ -> None

    let private buildMs1RtIndexStreamSql (logger: NLog.Logger) (reader: MzIO.MzSQL.MzSQL) (runID: string) =
        use cmd = new SQLiteCommand("SELECT SpectrumID, Description FROM Spectrum WHERE RunID = @runID", reader.Connection)
        cmd.Parameters.AddWithValue("@runID", runID) |> ignore
        use sqlReader = cmd.ExecuteReader()
        let entries = ResizeArray<MzIO.Processing.MzIOLinq.RtIndexEntry>()
        let mutable totalRows = 0
        let mutable ms1Rows = 0
        let mutable skippedRows = 0
        while sqlReader.Read() do
            totalRows <- totalRows + 1
            let spectrumID = sqlReader.GetString(0)
            let description = sqlReader.GetString(1)
            match tryGetMsLevelAndScanTime description with
            | Some (level, rt) when level = 1 ->
                entries.Add(new MzIO.Processing.MzIOLinq.RtIndexEntry(rt, spectrumID))
                ms1Rows <- ms1Rows + 1
            | _ ->
                skippedRows <- skippedRows + 1
        entries.Sort(fun a b -> a.Rt.CompareTo(b.Rt))
        logger.Trace (sprintf "RT index stream_sql: rows=%d, ms1Entries=%d, skipped=%d" totalRows ms1Rows skippedRows)
        MzIO.Commons.Arrays.MzIOArray.ToMzIOArray(entries)

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
    
    /// Calculates the Kullback-Leibler divergence Dkl(p||q) from q (theory, model, description, or approximation of p) 
    /// to p (the "true" distribution of data, observations, or a precisely calculated theoretical distribution).
    let klDiv (p:float []) (q:float []) = 
        Array.fold2 (fun acc p q -> (System.Math.Log(p/q)*p) + acc ) 0. p q
     
    ///
    let substractBaseLine (logger: NLog.Logger) (baseLineParams:Domain.BaseLineCorrection) (yData:float []) =
        if yData.Length > 500 then
            yData
        elif yData |> Array.forall (fun y -> y = 0.) then
            yData
        else
            let baseLine = FSharp.Stats.Signal.Baseline.baselineAls' baseLineParams.MaxIterations baseLineParams.Lambda baseLineParams.P yData |> Array.ofSeq
            Array.map2 (fun y b ->
                           let c = y - b
                           if c < 0. then 0. else c
                       ) yData baseLine

    ///
    let initGetProcessedXIC logger (baseLineCorrection:Domain.BaseLineCorrection option) getRtIntensity idx scanTimeWindow mzWindow_Da imWindow meanScanTime meanPrecMz meanIM =
        let rtQuery = Query.createRangeQuery meanScanTime scanTimeWindow
        let mzQuery = Query.createRangeQuery meanPrecMz mzWindow_Da
        let imQuery = Query.createRangeQuery meanIM imWindow
        let retData',itzData' =
            let (rtRaw: float[]), (intensityRaw: float[]) = getRtIntensity idx rtQuery mzQuery imQuery
            let n = rtRaw.Length
            let rtData = Array.zeroCreate<float> n
            let itzData = Array.zeroCreate<float> n
            let mutable outCount = 0
            for i = 0 to n - 1 do
                let intensity = intensityRaw.[i]
                let keep =
                    if i = 0 || i = n - 1 || intensity > 0. then
                        true
                    else
                        let prevIntensity = intensityRaw.[i-1]
                        if prevIntensity = 0. then
                            true
                        elif prevIntensity > (100. * (intensity + 1.)) then
                            false
                        else
                            true
                if keep then
                    rtData.[outCount] <- rtRaw.[i]
                    itzData.[outCount] <- intensity
                    outCount <- outCount + 1
            if outCount = n then
                rtData, itzData
            else
                rtData.[.. outCount - 1], itzData.[.. outCount - 1]
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

    ///
    let weightedMean (weights:seq<'T>) (items:seq<'T>) =
        let sum,n = Seq.fold2 (fun (sum,n) w i -> w*i+sum,n + w ) (0.,0.) weights items
        sum / n
        
    ///
    let average getXic scanTimeToMzCorrection theoMz (psms:(PSMStatisticsResultFragpipe*float) []) =
            //let meanPrecMz   = psms |> Seq.meanBy (fun (psm,m) -> psm.PrecursorMZ)
            
            let meanScanTime = psms |> Seq.meanBy (fun (psm,m) -> psm.ScanTime)
            let meanScore = psms |> Seq.averageBy (fun (psm,m) -> psm.Hyperscore)
            let psms' = 
                let tmp = Array.sortByDescending (fun (psm,m) -> m) psms
                if tmp.Length > 3 then tmp.[..2] else tmp 
            let weightedAvgScanTime =
                let scanTimes = psms' |> Array.map (fun (psm,m) -> psm.ScanTime)
                let weights = psms' |> Array.map snd
                weightedMean weights scanTimes
            let weightedAvgIM =
                let ims = psms' |> Array.map (fun (psm,m) -> psm.IonMobility)
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

    /// Calculates the difference between the scan time used to retreave the peak and the fitted peak midpoint.
    let searchRTMinusFittedRtInferred searchRT fit = 
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
    let saveChart sequence globalMod ch (xXic:float[]) (yXic:float[]) ms2s avgScanTime (xToQuantify:float[]) (ypToQuantify:float[]) (fitY:float[])
            (xXicInferred:float[]) (yXicinferred:float[]) (xInferred:float[]) (inferredFit:float[]) (*(xEnvelopeSum:float[]) (yEnvelopeSum:float[])*) (peaks:FSharp.Stats.Signal.PeakDetection.IdentifiedPeak []) (pattern:PeakComparison []) plotDirectory =
        let xic = 
            [
            Chart.Point(xXic, yXic)                     |> Chart.withTraceName "Target XIC"
            peaks
            |> Array.map (fun x -> Chart.Point(x.XData,x.YData)) 
            |> Chart.Combine
            Chart.Point(ms2s)                           |> Chart.withTraceName "MS2s with scores"
            Chart.Point([avgScanTime],[1.])             |> Chart.withTraceName "Weighted Mean of Ms2 scan times"
            Chart.Point((xToQuantify), (ypToQuantify))  |> Chart.withTraceName "Identified Target Peak"
            Chart.Line(xToQuantify,fitY)                |> Chart.withTraceName "Fit of target Peak"
            Chart.Point(xXicInferred, yXicinferred)     |> Chart.withTraceName "Inferred XIC"
            Chart.Line(xInferred,inferredFit)           |> Chart.withTraceName "Fit of inferred Peak"
            //Chart.Point(xEnvelopeSum, yEnvelopeSum)     |> Chart.withTraceName "Target Envelope Sum"
            ]
            |> Chart.Combine
        let pattern = 
            [
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensity))          |> Chart.withTraceName "Measured"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensityCorrected)) |> Chart.withTraceName "Measured Corrected"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.PredictedRelFrequency))      |> Chart.withTraceName "Predicted Relative Frequency"
            //Chart.Point(xEnvelopeSum, yEnvelopeSum)     |> Chart.withTraceName "Target Envelope Sum"
            ]
            |> Chart.Combine
        [xic;pattern]
        |> Chart.Stack(2, 0.1)
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i" sequence globalMod)
        |> Chart.withSize(2500.,800.)
        |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + globalMod.ToString() + "Ch" + ch.ToString())|])

    let saveErrorChart (xXic:float[]) (yXic:float[]) pepIon desc plotDirectory =      
        Chart.Point(xXic, yXic)
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i_%s" pepIon.Sequence pepIon.GlobalMod desc)
        |> Chart.withSize(1500.,800.)
        |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((pepIon.Sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + pepIon.GlobalMod.ToString() + "Ch" + pepIon.Charge.ToString() + "_notQuantified")|])
            
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
    let initComparePredictedAndMeasuredIsotopicCluster
        (getPredictedIsotopicPattern: int -> AminoAcids.AminoAcid list -> (float * float)[])
        (readSpecPeaks: string -> MzIO.Binary.Peak1DArray)
        ms1s
        ms1AccuracyEstimate
        (x_Xic:float[])
        (y_Xic:float[])
        y_Xic_uncorrected
        ch
        peptideSequence
        tarRt
        tarMz =
        /// IsotopicCluster
        let targetIsotopicPattern_predicted = 
            getPredictedIsotopicPattern ch peptideSequence
        let baseLineCorrectionF = getBaseLineCorrectionOffsetAt tarRt x_Xic y_Xic y_Xic_uncorrected
        let closestMS1SpectrumID = getClosestMs1 ms1s tarRt
        let peaksInWindow: ResizeArray<MzIO.Binary.Peak1D> =
            let peaks = (readSpecPeaks closestMS1SpectrumID).Peaks
            let lowMz = tarMz - 0.6
            let highMz = tarMz + 1.
            let candidates = ResizeArray<MzIO.Binary.Peak1D>()
            for i = 0 to peaks.Length - 1 do
                let peak = peaks.[i]
                if peak.Mz > lowMz && peak.Mz < highMz then
                    candidates.Add peak
            candidates
        let recordedVsPredictedPattern = 
            if peaksInWindow.Count = 0 then
                [||]
            else
                let groupedRelFreq = System.Collections.Generic.Dictionary<int, float>()
                let groupedOrder = ResizeArray<int>()
                for i = 0 to targetIsotopicPattern_predicted.Length - 1 do
                    let mz,relFreq = targetIsotopicPattern_predicted.[i]
                    let mutable bestIdx = -1
                    let mutable bestDelta = System.Double.PositiveInfinity
                    for pIdx = 0 to peaksInWindow.Count - 1 do
                        let peak = peaksInWindow.[pIdx]
                        let delta = abs (peak.Mz - mz)
                        if delta < bestDelta then
                            bestDelta <- delta
                            bestIdx <- pIdx
                    if bestIdx >= 0 && bestDelta < 4. * ms1AccuracyEstimate then
                        let mutable currentRelFreq = 0.
                        if groupedRelFreq.TryGetValue(bestIdx, &currentRelFreq) then
                            groupedRelFreq.[bestIdx] <- currentRelFreq + relFreq
                        else
                            groupedRelFreq.[bestIdx] <- relFreq
                            groupedOrder.Add(bestIdx)
                groupedOrder
                |> Seq.map (fun peakIdx ->
                    let peak = peaksInWindow.[peakIdx]
                    let predictedRelFrequency = groupedRelFreq.[peakIdx]
                    let mz = peak.Mz
                    let measuredIntensity = peak.Intensity
                    {
                        Mz = mz
                        MeasuredIntensity = measuredIntensity
                        MeasuredIntensityCorrected = measuredIntensity - baseLineCorrectionF
                        PredictedRelFrequency = predictedRelFrequency
                    }
                )
                |> Seq.filter (fun isoP -> isoP.MeasuredIntensityCorrected > 0.)
                |> Array.ofSeq
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
    let quantifyPeptides diagCharts zipCharts (processParams:Domain.QuantificationParams) (runtimeTuning:RuntimeTuning) (outputDir:string) (cn:SQLiteConnection) (instrumentOutput:string) (scoredPSMs:string)  =
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension scoredPSMs)
        let swOverall = System.Diagnostics.Stopwatch.StartNew()
        let swSetup = System.Diagnostics.Stopwatch.StartNew()
        let mutable stageCountMatchedMassesMs = 0L
        let mutable stageXicExtractionMs = 0L
        let mutable stageIsotopicCompareMs = 0L
        let mutable readSpectrumPeaksDbReads = 0L
        let mutable peakSpectrumCacheHits = 0L
        let mutable peakSpectrumCacheMisses = 0L
        let mutable peakSpectrumCacheEvictions = 0L
        let mutable peakSpectrumCacheResidentPeaks = 0L
        let mutable peakSpectrumCacheResidentSpectra = 0L
        let mutable isotopicPatternCacheHits = 0L
        let mutable isotopicPatternCacheMisses = 0L
        let getRuntimeSetting settingValue =
            match settingValue with
            | Some value when not (System.String.IsNullOrWhiteSpace value) -> Some value
            | _ -> None
        let parsePositiveInt settingName fallback settingValue =
            match settingValue with
            | Some value when value > 0 -> value
            | Some value ->
                logger.Trace (sprintf "Invalid runtime %s=%d. Using default %d." settingName value fallback)
                fallback
            | None -> fallback
        let parsePositiveInt64Opt settingName settingValue =
            match settingValue with
            | Some value when value > 0L -> Some value
            | Some value ->
                logger.Trace (sprintf "Invalid runtime %s=%d. Using unbounded peak count." settingName value)
                None
            | None -> None
        logger.Trace (sprintf "Input file: %s" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        logger.Trace (sprintf "RuntimeTuning: %A" runtimeTuning)
        logger.Trace (sprintf "Now performing Quantification using: %s and %s, Results will be written to: %s" instrumentOutput scoredPSMs outputDir)
        // initialize Reader and Transaction
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".quant"
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
        logger.Trace (sprintf "plotDirectory:%s" plotDirectory)
        let peptideDbMode = getRuntimeSetting runtimeTuning.PeptideDbMode
        let useMemoryDbCopy =
            match peptideDbMode with
            | Some mode ->
                System.String.Equals(mode, "memory", System.StringComparison.OrdinalIgnoreCase)
                || System.String.Equals(mode, "memory_copy", System.StringComparison.OrdinalIgnoreCase)
            | None -> false
        let peptideDB =
            if useMemoryDbCopy then
                logger.Trace "Peptide DB mode: memory_copy"
                logger.Trace "Copy peptide DB into Memory"
                let memoryDb = SearchDB.copyDBIntoMemory cn
                logger.Trace "Copy peptide DB into Memory: finished"
                memoryDb
            else
                logger.Trace "Peptide DB mode: file_backed"
                let fileDb = new SQLiteConnection(cn.ConnectionString)
                fileDb.Open()
                fileDb
        logger.Trace "Get peptide lookUp function"
        let dBParams     = getSDBParams peptideDB
        //let massLookUp = prepareSelectMassByModSequenceAndGlobalMod peptideDB
        let peptideLookUp = getThreadSafePeptideLookUpFromFileBySequenceAndGMod peptideDB dBParams
        let calcIonSeries aal = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal
        logger.Trace "Get peptide lookUp function: finished"
        // initialize Reader and Transaction
        logger.Trace "Init connection to mass spectrum data."
        let inReader = Core.MzIO.Reader.getReader instrumentOutput :?> MzIO.MzSQL.MzSQL
        inReader.Connection.Open()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader       
        let inTr = inReader.BeginTransaction()
        let readSpectrumPeaksCounted spectrumID =
            System.Threading.Interlocked.Increment(&readSpectrumPeaksDbReads) |> ignore
            inReader.ReadSpectrumPeaks spectrumID
        let peakSpectrumCacheMax = parsePositiveInt "runtime-peak-cache-max" 700 runtimeTuning.PeakCacheMax
        let peakSpectrumCacheMaxPeaksOpt = parsePositiveInt64Opt "runtime-peak-cache-max-peaks" runtimeTuning.PeakCacheMaxPeaks
        let peakSpectrumCacheMaxPeaksLabel =
            match peakSpectrumCacheMaxPeaksOpt with
            | Some maxPeaks -> maxPeaks.ToString()
            | None -> "unbounded"
        let peakSpectrumCacheMode = getRuntimeSetting runtimeTuning.PeakCacheMode
        let readSpecPeaksWithMem, peakSpectrumCacheModeUsed =
            if
                match peakSpectrumCacheMode with
                | Some mode -> System.String.Equals(mode, "unbounded", System.StringComparison.OrdinalIgnoreCase)
                | None -> false
            then
                logger.Trace "Peak spectrum cache mode: unbounded memoize"
                FSharpAux.Memoization.memoize readSpectrumPeaksCounted, "unbounded"
            else
                logger.Trace (sprintf "Peak spectrum cache mode: bounded_lru (maxSpectra=%d, maxPeaks=%s)" peakSpectrumCacheMax peakSpectrumCacheMaxPeaksLabel)
                let lruList = System.Collections.Generic.LinkedList<string * MzIO.Binary.Peak1DArray>()
                let cacheMap = System.Collections.Generic.Dictionary<string, System.Collections.Generic.LinkedListNode<string * MzIO.Binary.Peak1DArray>>(System.StringComparer.Ordinal)
                let readSpecPeaksWithLru spectrumID =
                    let mutable node = Unchecked.defaultof<System.Collections.Generic.LinkedListNode<string * MzIO.Binary.Peak1DArray>>
                    if cacheMap.TryGetValue(spectrumID, &node) then
                        peakSpectrumCacheHits <- peakSpectrumCacheHits + 1L
                        lruList.Remove(node)
                        lruList.AddFirst(node)
                        snd node.Value
                    else
                        peakSpectrumCacheMisses <- peakSpectrumCacheMisses + 1L
                        let loaded = readSpectrumPeaksCounted spectrumID
                        peakSpectrumCacheResidentPeaks <- peakSpectrumCacheResidentPeaks + (loaded.Peaks.Length |> int64)
                        let newNode = System.Collections.Generic.LinkedListNode<string * MzIO.Binary.Peak1DArray>((spectrumID, loaded))
                        lruList.AddFirst(newNode)
                        cacheMap.[spectrumID] <- newNode
                        let shouldEvictByPeaks =
                            match peakSpectrumCacheMaxPeaksOpt with
                            | Some maxPeaks -> peakSpectrumCacheResidentPeaks > maxPeaks
                            | None -> false
                        let mutable shouldEvict = cacheMap.Count > peakSpectrumCacheMax || shouldEvictByPeaks
                        while shouldEvict do
                            let tailNode = lruList.Last
                            if not (obj.ReferenceEquals(tailNode, null)) then
                                lruList.RemoveLast()
                                peakSpectrumCacheResidentPeaks <- peakSpectrumCacheResidentPeaks - ((snd tailNode.Value).Peaks.Length |> int64)
                                cacheMap.Remove(fst tailNode.Value) |> ignore
                                peakSpectrumCacheEvictions <- peakSpectrumCacheEvictions + 1L
                            let shouldEvictByPeaks' =
                                match peakSpectrumCacheMaxPeaksOpt with
                                | Some maxPeaks -> peakSpectrumCacheResidentPeaks > maxPeaks
                                | None -> false
                            shouldEvict <- cacheMap.Count > peakSpectrumCacheMax || shouldEvictByPeaks'
                        peakSpectrumCacheResidentSpectra <- cacheMap.Count |> int64
                        loaded
                readSpecPeaksWithLru, "bounded_lru"
        let fragmentMatchMode = getRuntimeSetting runtimeTuning.FragmentMatchMode
        let useLegacyFragmentMatch =
            match fragmentMatchMode with
            | Some mode -> System.String.Equals(mode, "legacy", System.StringComparison.OrdinalIgnoreCase)
            | None -> false
        if useLegacyFragmentMatch then
            logger.Trace "Fragment match mode: legacy"
        else
            logger.Trace "Fragment match mode: optimized_binary_search"
        let isotopicPatternCacheMode = getRuntimeSetting runtimeTuning.IsotopicPatternCacheMode
        let useIsotopicPatternCache =
            match isotopicPatternCacheMode with
            | Some mode -> System.String.Equals(mode, "on", System.StringComparison.OrdinalIgnoreCase)
            | None -> false
        if useIsotopicPatternCache then
            logger.Trace "Isotopic pattern cache mode: enabled"
        else
            logger.Trace "Isotopic pattern cache mode: disabled"
        let isotopicPatternCache: System.Collections.Generic.Dictionary<int * (AminoAcids.AminoAcid list), (float * float)[]> =
            System.Collections.Generic.Dictionary<int * (AminoAcids.AminoAcid list), (float * float)[]>()
        let getPredictedIsotopicPattern (ch:int) (peptideSequence:AminoAcids.AminoAcid list) =
            if useIsotopicPatternCache then
                let key = ch, peptideSequence
                let mutable cached = Unchecked.defaultof<(float * float)[]>
                if isotopicPatternCache.TryGetValue(key, &cached) then
                    isotopicPatternCacheHits <- isotopicPatternCacheHits + 1L
                    cached
                else
                    let generated = generateIsotopicDistributionOfFormulaBySum ch peptideSequence
                    isotopicPatternCacheMisses <- isotopicPatternCacheMisses + 1L
                    isotopicPatternCache.Add(key, generated)
                    generated
            else
                let generated = generateIsotopicDistributionOfFormulaBySum ch peptideSequence
                isotopicPatternCacheMisses <- isotopicPatternCacheMisses + 1L
                generated
        let hasFragmentMatchWithinPpm100 (sortedFragments: float[]) peakMz =
            if sortedFragments.Length = 0 then
                false
            else
                let tolerance = Mass.deltaMassByPpm 100. peakMz
                let searchResult = System.Array.BinarySearch(sortedFragments, peakMz)
                if searchResult >= 0 then
                    true
                else
                    let insertionIdx = ~~~searchResult
                    let mutable leftIdx = insertionIdx - 1
                    let mutable rightIdx = insertionIdx
                    let mutable found = false
                    while not found && leftIdx >= 0 && peakMz - sortedFragments.[leftIdx] <= tolerance do
                        if abs (sortedFragments.[leftIdx] - peakMz) <= tolerance then
                            found <- true
                        leftIdx <- leftIdx - 1
                    while not found && rightIdx < sortedFragments.Length && sortedFragments.[rightIdx] - peakMz <= tolerance do
                        if abs (sortedFragments.[rightIdx] - peakMz) <= tolerance then
                            found <- true
                        rightIdx <- rightIdx + 1
                    found
        let countMatchedMasses (peptide: LookUpResult<AminoAcids.AminoAcid>)(psms: PSMStatisticsResultFragpipe []) =
            let sw = System.Diagnostics.Stopwatch.StartNew()
            let result =
                let fragList = 
                    let ionSeries = (calcIonSeries peptide.BioSequence).TargetMasses
                    [1. .. 2.]
                    |> List.collect (fun ch -> 
                        ionSeries 
                        |> List.map (fun x -> Mass.toMZ x.MainPeak.Mass ch)
                    )
                psms
                |> Array.map (fun psm -> 
                    let spec = readSpecPeaksWithMem psm.PSMId
                    let sum = 
                        if useLegacyFragmentMatch then
                            spec.Peaks 
                            |> Seq.filter (fun peak -> 
                                fragList
                                |> List.exists (fun ion -> abs (ion - peak.Mz) <= (Mass.deltaMassByPpm 100. peak.Mz))
                            )
                            |> Seq.sumBy (fun x -> x.Intensity)
                        else
                            let frag = fragList |> Array.ofList
                            Array.sortInPlace frag
                            let peaks = spec.Peaks
                            let mutable intensitySum = 0.
                            for i = 0 to peaks.Length - 1 do
                                let peak = peaks.[i]
                                if hasFragmentMatchWithinPpm100 frag peak.Mz then
                                    intensitySum <- intensitySum + peak.Intensity
                            intensitySum
                    psm,sum
                )
            sw.Stop()
            stageCountMatchedMassesMs <- stageCountMatchedMassesMs + sw.ElapsedMilliseconds
            result
        logger.Trace "Create RetentionTime index"
        let rtIndexMode =
            getRuntimeSetting runtimeTuning.RtIndexMode
            |> Option.defaultValue "stream_sql"
        let retTimeIdxed =
            if System.String.Equals(rtIndexMode, "legacy", System.StringComparison.OrdinalIgnoreCase) then
                logger.Trace "RT index mode: legacy"
                Query.getMS1RTIdx inReader inRunID
            else
                logger.Trace "RT index mode: stream_sql"
                try
                    let streamed = buildMs1RtIndexStreamSql logger inReader inRunID
                    if streamed |> Seq.isEmpty |> not then
                        streamed
                    else
                        logger.Trace "RT index stream_sql returned no entries. Falling back to legacy builder."
                        Query.getMS1RTIdx inReader inRunID
                with
                | ex ->
                    logger.Trace (sprintf "RT index stream_sql failed (%s). Falling back to legacy builder." ex.Message)
                    Query.getMS1RTIdx inReader inRunID
        logger.Trace "Create RetentionTime index:finished"
        
        logger.Trace "Read and sort ms1s"
        /// Build sorted ms1 index directly from rt index entries to avoid a second full metadata read.
        let ms1SortedByScanTime =
            retTimeIdxed
            |> Seq.map (fun entry -> entry.Rt, entry.SpectrumID)
            |> Array.ofSeq
        logger.Trace "Read and sort ms1s:finished"
        
        logger.Trace "Read scored PSMs"
        ///
        let qpsms =
            Csv.CsvReader<PSMStatisticsResultFragpipe>(SchemaMode=Csv.Fill).ReadFile(scoredPSMs,'\t',false,1)
            |> Array.ofSeq
        logger.Trace "Read scored PSMs:finished"
        
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
            stDev, scanTimeToMzCorrection 
        logger.Trace (sprintf "Estimate precursor mz standard deviation and mz correction.:finished, standard deviation: %f" ms1AccuracyEstimate) 

        let qpsmsMzRefined = 
            let refined = 
                qpsms
                |> Array.map (fun x -> 
                    let precMz = x.PrecursorMZ
                    let theoMz  = Mass.toMZ x.TheoMass (float x.Charge)
                    let diffPToTheo = precMz - theoMz 
                    let absDiffToPredDiff = 
                        (abs (scanTimeToMzCorrection x.ScanTime)) - (abs diffPToTheo)
                        |> abs
                    if absDiffToPredDiff > 4.*ms1AccuracyEstimate then 
                        let precMz' = theoMz + scanTimeToMzCorrection x.ScanTime
                        {x with PrecursorMZ = precMz'}
                    else    
                        x
                    )
            if diagCharts then 
                [
                qpsms
                |> Array.map (fun x -> 
                    let precMz = x.PrecursorMZ
                    let theoMz  = Mass.toMZ x.TheoMass (float x.Charge)
                    let diff = precMz - theoMz 
                    x.ScanTime,diff
                    )
                |> Chart.Point
                |> Chart.withTraceName "Raw"
                refined
                |> Array.map (fun x -> 
                    let precMz = x.PrecursorMZ
                    let theoMz  = Mass.toMZ x.TheoMass (float x.Charge)
                    let diff = precMz - theoMz 
                    x.ScanTime,diff
                    )
                |> Chart.Point
                |> Chart.withTraceName "Corrected"
                ]
                |> Chart.Combine
                |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; "precMzCorrected"|])
            refined
                
           
        logger.Trace "init lookup functions"        
        ///
        let comparePredictedAndMeasuredIsotopicClusterRaw = initComparePredictedAndMeasuredIsotopicCluster getPredictedIsotopicPattern readSpecPeaksWithMem ms1SortedByScanTime ms1AccuracyEstimate
        let comparePredictedAndMeasuredIsotopicCluster xXic yXic yXicUncorrected ch peptideSequence tarRt tarMz =
            let sw = System.Diagnostics.Stopwatch.StartNew()
            let result = comparePredictedAndMeasuredIsotopicClusterRaw xXic yXic yXicUncorrected ch peptideSequence tarRt tarMz
            sw.Stop()
            stageIsotopicCompareMs <- stageIsotopicCompareMs + sw.ElapsedMilliseconds
            result
        
        let mzWindow = 
            match processParams.XicExtraction.MzWindow_Da with 
            | Domain.Window.Fixed v -> v 
            | Domain.Window.Estimate -> 
                let mzW = ms1AccuracyEstimate*4.
                logger.Trace (sprintf "optimal mz Window for XIC look up found by estimation :%f Da" mzW)  
                mzW

        ///
        let getRtIntensity = Query.initRTIntensityProfile readSpecPeaksWithMem
        let getXICRaw = initGetProcessedXIC logger processParams.BaseLineCorrection getRtIntensity retTimeIdxed processParams.XicExtraction.ScanTimeWindow mzWindow 0.05
        let getXIC meanScanTime meanPrecMz meanIM =
            let sw = System.Diagnostics.Stopwatch.StartNew()
            let result = getXICRaw meanScanTime meanPrecMz meanIM
            sw.Stop()
            stageXicExtractionMs <- stageXicExtractionMs + sw.ElapsedMilliseconds
            result
        ///
        let identifyPeaks = initIdentifyPeaks processParams.XicExtraction.XicProcessing
        logger.Trace "init lookup functions:finished"
        
        logger.Trace "init quantification functions"
        ///
        let labledQuantification (pepIon:PeptideIon) (psms:PSMStatisticsResultFragpipe []) = 
            try
            let bestQValue,prots = psms |> Array.minBy (fun x -> x.Expectscore) |> fun x -> x.Expectscore,x.ProteinNames
            let unlabledPeptide = peptideLookUp pepIon.Sequence 0
            let labeledPeptide  = peptideLookUp pepIon.Sequence 1
            let targetPeptide = if pepIon.GlobalMod = 0 then unlabledPeptide else labeledPeptide            
            let psmsWithMatchedSums = countMatchedMasses targetPeptide psms 
            let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
            let theoMz = Mass.toMZ targetPeptide.Mass (float pepIon.Charge)
            let averagePSM = average getXIC scanTimeToMzCorrection theoMz psmsWithMatchedSums
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
            let clusterComparisonTarget = comparePredictedAndMeasuredIsotopicCluster averagePSM.X_Xic averagePSM.Y_Xic averagePSM.Y_Xic_uncorrected pepIon.Charge targetPeptide.BioSequence quantP.EstimatedParams.[1] averagePSM.MeanPrecMz            
            if pepIon.GlobalMod = 0 then
                let mzHeavy = 
                    let mz = Mass.toMZ (labeledPeptide.Mass) (pepIon.Charge|> float)
                    let correctedMz = scanTimeToMzCorrection inferredScanTime + mz
                    correctedMz 
                let inferredQuant = 
                    let inferredXicHeavy = getInferredXic getXIC averagePSM.WeightedAvgScanTime mzHeavy averagePSM.WeightedAvgIM
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
                    let clusterComparisonHeavy = comparePredictedAndMeasuredIsotopicCluster inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic inferredXicHeavy.Y_Xic_uncorrected pepIon.Charge labeledPeptide.BioSequence quantP.EstimatedParams.[1] mzHeavy
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
                    PEPValue                                    = nan
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
                    AlignmentQValue                             = nan
                    IonMobility                                 = averagePSM.WeightedAvgIM
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
                    PEPValue                                    = nan
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
                    AlignmentQValue                             = nan
                    IonMobility                                 = averagePSM.WeightedAvgIM
                    }
                    |> Option.Some
            else
                let mzLight = 
                    let mz = Mass.toMZ (unlabledPeptide.Mass) (pepIon.Charge|> float)
                    let correctedMz = scanTimeToMzCorrection inferredScanTime + mz
                    correctedMz   
                let inferredQuant = 
                    let inferredXicLight = getInferredXic getXIC averagePSM.WeightedAvgScanTime mzLight averagePSM.WeightedAvgIM
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
                    let clusterComparisonLight = comparePredictedAndMeasuredIsotopicCluster inferredXicLight.X_Xic inferredXicLight.Y_Xic inferredXicLight.Y_Xic_uncorrected pepIon.Charge unlabledPeptide.BioSequence quantP.EstimatedParams.[1] mzLight
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
                    PEPValue                                    = nan
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
                    AlignmentQValue                             = nan
                    IonMobility                                 = averagePSM.WeightedAvgIM
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
                    PEPValue                                    = nan
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
                    AlignmentQValue                             = nan
                    IonMobility                                 = averagePSM.WeightedAvgIM
                    }
                    |> Option.Some
            with
            | ex ->             
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None

        let lableFreeQuantification (pepIon:PeptideIon)  (psms:PSMStatisticsResultFragpipe []) = 
            try
            if pepIon.GlobalMod <> 0 then None 
            else
            let bestQValue,prots = psms |> Array.minBy (fun x -> x.Expectscore) |> fun x -> x.Expectscore,x.ProteinNames
            let unlabledPeptide = peptideLookUp pepIon.Sequence 0
            let psmsWithMatchedSums = countMatchedMasses unlabledPeptide psms 
            let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
            let theoMz = Mass.toMZ unlabledPeptide.Mass (float pepIon.Charge)
            let averagePSM = average getXIC scanTimeToMzCorrection theoMz psmsWithMatchedSums
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
            let clusterComparisonTarget = comparePredictedAndMeasuredIsotopicCluster averagePSM.X_Xic averagePSM.Y_Xic averagePSM.Y_Xic_uncorrected pepIon.Charge unlabledPeptide.BioSequence quantP.EstimatedParams.[1] averagePSM.MeanPrecMz
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
            PEPValue                                    = nan
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
            AlignmentQValue                             = nan
            IonMobility                                 = averagePSM.WeightedAvgIM
            }
            |> Option.Some
            with
            | ex ->
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None
        logger.Trace "init quantification functions:finished"
        swSetup.Stop()
        logger.Trace (sprintf "Perf: setupMs=%d" swSetup.ElapsedMilliseconds)
        
        logger.Trace "executing quantification"
        let swQuantLoop = System.Diagnostics.Stopwatch.StartNew()
        let quantResults = 
            qpsmsMzRefined 
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
                        |> Array.sortByDescending (fun x -> x.Hyperscore)
                        |> Array.take x
                    | _ -> 
                        pepIon,
                        psms                  
                )
            |> Array.mapi (fun i (pepIon,psms) -> 
                if i % 100 = 0 then
                    let proc = System.Diagnostics.Process.GetCurrentProcess()
                    let workingSetMB = proc.WorkingSet64 / (1024L * 1024L)
                    let privateBytesMB = proc.PrivateMemorySize64 / (1024L * 1024L)
                    let managedHeapMB = System.GC.GetTotalMemory(false) / (1024L * 1024L)
                    let gc0 = System.GC.CollectionCount 0
                    let gc1 = System.GC.CollectionCount 1
                    let gc2 = System.GC.CollectionCount 2
                    logger.Trace (sprintf "%i peptides quantified" i)
                    logger.Trace (
                        sprintf
                            "PerfProgress peptides=%d quantLoopMs=%d countMatchedMassesMs=%d xicExtractionMs=%d isotopicCompareMs=%d readSpectrumPeaksDbReads=%d peakCacheHits=%d peakCacheMisses=%d peakCacheEvictions=%d peakCacheResidentSpectra=%d peakCacheResidentPeaks=%d managedHeapMB=%d workingSetMB=%d privateBytesMB=%d gc0=%d gc1=%d gc2=%d"
                            i
                            swQuantLoop.ElapsedMilliseconds
                            stageCountMatchedMassesMs
                            stageXicExtractionMs
                            stageIsotopicCompareMs
                            readSpectrumPeaksDbReads
                            peakSpectrumCacheHits
                            peakSpectrumCacheMisses
                            peakSpectrumCacheEvictions
                            peakSpectrumCacheResidentSpectra
                            peakSpectrumCacheResidentPeaks
                            managedHeapMB
                            workingSetMB
                            privateBytesMB
                            gc0
                            gc1
                            gc2
                    )
                
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
        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()
        peptideDB.Dispose()
        swOverall.Stop()
        logger.Trace (
            sprintf
                "PerfSummary totalMs=%d setupMs=%d countMatchedMassesMs=%d xicExtractionMs=%d isotopicCompareMs=%d readSpectrumPeaksDbReads=%d peakCacheResidentSpectra=%d peakCacheResidentPeaks=%d managedHeapMB=%d"
                swOverall.ElapsedMilliseconds
                swSetup.ElapsedMilliseconds
                stageCountMatchedMassesMs
                stageXicExtractionMs
                stageIsotopicCompareMs
                readSpectrumPeaksDbReads
                peakSpectrumCacheResidentSpectra
                peakSpectrumCacheResidentPeaks
                (System.GC.GetTotalMemory(false) / (1024L * 1024L))
        )
        logger.Trace (
            sprintf
                "PerfCaches peakSpectrumCacheMode=%s peakSpectrumCacheMax=%d peakSpectrumCacheMaxPeaks=%s peakSpectrumCacheHits=%d peakSpectrumCacheMisses=%d peakSpectrumCacheEvictions=%d isotopicPatternCacheHits=%d isotopicPatternCacheMisses=%d isotopicPatternCacheCount=%d"
                peakSpectrumCacheModeUsed
                peakSpectrumCacheMax
                peakSpectrumCacheMaxPeaksLabel
                peakSpectrumCacheHits
                peakSpectrumCacheMisses
                peakSpectrumCacheEvictions
                isotopicPatternCacheHits
                isotopicPatternCacheMisses
                isotopicPatternCache.Count
        )
        logger.Trace "executing quantification:finished"
        
