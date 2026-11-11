namespace ProteomIQon

open System
open System.IO
open System.Collections.Generic
open FSharpAux.IO.SchemaReader
open FSharp.Stats
open Plotly.NET
open BioFSharp
open BioFSharp.Mz
open BioFSharp.Mz.Quantification
open MzIO.Processing
open Dto
open Dto.QuantificationResult
open Xic
open PeakQuantification


/// MS1 access: every spectrum is read once, its peaks sorted by m/z, and the points of all XICs whose scan time window contains it are taken from it.
module XicSweep =

    open Numerics

    /// How the points of the XIC of a request are taken.
    type XicPoint =
        /// The peak closest to the lock m/z (Xic.closestPoint).
        | ClosestPeak
        /// The sum of the peaks inside the m/z and ion mobility window (Xic.extractPoint).
        | Sum
        /// The same sum; the mobility profile of every trace is recorded with it (Xic.extractPointWithProfile).
        | SumWithProfiles

    /// The XIC traces of one ion: one trace per lock m/z, all with the same scan time and ion mobility window.
    type XicRequest =
        {
            RtLow       : float
            RtHigh      : float
            LockMz      : float []
            MzWindow    : float
            /// nan without ion mobility window.
            ImLow       : float
            ImHigh      : float
            Point       : XicPoint
        }

    /// One point per MS1 spectrum in the scan time window; arrays per trace are indexed [trace][point].
    type RawXic =
        {
            Rt          : float []
            Intensity   : float [] []
            ObservedMz  : float [] []
            IonMobility : float [] []
        }

    type private Accumulator(request: XicRequest) =
        let nTraces = request.LockMz.Length
        member val Rt          = ResizeArray<float>()
        member val Intensity   = Array.init nTraces (fun _ -> ResizeArray<float>())
        member val ObservedMz  = Array.init nTraces (fun _ -> ResizeArray<float>())
        member val IonMobility = Array.init nTraces (fun _ -> ResizeArray<float>())
        /// One mobility profile per trace if the request records them.
        member val Profiles    =
            match request.Point with
            | SumWithProfiles -> Some (Array.init nTraces (fun _ -> Xic.MobilityProfile()))
            | _ -> None
        member this.ToXic() =
            {
                Rt          = this.Rt.ToArray()
                Intensity   = this.Intensity |> Array.map (fun a -> a.ToArray())
                ObservedMz  = this.ObservedMz |> Array.map (fun a -> a.ToArray())
                IonMobility = this.IonMobility |> Array.map (fun a -> a.ToArray())
            }

    /// Walks once through the MS1 spectra (sorted by scan time), reads only spectra inside at least one window
    /// and hands every XIC to onComplete as soon as its window is passed. Only the current spectrum and the XICs
    /// of open windows are held in memory. Every request takes its points by the rule in its Point field. onComplete gets the profiles of the traces of
    /// a request that records them, None otherwise. Returns the number of spectra read.
    let sweep (readSpectrum: string -> SortedPeaks.MzSortedPeaks) (ms1: (float*string) []) (requests: XicRequest []) (onComplete: int -> RawXic -> Xic.MobilityProfile [] option -> unit) =
        let n = requests.Length
        let order =
            Array.init n id
            |> Array.filter (fun r -> isFinite requests.[r].RtLow && isFinite requests.[r].RtHigh)
            |> Array.sortBy (fun r -> requests.[r].RtLow)
        let accumulators = Dictionary<int,Accumulator>()
        let complete r =
            onComplete r (accumulators.[r].ToXic()) accumulators.[r].Profiles
            accumulators.Remove r |> ignore
        let mutable next = 0
        let mutable spectraRead = 0
        for (rt, spectrumID) in ms1 do
            while next < order.Length && requests.[order.[next]].RtLow <= rt do
                accumulators.[order.[next]] <- Accumulator(requests.[order.[next]])
                next <- next + 1
            let finished = accumulators.Keys |> Seq.filter (fun r -> requests.[r].RtHigh < rt) |> Array.ofSeq
            finished |> Array.sort |> Array.iter complete
            if accumulators.Count > 0 then
                let spectrum = readSpectrum spectrumID
                spectraRead <- spectraRead + 1
                for KeyValue(r, acc) in accumulators do
                    let request = requests.[r]
                    acc.Rt.Add rt
                    for t in 0 .. request.LockMz.Length - 1 do
                        let intensity, observedMz, ionMobility =
                            match request.Point with
                            | ClosestPeak -> Xic.closestPoint spectrum request.LockMz.[t] request.MzWindow request.ImLow request.ImHigh
                            | Sum -> Xic.extractPoint spectrum request.LockMz.[t] request.MzWindow request.ImLow request.ImHigh
                            | SumWithProfiles -> Xic.extractPointWithProfile acc.Profiles.Value.[t] spectrum request.LockMz.[t] request.MzWindow request.ImLow request.ImHigh
                        acc.Intensity.[t].Add intensity
                        acc.ObservedMz.[t].Add observedMz
                        acc.IonMobility.[t].Add ionMobility
        // windows after the last spectrum and the ones still open
        while next < order.Length do
            accumulators.[order.[next]] <- Accumulator(requests.[order.[next]])
            next <- next + 1
        accumulators.Keys |> Array.ofSeq |> Array.sort |> Array.iter complete
        spectraRead

    /// The raw XIC of an ion summed around the ion mobility centre of the ion, and that centre. The centre is Xic.tryMobilityCentre of
    /// the monoisotopic target trace at the scan time, around the expected ion mobility. With a centre every trace is the sum of the
    /// pairs of its profile within halfWidth around it (Xic.summedTrace: intensity and ion mobility of a point); the observed m/z of a
    /// point stays that of the point of the sweep, except where the summed intensity is 0, where it is nan because the m/z of the
    /// sweep comes from signal outside the summed window. Without a centre the XIC is returned unchanged.
    let sumAroundMobilityCentre (expected: float) (halfWidth: float) (scanTime: float) (raw: RawXic) (profiles: Xic.MobilityProfile []) =
        match Xic.tryMobilityCentre profiles.[0] raw.Rt raw.Intensity.[0] scanTime expected halfWidth with
        | None -> raw, None
        | Some centre ->
            let traces = profiles |> Array.map (fun profile -> Xic.summedTrace profile centre halfWidth)
            let intensity = traces |> Array.map fst
            { raw with
                Intensity   = intensity
                ObservedMz  = Array.map2 (fun (mz: float []) (y: float []) -> Array.map2 (fun m v -> if v > 0. then m else nan) mz y) raw.ObservedMz intensity
                IonMobility = traces |> Array.map snd },
            Some centre


/// Processing of the raw XIC and candidate peaks with their selection features.
module PeakSelection =

    open Numerics

    /// The traces of the ion (or of its labeled partner) after removing drop-out points; index 0 is the monoisotopic trace.
    type TraceSet =
        {
            Rt          : float []
            /// Baseline corrected monoisotopic trace.
            Corrected   : float []
            /// Raw intensity per isotope trace.
            Intensity   : float [] []
            ObservedMz  : float [] []
            IonMobility : float [] []
        }

    /// The traces with the given indices of the raw XIC; the first index is the monoisotopic trace.
    let createTraceSet (baseline: Domain.BaseLineCorrection option) (keepIonMobility: bool) (raw: XicSweep.RawXic) (traces: int []) =
        let keep = keepMask raw.Intensity.[traces.[0]]
        let pick (a: float []) = [| for i in 0 .. a.Length - 1 do if keep.[i] then yield a.[i] |]
        let intensity = traces |> Array.map (fun t -> pick raw.Intensity.[t])
        {
            Rt          = pick raw.Rt
            Corrected   = subtractBaseline baseline intensity.[0]
            Intensity   = intensity
            ObservedMz  = traces |> Array.map (fun t -> pick raw.ObservedMz.[t])
            // only the monoisotopic ion mobility of the target is read afterwards
            IonMobility = if keepIonMobility then [| pick raw.IonMobility.[traces.[0]] |] else [||]
        }

    /// Names of the candidate features, in the order of Candidate.Features.
    let featureNames = [|"RtZ"; "RtZ2"; "ShapeCorr"; "LogRelHeight"; "AbsLog2FwhmRatio"; "Log2NCandidates"; "IsoCorr"; "IsoRatioLog2"; "AbsPpmMono"; "ImZ"|]
    let [<Literal>] RtZ              = 0
    let [<Literal>] RtZ2             = 1
    let [<Literal>] ShapeCorr        = 2
    let [<Literal>] LogRelHeight     = 3
    let [<Literal>] AbsLog2FwhmRatio = 4
    let [<Literal>] Log2NCandidates  = 5
    let [<Literal>] IsoCorr          = 6
    let [<Literal>] IsoRatioLog2     = 7
    let [<Literal>] AbsPpmMono       = 8
    let [<Literal>] ImZ              = 9

    /// A peak of the monoisotopic XIC that may be the ion.
    type Candidate =
        {
            Peak     : FSharp.Stats.Signal.PeakDetection.IdentifiedPeak
            /// Peak position of the detection (minutes).
            Center   : float
            /// First, last and highest point of the candidate in the trace set.
            First    : int
            Last     : int
            Apex     : int
            Features : float []
        }

    /// What the candidate features compare the target XIC with.
    type IonContext =
        {
            PredictedScanTime    : float
            /// Candidates are peaks with their position within +-SelectionWindow of the predicted scan time (minutes).
            SelectionWindow      : float
            /// Local scale of the scan time prediction (minutes).
            ScanTimeScale        : float
            /// Source trace sorted by scan time.
            SourceRt             : float []
            SourceIntensity      : float []
            SourceApexRt         : float
            SourceFwhm           : float
            LockMz               : float
            /// nan without ion mobility.
            PredictedIonMobility : float
            IonMobilityScale     : float
        }

    /// Full width at half maximum of the peak around its highest point among first..last: the contiguous points at or above
    /// half of it, one point spacing added.
    let fwhm (rt: float []) (y: float []) (first: int) (last: int) (spacing: float) =
        if last < first then nan
        else
            let apex = [|first .. last|] |> Array.maxBy (fun i -> y.[i])
            let half = 0.5 * y.[apex]
            if half <= 0. then nan
            else
                let mutable lo = apex
                while lo > first && y.[lo - 1] >= half do lo <- lo - 1
                let mutable hi = apex
                while hi < last && y.[hi + 1] >= half do hi <- hi + 1
                rt.[hi] - rt.[lo] + (if isFinite spacing then spacing else 0.)

    /// The source trace is used without its points with non-finite scan time or intensity.
    let createIonContext predictedScanTime selectionWindow scanTimeScale (sourceRt: float []) (sourceIntensity: float []) lockMz predictedIonMobility ionMobilityScale =
        let order =
            if sourceIntensity.Length <> sourceRt.Length then [||]
            else Array.init sourceRt.Length id |> Array.filter (fun i -> isFinite sourceRt.[i] && isFinite sourceIntensity.[i]) |> Array.sortBy (fun i -> sourceRt.[i])
        let rt = order |> Array.map (fun i -> sourceRt.[i])
        let intensity = order |> Array.map (fun i -> sourceIntensity.[i])
        let usable = rt.Length >= 4 && (intensity |> Array.exists (fun v -> v > 0.))
        {
            PredictedScanTime    = predictedScanTime
            SelectionWindow      = selectionWindow
            ScanTimeScale        = scanTimeScale
            SourceRt             = if usable then rt else [||]
            SourceIntensity      = if usable then intensity else [||]
            SourceApexRt         = if usable then rt.[Array.maxBy (fun i -> intensity.[i]) [|0 .. intensity.Length - 1|]] else nan
            SourceFwhm           = if usable then fwhm rt intensity 0 (rt.Length - 1) (medianSpacing rt) else nan
            LockMz               = lockMz
            PredictedIonMobility = predictedIonMobility
            IonMobilityScale     = ionMobilityScale
        }

    /// Intensity weighted mean ion mobility of the monoisotopic points with the given indices.
    let weightedIonMobility (traces: TraceSet) (indices: seq<int>) =
        let mutable sum = 0.
        let mutable sumIm = 0.
        for i in indices do
            let w = traces.Intensity.[0].[i]
            let im = traces.IonMobility.[0].[i]
            if w > 0. && isFinite im then
                sum <- sum + w
                sumIm <- sumIm + w * im
        if sum > 0. then sumIm / sum else nan

    /// The indices of the candidate points at or above half of its apex.
    let topHalf (traces: TraceSet) (c: Candidate) =
        let half = 0.5 * traces.Corrected.[c.Apex]
        [| for i in c.First .. c.Last do if traces.Corrected.[i] >= half then yield i |]

    /// Candidates with their features from the detected peaks of the monoisotopic trace.
    let createCandidates (ctx: IonContext) (traces: TraceSet) (peaks: FSharp.Stats.Signal.PeakDetection.IdentifiedPeak []) =
        let rt = traces.Rt
        let y = traces.Corrected
        let n = rt.Length
        let spacing = medianSpacing rt
        let located =
            peaks
            |> Array.choose (fun peak ->
                let first = SortedSearch.lowerBound rt peak.LeftEnd.XVal |> min (n - 1)
                let last = (SortedSearch.upperBound rt peak.RightEnd.XVal - 1) |> max first
                let apex = [|first .. last|] |> Array.maxBy (fun i -> y.[i])
                if y.[apex] > 0. && isFinite peak.Apex.XVal && abs (peak.Apex.XVal - ctx.PredictedScanTime) <= ctx.SelectionWindow then Some (peak, first, last, apex) else None
                )
        if located.Length = 0 then
            [||]
        else
        let maxHeight = located |> Array.map (fun (_, _, _, apex) -> y.[apex]) |> Array.max
        let hasIsotope = traces.Intensity.Length > 1
        located
        |> Array.map (fun (peak, first, last, apex) ->
            let center = peak.Apex.XVal
            let features = Array.create featureNames.Length nan
            let d = abs (center - ctx.PredictedScanTime) / ctx.ScanTimeScale
            features.[RtZ]  <- d
            features.[RtZ2] <- d * d
            // source trace shifted with its apex onto the candidate, over the candidate +- 3 points
            if ctx.SourceRt.Length > 0 then
                let a = max 0 (first - 3)
                let b = min (n - 1) (last + 3)
                let shifted = [| for i in a .. b do yield interpolateOrZero ctx.SourceRt ctx.SourceIntensity (rt.[i] - center + ctx.SourceApexRt) |]
                features.[ShapeCorr] <- pearson y.[a .. b] shifted
                let width = fwhm rt y first last spacing
                if width > 0. && ctx.SourceFwhm > 0. then
                    features.[AbsLog2FwhmRatio] <- abs (Math.Log2 (width / ctx.SourceFwhm))
            features.[LogRelHeight] <- Math.Log10 (max (y.[apex] / maxHeight) 1e-3)
            features.[Log2NCandidates] <- Math.Log2 (float located.Length)
            if hasIsotope then
                let m0 = y.[first .. last]
                let m1 = traces.Intensity.[1].[first .. last]
                features.[IsoCorr] <- pearson m0 m1
                features.[IsoRatioLog2] <- Math.Log2 ((Array.sum m1 + 1.) / (Array.sum m0 + 1.))
            let ppm =
                [| for i in max 0 (apex - 2) .. min (n - 1) (apex + 2) do
                    let observed = traces.ObservedMz.[0].[i]
                    if isFinite observed then yield (observed - ctx.LockMz) / ctx.LockMz * 1e6 |]
            features.[AbsPpmMono] <- abs (median ppm)
            let candidate = { Peak = peak; Center = center; First = first; Last = last; Apex = apex; Features = features }
            if isFinite ctx.PredictedIonMobility then
                let im = weightedIonMobility traces (topHalf traces candidate)
                features.[ImZ] <- abs (im - ctx.PredictedIonMobility) / ctx.IonMobilityScale
            candidate
            )

    /// Score of the candidate when no model can be trained: Gaussian prior on the scan time distance, shape correlation with the source and relative height.
    let fixedScore (features: float []) =
        let shape = if isFinite features.[ShapeCorr] then features.[ShapeCorr] else 0.
        -0.5 * min features.[RtZ2] 50. + 4. * shape + 2. * features.[LogRelHeight]

    /// L2-penalised logistic regression on standardised features (nan replaced by the training median), fitted by Newton's method.
    type LogisticModel =
        {
            /// Indices into Candidate.Features.
            Features : int []
            Fill     : float []
            Mean     : float []
            Sd       : float []
            /// Intercept, then one weight per standardised feature.
            Weights  : float []
        }

    let private standardise (model: LogisticModel) (features: float []) =
        model.Features
        |> Array.mapi (fun j f ->
            let v = if isFinite features.[f] then features.[f] else model.Fill.[j]
            (v - model.Mean.[j]) / model.Sd.[j]
            )

    let private sigmoid x = 1. / (1. + exp (-x))

    let probability (model: LogisticModel) (features: float []) =
        let z = standardise model features
        let mutable eta = model.Weights.[0]
        for j in 0 .. z.Length - 1 do
            eta <- eta + model.Weights.[j + 1] * z.[j]
        sigmoid eta

    /// Solves a x = b by Gaussian elimination with partial pivoting (a and b are overwritten).
    let private solve (a: float [,]) (b: float []) =
        let n = b.Length
        for col in 0 .. n - 1 do
            let pivot = [|col .. n - 1|] |> Array.maxBy (fun r -> abs a.[r, col])
            if pivot <> col then
                for c in 0 .. n - 1 do
                    let t = a.[col, c]
                    a.[col, c] <- a.[pivot, c]
                    a.[pivot, c] <- t
                let t = b.[col]
                b.[col] <- b.[pivot]
                b.[pivot] <- t
            for r in col + 1 .. n - 1 do
                let f = a.[r, col] / a.[col, col]
                for c in col .. n - 1 do
                    a.[r, c] <- a.[r, c] - f * a.[col, c]
                b.[r] <- b.[r] - f * b.[col]
        let x = Array.zeroCreate n
        for r in n - 1 .. -1 .. 0 do
            let mutable s = b.[r]
            for c in r + 1 .. n - 1 do
                s <- s - a.[r, c] * x.[c]
            x.[r] <- s / a.[r, r]
        x

    let fitLogistic (lambda: float) (features: int []) (x: float [] []) (y: float []) =
        let p = features.Length
        let column j = x |> Array.map (fun row -> row.[features.[j]]) |> Array.filter isFinite
        let fill = Array.init p (fun j -> let c = column j in if c.Length = 0 then 0. else median c)
        let filled = x |> Array.map (fun row -> Array.init p (fun j -> let v = row.[features.[j]] in if isFinite v then v else fill.[j]))
        let mean = Array.init p (fun j -> filled |> Array.averageBy (fun row -> row.[j]))
        let sd =
            Array.init p (fun j ->
                let s = sqrt (filled |> Array.averageBy (fun row -> (row.[j] - mean.[j]) ** 2.))
                if s > 0. then s else 1.)
        let z = filled |> Array.map (fun row -> Array.init (p + 1) (fun j -> if j = 0 then 1. else (row.[j - 1] - mean.[j - 1]) / sd.[j - 1]))
        let w = Array.zeroCreate (p + 1)
        let mutable iteration = 0
        let mutable converged = false
        while not converged && iteration < 100 do
            let gradient = Array.zeroCreate (p + 1)
            let hessian = Array2D.zeroCreate (p + 1) (p + 1)
            for i in 0 .. z.Length - 1 do
                let row = z.[i]
                let mutable eta = 0.
                for j in 0 .. p do
                    eta <- eta + w.[j] * row.[j]
                let mu = sigmoid eta
                let r = mu - y.[i]
                let s = mu * (1. - mu)
                for a in 0 .. p do
                    gradient.[a] <- gradient.[a] + r * row.[a]
                    for b in a .. p do
                        hessian.[a, b] <- hessian.[a, b] + s * row.[a] * row.[b]
            for a in 0 .. p do
                for b in 0 .. a - 1 do
                    hessian.[a, b] <- hessian.[b, a]
            // no penalty on the intercept
            for j in 1 .. p do
                gradient.[j] <- gradient.[j] + lambda * w.[j]
                hessian.[j, j] <- hessian.[j, j] + lambda
            let step = solve hessian gradient
            for j in 0 .. p do
                w.[j] <- w.[j] - step.[j]
            converged <- step |> Array.forall (fun s -> abs s < 1e-8)
            iteration <- iteration + 1
        {
            Features = features
            Fill     = fill
            Mean     = mean
            Sd       = sd
            Weights  = w
        }


module AlignmentBasedQuantification =

    open Numerics
    open XicSweep
    open PeakSelection

    /// Lower bound of the local scan time scale in the candidate distance RtZ (minutes).
    let minScanTimeScale = 0.03
    /// IonMobilityWindow Estimate: +-max(minIonMobilityWindow, ionMobilityWindowSdFactor * PredictedIonMobilitySd). The minimum is the
    /// default window of PSMBasedQuantificationTIMs; with +-0.02 the quantities of the ions the run quantified itself were 10 % lower.
    let minIonMobilityWindow = 0.05
    let ionMobilityWindowSdFactor = 4.
    /// Ions without predicted but with a source ion mobility get a window of +-0.05 around it.
    let sourceIonMobilityWindow = 0.05
    /// Lower bound of the ion mobility scale in the candidate distance ImZ.
    let minIonMobilityScale = 0.001
    /// The m/z calibration is measured on every k-th ion the run quantified itself, k = max 1 (ions / maxCalibrationIons) (spread over
    /// the scan time, so up to 2 * maxCalibrationIons - 1 ions),
    /// with the peak closest to the theoretical m/z plus the initial correction within +-calibrationMzWindow, and needs at least
    /// minCalibrationIons points in the fit. With fewer, the initial correction centres the XICs and the window is
    /// MzCalibration.fallbackMzWindow.
    let maxCalibrationIons = 2000
    let minCalibrationIons = 50
    let calibrationMzWindow = 0.05
    /// The selection model is trained when at least this many clean overlap ions have a candidate; else the fixed score is used.
    let minTrainingIons = 200
    let logisticLambda = 1.
    /// Number of cross-fitting folds of the selection model.
    let selectionFolds = 5
    /// Decoy transfers: the predicted scan time is shifted by decoyShiftMin to decoyShiftMax minutes with random sign (the other sign if
    /// the first leaves the decoyGradientQuantile to 1 - decoyGradientQuantile quantiles of the run's fitted scan times). A decoy is left
    /// out if the scan time of the ion (its chosen candidate, else its prediction) is within XicExtraction.ScanTimeWindow + decoyMargin
    /// of the decoy, so that the ion's own peak is not in the decoy XIC.
    let decoyShiftMin = 4.
    let decoyShiftMax = 6.
    let decoyGradientQuantile = 0.05
    let decoyMargin = 1.5

    /// Assigns every peptide (PepSequenceID, so all modified forms of a sequence) to one of nFolds folds by a seeded permutation
    /// of the sorted ids, as QuantBasedAlignment.
    let assignFolds (seed:int) (nFolds:int) (peptides: seq<int>) =
        let sorted =
            peptides
            |> Seq.distinct
            |> Array.ofSeq
            |> Array.sort
        let rnd = Random(seed)
        for i in sorted.Length - 1 .. -1 .. 1 do
            let j = rnd.Next(i + 1)
            let tmp = sorted.[i]
            sorted.[i] <- sorted.[j]
            sorted.[j] <- tmp
        let foldOf = Dictionary<int,int>()
        sorted |> Array.iteri (fun i peptide -> foldOf.[peptide] <- i % nFolds)
        foldOf
    /// A predicted isotope peak, the monoisotopic peak first.
    type IsotopePeak =
        {
            Mz           : float
            RelFrequency : float
        }

    /// The predicted isotope peaks within the isotope range, summed per isotope; monoisotopic peak first, M+1 second.
    let isotopePeaks (getPattern: SearchDB.LookUpResult<AminoAcids.AminoAcid> -> int -> (float * float) []) (charge: int) (peptide: SearchDB.LookUpResult<AminoAcids.AminoAcid>) =
        let monoMz = Mass.toMZ peptide.Mass (float charge)
        let spacing = 1.0033548378 / float charge
        let byIsotope =
            getPattern peptide charge
            |> Array.groupBy (fun (mz, _) -> int (Math.Round ((mz - monoMz) / spacing)))
            |> Array.map (fun (k, peaks) ->
                let frequency = peaks |> Array.sumBy snd
                let mz = if k = 0 then monoMz else (peaks |> Array.sumBy (fun (mz, f) -> mz * f)) / frequency
                k, { Mz = mz; RelFrequency = frequency }
                )
            |> dict
        let get k = match byIsotope.TryGetValue k with | true, p -> p | _ -> { Mz = monoMz + float k * spacing; RelFrequency = 0. }
        let others =
            byIsotope.Keys
            |> Seq.filter (fun k -> k <> 0 && k <> 1 && float k * spacing > isotopeRangeLow && float k * spacing <= isotopeRangeHigh)
            |> Seq.sort
            |> Seq.map get
            |> Array.ofSeq
        Array.append [|get 0; get 1|] others

    /// Measured against predicted isotope pattern at the XIC point closest to the scan time, as AlignmentBasedQuantification:
    /// the baseline offset of the monoisotopic trace at that point is subtracted from every isotope peak.
    let compareIsotopicCluster (traces: TraceSet) (isotopes: IsotopePeak []) (scanTime: float) =
        if traces.Rt.Length = 0 || not (isFinite scanTime) then
            { PeakComparisons = [||]; KLDiv_UnCorrected = nan; KLDiv_Corrected = nan }
        else
        let c = Array.init traces.Rt.Length id |> Array.minBy (fun i -> abs (traces.Rt.[i] - scanTime))
        let offset = traces.Intensity.[0].[c] - traces.Corrected.[c]
        let measured =
            isotopes
            |> Array.mapi (fun k isotope -> k, isotope)
            |> Array.choose (fun (k, isotope) ->
                let intensity = traces.Intensity.[k].[c]
                let mz = traces.ObservedMz.[k].[c]
                if intensity > 0. && isFinite mz && isotope.RelFrequency > 0. && intensity - offset > 0. then
                    Some { Mz = mz; MeasuredIntensity = intensity; MeasuredIntensityCorrected = intensity - offset; PredictedRelFrequency = isotope.RelFrequency }
                else
                    None
                )
        if measured.Length = 0 then
            { PeakComparisons = [||]; KLDiv_UnCorrected = nan; KLDiv_Corrected = nan }
        else normalizeClusterComparison measured

    /// The precursor m/z written by the PSM based quantification with the theoretical m/z.
    let precursorObservations (psmQuant: QuantificationResult []) =
        psmQuant
        |> Array.map (fun qp ->
            let observation : MzCalibration.Observation =
                {
                    ModSequenceID = qp.ModSequenceID
                    Charge        = qp.Charge
                    GlobalMod     = qp.GlobalMod
                    ScanTime      = getTargetScanTime qp
                    TheoreticalMz = Mass.toMZ qp.TheoMass (float qp.Charge)
                    MeasuredMz    = qp.PrecursorMZ
                }
            observation
            )

    /// An ion of the .align file with its extraction windows.
    type PreparedIon =
        {
            Row               : int
            Align             : AlignmentResult
            Target            : SearchDB.LookUpResult<AminoAcids.AminoAcid>
            Partner           : SearchDB.LookUpResult<AminoAcids.AminoAcid> option
            TargetIsotopes    : IsotopePeak []
            PartnerIsotopes   : IsotopePeak []
            /// Corrected monoisotopic m/z of the target.
            LockMz            : float
            PartnerLockMz     : float
            /// The labeled partner of an overlap ion with a fitted PSM based peak is not extracted nor quantified.
            PartnerSkipped    : bool
            ScanTimeWindow    : float
            IonMobilityWindow : float
            /// The predicted, else the source, ion mobility of the ion; nan without ion mobility window.
            ExpectedIonMobility : float
            Request           : XicRequest
        }

    /// The ion after XIC extraction and peak detection.
    type AnalysedIon =
        {
            Ion        : PreparedIon
            Context    : IonContext
            Target     : TraceSet
            Partner    : TraceSet option
            Candidates : Candidate []
            /// The ion mobility the traces were summed around: the centre from the signal for a fixed window, else the expected ion
            /// mobility of the ion; nan without ion mobility window.
            IonMobilityCentre : float
        }

    /// Selection features of an ion, written to <run>.alignfeatures.
    type AlignFeatures =
        {
            StringSequence       : string
            GlobalMod            : int
            Charge               : int
            PepSequenceID        : int
            ModSequenceID        : int
            /// Model probability of the chosen candidate (nan for the fixed score).
            CandidateProbability : float
            /// Probability (or fixed score) of the chosen minus the second best candidate; for one candidate its probability (nan for the fixed score).
            CandidateMargin      : float
            NCandidates          : int
            AbsPpmMono           : float
            IsoCorr              : float
            IsoRatioLog2         : float
            ShapeCorr            : float
            RtZ                  : float
            ImZ                  : float
            LogRelHeight         : float
            AbsLog2FwhmRatio     : float
            CandidateScore       : float
            SelectionModel       : string
            CandidateScanTime    : float
            ScanTimeWindow       : float
            MzWindow             : float
            IonMobilityWindow    : float
            Status               : string
        }

    /// One side (light or heavy) of a quantification row.
    type QuantSide =
        {
            QuantMz      : float
            Quant        : float
            MeasuredApex : float
            Seo          : float
            Params       : float []
            Difference   : float
            Cluster      : ClusterComparison
            Traces       : TraceSet option
        }

    let emptySide quantMz (traces: TraceSet option) =
        {
            QuantMz      = quantMz
            Quant        = nan
            MeasuredApex = nan
            Seo          = nan
            Params       = [||]
            Difference   = nan
            Cluster      = { PeakComparisons = [||]; KLDiv_UnCorrected = nan; KLDiv_Corrected = nan }
            Traces       = traces
        }

    let createRow (ion: PreparedIon) (measuredMz: float) (target: QuantSide) (partner: QuantSide) (correlation: float) (ionMobility: float) (ionMobilityCentre: float) : QuantificationResult =
        let align = ion.Align
        let light, heavy = if align.GlobalMod = 0 then target, partner else partner, target
        let measuredMass = Mass.ofMZ (if isFinite measuredMz then measuredMz else ion.LockMz) (float align.Charge)
        let traceOf (side: QuantSide) f = match side.Traces with Some t -> f t | None -> [||]
        {
            StringSequence                              = align.StringSequence
            GlobalMod                                   = align.GlobalMod
            Charge                                      = align.Charge
            PepSequenceID                               = align.PepSequenceID
            ModSequenceID                               = align.ModSequenceID
            PrecursorMZ                                 = ion.LockMz
            MeasuredMass                                = measuredMass
            TheoMass                                    = ion.Target.Mass
            AbsDeltaMass                                = abs (measuredMass - ion.Target.Mass)
            MeanPercolatorScore                         = nan
            QValue                                      = nan
            PEPValue                                    = nan
            ProteinNames                                = align.ProteinNames
            QuantMz_Light                               = light.QuantMz
            Quant_Light                                 = light.Quant
            MeasuredApex_Light                          = light.MeasuredApex
            Seo_Light                                   = light.Seo
            Params_Light                                = light.Params
            Difference_SearchRT_FittedRT_Light          = light.Difference
            KLDiv_Observed_Theoretical_Light            = light.Cluster.KLDiv_UnCorrected
            KLDiv_CorrectedObserved_Theoretical_Light   = light.Cluster.KLDiv_Corrected
            QuantMz_Heavy                               = heavy.QuantMz
            Quant_Heavy                                 = heavy.Quant
            MeasuredApex_Heavy                          = heavy.MeasuredApex
            Seo_Heavy                                   = heavy.Seo
            Params_Heavy                                = heavy.Params
            Difference_SearchRT_FittedRT_Heavy          = heavy.Difference
            KLDiv_Observed_Theoretical_Heavy            = heavy.Cluster.KLDiv_UnCorrected
            KLDiv_CorrectedObserved_Theoretical_Heavy   = heavy.Cluster.KLDiv_Corrected
            Correlation_Light_Heavy                     = correlation
            QuantificationSource                        = QuantificationSource.Alignment
            IsotopicPatternMz_Light                     = light.Cluster.PeakComparisons |> Array.map (fun x -> x.Mz)
            IsotopicPatternIntensity_Observed_Light     = light.Cluster.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
            IsotopicPatternIntensity_Corrected_Light    = light.Cluster.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
            RtTrace_Light                               = traceOf light (fun t -> t.Rt)
            IntensityTrace_Observed_Light               = traceOf light (fun t -> t.Intensity.[0])
            IntensityTrace_Corrected_Light              = traceOf light (fun t -> t.Corrected)
            IsotopicPatternMz_Heavy                     = heavy.Cluster.PeakComparisons |> Array.map (fun x -> x.Mz)
            IsotopicPatternIntensity_Observed_Heavy     = heavy.Cluster.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
            IsotopicPatternIntensity_Corrected_Heavy    = heavy.Cluster.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
            RtTrace_Heavy                               = traceOf heavy (fun t -> t.Rt)
            IntensityTrace_Observed_Heavy               = traceOf heavy (fun t -> t.Intensity.[0])
            IntensityTrace_Corrected_Heavy              = traceOf heavy (fun t -> t.Corrected)
            AlignmentScore                              = nan
            AlignmentQValue                             = nan
            IonMobility                                 = ionMobility
            IonMobilityCentre                           = ionMobilityCentre
        }

    /// Counts per reason, for the run summary.
    type private ReasonCounter() =
        let counts = SortedDictionary<string,int>()
        member this.Add (reason: string) =
            counts.[reason] <- (match counts.TryGetValue reason with | true, c -> c + 1 | _ -> 1)
        member this.Summary =
            if counts.Count = 0 then "none"
            else counts |> Seq.map (fun kv -> sprintf "%s: %i" kv.Key kv.Value) |> String.concat ", "

    let private ionName (align: AlignmentResult) = sprintf "%s, GlobalMod %i, charge %i" align.StringSequence align.GlobalMod align.Charge

    /// What the run needs from the PSM based quantification. The rows keep their traces only in PsmRowsHeld, when they could not be written early.
    type private PsmBasedRun =
        { PsmCount: int
          PsmFitted: int
          PsmWrittenEarly: bool
          PsmWritten: int
          PsmRowsHeld: QuantificationResult []
          OwnQuant: Dictionary<(string*int*int),QuantificationResult>
          PsmPeakIons: HashSet<(string*int*int)>
          PsmIons: HashSet<(string*int*int)>
          PsmMzObservations: MzCalibration.Observation [] }

    /// Reads the PSM based rows once and collects what the run needs from them. If every row is its own group in the output,
    /// the rows with a fitted peak are written to outTmpPath right away and the alignment rows are appended later, so the
    /// traces of all PSM based rows are not held while the run is processed.
    let private readPsmBasedQuantification (psmbasedQuant: string) (outTmpPath: string) =
        let psmQuant =
            Csv.CsvReader<QuantificationResult>().ReadFile(psmbasedQuant,'\t',false,1)
            |> Array.ofSeq
        // own PSM quantification of the ions, the label of the selection model; without the traces
        let withoutTraces (qp: QuantificationResult) =
            { qp with
                RtTrace_Light = [||]; IntensityTrace_Observed_Light = [||]; IntensityTrace_Corrected_Light = [||]
                RtTrace_Heavy = [||]; IntensityTrace_Observed_Heavy = [||]; IntensityTrace_Corrected_Heavy = [||] }
        let ownQuant = Dictionary<(string*int*int),QuantificationResult>(HashIdentity.Structural)
        for qp in psmQuant do
            let scanTime = getTargetPeakApex qp
            let sd = getTargetPeakWidth qp
            if isFinite scanTime && isFinite sd && sd > 0. then
                ownQuant.[(qp.StringSequence, qp.GlobalMod, qp.Charge)] <- withoutTraces qp
        // ions whose PSM based row has a fitted scan time: AlignmentBasedQuantStatistics keeps this row, the alignment row is a training example and the truth
        let psmPeakIons = HashSet<(string*int*int)>(HashIdentity.Structural)
        for qp in psmQuant do
            let key = (qp.StringSequence, qp.GlobalMod, qp.Charge)
            if isFinite (getTargetPeakApex qp) then psmPeakIons.Add key |> ignore else psmPeakIons.Remove key |> ignore
        let psmIons = HashSet<(string*int*int)>(psmQuant |> Seq.map (fun qp -> qp.StringSequence, qp.GlobalMod, qp.Charge), HashIdentity.Structural)
        let psmMzObservations = precursorObservations psmQuant
        let psmRows = psmQuant |> Array.filter (fun qp -> (tryTargetGetScanTime qp).IsSome)
        let keyOf (qp: QuantificationResult) = qp.StringSequence, qp.ModSequenceID, qp.Charge, qp.GlobalMod, qp.QuantificationSource
        let ownGroups =
            psmRows |> Array.forall (fun qp -> qp.QuantificationSource = QuantificationSource.PSM)
            && (psmRows |> Array.distinctBy keyOf).Length = psmRows.Length
        if ownGroups then
            File.WriteAllLines(outTmpPath, psmRows |> SeqIO'.csv "\t" true false)
        { PsmCount = psmQuant.Length
          PsmFitted = psmRows.Length
          PsmWrittenEarly = ownGroups
          PsmWritten = (if ownGroups then psmRows.Length else 0)
          PsmRowsHeld = (if ownGroups then [||] else psmRows)
          OwnQuant = ownQuant
          PsmPeakIons = psmPeakIons
          PsmIons = psmIons
          PsmMzObservations = psmMzObservations }

    /// Quantifies the ions of one .align file in one run and writes <run>.quant and <run>.alignfeatures.
    let quantifyPeptides diagCharts (processParams:Domain.AlignmentBasedQuantificationParams) (outputDir:string) (d:string) (instrumentOutput:string) (psmbasedQuant:string) (alignmentResults:string) =
        let runName = Path.GetFileNameWithoutExtension instrumentOutput
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension alignmentResults)
        let stopwatch = Diagnostics.Stopwatch.StartNew()
        let elapsed () = stopwatch.Elapsed.TotalSeconds
        logger.Info (sprintf "%s: quantifying the ions of %s in %s" runName alignmentResults instrumentOutput)
        logger.Trace (sprintf "Parameters: %A" processParams)
        let outFilePath = Path.Combine [|outputDir; runName + ".quant"|]
        let featureFilePath = Path.Combine [|outputDir; runName + ".alignfeatures"|]
        let decoyFilePath = Path.Combine [|outputDir; runName + ".aligndecoy"|]

        // inputs
        let alignments =
            Csv.CsvReader<AlignmentResult>().ReadFile(alignmentResults,'\t',false,1)
            |> Array.ofSeq
        let outTmpPath = outFilePath + ".tmp"
        if File.Exists outTmpPath then File.Delete outTmpPath
        // removes the file of the early write if the run fails; after the move there is nothing left to remove
        use _tmpGuard = { new IDisposable with member _.Dispose() = try (if File.Exists outTmpPath then File.Delete outTmpPath) with _ -> () }
        let psm = readPsmBasedQuantification psmbasedQuant outTmpPath
        logger.Info (sprintf "%s: %i .align rows, %i PSM based rows (%i ions with a fitted peak)" runName alignments.Length psm.PsmCount psm.OwnQuant.Count)

        // the peptide data base is read from disk: only the ions of the .align file are looked up, and a copy in memory
        // would cost the size of the data base per run (5.5 GB for the TIMs data base)
        use peptideDB = SearchDB.getDBConnection d
        let dBParams = SearchDB.getSDBParamsByCn peptideDB
        let peptideLookUp = SearchDB.getThreadSafePeptideLookUpFromFileBySequenceAndGMod peptideDB dBParams
        use isotopePatterns = new IsotopePatterns.IsotopePatternLookUp(logger, peptideDB)

        use inReader = Core.MzIO.Reader.getReader instrumentOutput :?> MzIO.MzSQL.MzSQL
        inReader.Connection.Open()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader
        let ms1 = SpectrumIndex.readMs1Index inReader inRunID (fun message -> logger.Warn message) ignore
        let readSpectrum = Xic.readSortedPeaks inReader
        let spectraHaveIonMobility =
            ms1
            |> Seq.skip (ms1.Length / 2)
            |> Seq.truncate 20
            |> Seq.map (snd >> readSpectrum)
            |> Seq.tryFind (fun s -> s.Mz.Length > 0)
            |> Option.map (fun s -> s.IonMobility.Length > 0)
            |> Option.defaultValue false
        let alignHasIonMobility = alignments |> Array.exists (fun a -> isFinite a.PredictedIonMobility || isFinite a.IonMobility)
        let performLabeled = processParams.PerformLabeledQuantification <> Domain.Labeling.Unlabeled
        let useIonMobility = spectraHaveIonMobility && alignHasIonMobility
        logger.Info (
            sprintf "%s: %i MS1 spectra, scan time %.2f - %.2f min; ion mobility in the spectra: %b, in the .align file: %b, used: %b (%.1f s)"
                runName ms1.Length (if ms1.Length > 0 then fst ms1.[0] else nan) (if ms1.Length > 0 then fst ms1.[ms1.Length - 1] else nan)
                spectraHaveIonMobility alignHasIonMobility useIonMobility (elapsed ())
            )

        // m/z calibration measured in the MS1 spectra: offset of the peak (with ion mobility the m/z group) closest to the theoretical m/z (within +-calibrationMzWindow)
        // at the fitted scan time of the ions the run quantified itself. The precursor m/z of the .quant rows only places the lock m/z.
        let initialCalibration = MzCalibration.fit psm.PsmMzObservations
        let calibrationIons =
            let ions = psm.OwnQuant.Values |> Array.ofSeq |> Array.sortBy getTargetPeakApex
            let step = max 1 (ions.Length / maxCalibrationIons)
            ions |> Array.indexed |> Array.filter (fun (i, _) -> i % step = 0) |> Array.map snd
        let calibration =
            let ms1ScanTimes = ms1 |> Array.map fst
            // the MS1 spectrum closest to the scan time
            let closestScanTime (scanTime: float) = ms1ScanTimes.[SortedSearch.closestIndex ms1ScanTimes scanTime]
            let requests =
                calibrationIons
                |> Array.map (fun qp ->
                    let scanTime = getTargetPeakApex qp
                    let theoMz = Mass.toMZ qp.TheoMass (float qp.Charge)
                    let imCenter = if useIonMobility then qp.IonMobility else nan
                    {
                        RtLow    = if ms1ScanTimes.Length > 0 then closestScanTime scanTime else nan
                        RtHigh   = if ms1ScanTimes.Length > 0 then closestScanTime scanTime else nan
                        LockMz   = [|theoMz + initialCalibration.CorrectionAt scanTime theoMz|]
                        MzWindow = calibrationMzWindow
                        ImLow    = imCenter - sourceIonMobilityWindow
                        ImHigh   = imCenter + sourceIonMobilityWindow
                        Point    = ClosestPeak
                    }
                    )
            let observations = ResizeArray<MzCalibration.Observation>()
            sweep readSpectrum ms1 requests (fun r raw _ ->
                if raw.Rt.Length > 0 && raw.Intensity.[0].[0] > 0. then
                    let qp = calibrationIons.[r]
                    observations.Add
                        {
                            ModSequenceID = qp.ModSequenceID
                            Charge        = qp.Charge
                            GlobalMod     = qp.GlobalMod
                            ScanTime      = getTargetPeakApex qp
                            TheoreticalMz = Mass.toMZ qp.TheoMass (float qp.Charge)
                            MeasuredMz    = raw.ObservedMz.[0].[0]
                        }
                )
            |> ignore
            MzCalibration.fit (observations.ToArray())
        let measuredUsable = calibration.PointsFence2 >= minCalibrationIons
        // the calibration that centres the XICs, and its window
        let centre, halfWindow =
            if measuredUsable then
                calibration, MzCalibration.residualHalfWindow calibration
            else
                logger.Warn (sprintf "%s: the m/z calibration has %i points after the fences (at least %i are needed); the precursor m/z of the .quant file centres the XICs and the m/z window is +-%.3f Da." runName calibration.PointsFence2 minCalibrationIons MzCalibration.fallbackMzWindow)
                initialCalibration, MzCalibration.fallbackMzWindow
        let mzWindow =
            match processParams.XicExtraction.MzWindow_Da with
            | Domain.Window.Fixed v -> v
            | Domain.Window.Estimate -> halfWindow
        logger.Info (sprintf "%s: XIC %s" runName (MzCalibration.describeWindow processParams.XicExtraction.MzWindow_Da mzWindow (not measuredUsable)))
        let within =
            if centre.Points.Length = 0 then nan
            else
                let inside = centre.Points |> Array.filter (fun p -> abs (p.OffsetPpm * p.Mz * 1e-6 - centre.CorrectionAt p.Time p.Mz) <= mzWindow)
                100. * float inside.Length / float centre.Points.Length
        logger.Info (
            sprintf "%s: m/z calibration measured in the MS1 spectra: %s; precursor m/z of the .quant file: %s; the m/z window contains %.1f %% of the points of the correction in use (%.1f s)"
                runName (MzCalibration.describe calibration) (MzCalibration.describe initialCalibration) within (elapsed ())
            )
        if diagCharts then
            MzCalibration.calibrationChart centre
            |> Chart.withTitle (sprintf "%s: m/z correction, m/z window +-%.4f Da" runName mzWindow)
            |> Chart.saveHtml (Path.Combine [|outputDir; runName + "_mzCorrection"|])

        let identifyPeaks = initIdentifyPeaks processParams.XicExtraction.XicProcessing
        // rows 0 .. n - 1 are the .align rows, row n + r is the decoy transfer of row r
        let n = alignments.Length
        let decoyAlign : AlignmentResult option [] = Array.create n None
        let alignOf (row: int) = if row < n then alignments.[row] else decoyAlign.[row - n].Value
        let nameOf (row: int) = if row < n then ionName alignments.[row] else ionName (alignOf row) + " (decoy)"
        let drops = ReasonCounter()
        let decoySkips = ReasonCounter()
        let features = Array.create (2 * n) None
        let status = Array.create (2 * n) "quantified"
        let dropIon (row: int) (reason: string) =
            logger.Debug (sprintf "Ion not quantified, %s: %s" reason (nameOf row))
            if row < n then drops.Add reason
            status.[row] <- reason
        // runs the processing of one ion; an exception drops the ion instead of the run
        let guarded (row: int) (step: string) (f: unit -> 'T option) =
            try
                f ()
            with ex ->
                logger.Debug (sprintf "%s failed (%s): %s" step ex.Message (nameOf row))
                dropIon row (sprintf "error in %s" step)
                None

        // the light and heavy row of a peptide and its charges share the data base lookup and the isotope pattern
        let lookUps = Dictionary<(string*int), SearchDB.LookUpResult<AminoAcids.AminoAcid> option>(HashIdentity.Structural)
        let lookUpCalls = ref 0
        let isotopeCalls = ref 0
        let lookUp sequence gMod =
            lookUpCalls.Value <- lookUpCalls.Value + 1
            match lookUps.TryGetValue ((sequence, gMod)) with
            | true, found -> found
            | _ ->
                // a sequence without a row in the data base comes back with the ID -1
                let found = Some (peptideLookUp sequence gMod) |> Option.filter (fun peptide -> peptide.ModSequenceID >= 0)
                lookUps.[(sequence, gMod)] <- found
                found
        let isotopes = Dictionary<(string*int*int), IsotopePeak []>(HashIdentity.Structural)
        let isotopePeaksOf (align: AlignmentResult) gMod (peptide: SearchDB.LookUpResult<AminoAcids.AminoAcid>) =
            let key = (align.StringSequence, gMod, align.Charge)
            isotopeCalls.Value <- isotopeCalls.Value + 1
            match isotopes.TryGetValue key with
            | true, peaks -> peaks
            | _ ->
                let peaks = isotopePeaks (fun lookUpResult charge -> isotopePatterns.Get (lookUpResult, charge)) align.Charge peptide
                isotopes.[key] <- peaks
                peaks

        // an ion without a fitted quantification of its own in the run (no PSM based row, or one without a fitted peak) whose labeled partner
        // the run quantified through its own PSM based row with a fitted peak, with a value of the label of the ion, is not transferred:
        // the PSM based row has more evidence than a match through alignment
        let labeledPartnerQuantified (align: AlignmentResult) =
            performLabeled
            && not (psm.OwnQuant.ContainsKey ((align.StringSequence, align.GlobalMod, align.Charge)))
            && (match psm.OwnQuant.TryGetValue ((align.StringSequence, (if align.GlobalMod = 0 then 1 else 0), align.Charge)) with
                | true, qp ->
                    let value = if align.GlobalMod = 0 then qp.Quant_Light else qp.Quant_Heavy
                    isFinite value && value > 0.
                | _ -> false)

        // windows per ion
        let prepare (row: int) (align: AlignmentResult) =
            guarded row "window preparation" <| fun () ->
                if row < n && labeledPartnerQuantified align then
                    dropIon row LabeledPartnerQuantifiedStatus
                    None
                elif not performLabeled && align.GlobalMod <> 0 then
                    dropIon row "label-free quantification of an ion with GlobalMod <> 0"
                    None
                elif not (isFinite align.PredictedScanTime) then
                    dropIon row "no predicted scan time"
                    None
                else
                match lookUp align.StringSequence align.GlobalMod with
                | None ->
                    dropIon row "peptide not in the database"
                    None
                | Some target ->
                // the alignment row of an overlap ion with a fitted PSM based peak is a training example and the truth of
                // AlignmentBasedQuantStatistics only, which reads the quantification of the ion itself
                let partnerSkipped = performLabeled && psm.PsmPeakIons.Contains ((align.StringSequence, align.GlobalMod, align.Charge))
                let partner =
                    if performLabeled && not partnerSkipped then
                        let p = lookUp align.StringSequence (if align.GlobalMod = 0 then 1 else 0)
                        if p.IsNone then logger.Debug (sprintf "Labeled partner not in the database: %s" (ionName align))
                        p
                    else None
                let correctionAt mz = centre.CorrectionAt align.PredictedScanTime mz
                let targetIsotopes = isotopePeaksOf align align.GlobalMod target
                let partnerIsotopes = partner |> Option.map (isotopePeaksOf align (if align.GlobalMod = 0 then 1 else 0)) |> Option.defaultValue [||]
                let scanTimeWindow =
                    let sd = align.PredictedScanTimeSd
                    if isFinite sd && sd > 0. then
                        min processParams.XicExtraction.ScanTimeWindow (max processParams.MinScanTimeWindow (processParams.ScanTimeWindowSdFactor * sd))
                    else
                        processParams.XicExtraction.ScanTimeWindow
                let imCenter, imWindow =
                    if not useIonMobility then nan, nan
                    elif isFinite align.PredictedIonMobility then
                        match processParams.XicExtraction.IonMobilityWindow with
                        | Domain.Window.Fixed v -> align.PredictedIonMobility, v
                        | Domain.Window.Estimate ->
                            let sd = align.PredictedIonMobilitySd
                            align.PredictedIonMobility, (if isFinite sd then max minIonMobilityWindow (ionMobilityWindowSdFactor * sd) else sourceIonMobilityWindow)
                    elif isFinite align.IonMobility then
                        match processParams.XicExtraction.IonMobilityWindow with
                        | Domain.Window.Fixed v -> align.IonMobility, v
                        | Domain.Window.Estimate -> align.IonMobility, sourceIonMobilityWindow
                    else nan, nan
                // every trace is centred and windowed at its own m/z
                let targetMz = Mass.toMZ target.Mass (float align.Charge)
                let lockMz = targetMz + correctionAt targetMz
                let partnerMz = match partner with Some p -> Mass.toMZ p.Mass (float align.Charge) | None -> nan
                let partnerLockMz = partnerMz + correctionAt partnerMz
                let theoreticalMz = Array.append targetIsotopes partnerIsotopes |> Array.map (fun p -> p.Mz)
                Some {
                    Row               = row
                    Align             = align
                    Target            = target
                    Partner           = partner
                    TargetIsotopes    = targetIsotopes
                    PartnerIsotopes   = partnerIsotopes
                    LockMz            = lockMz
                    PartnerLockMz     = partnerLockMz
                    PartnerSkipped    = partnerSkipped
                    ScanTimeWindow    = scanTimeWindow
                    IonMobilityWindow = imWindow
                    ExpectedIonMobility = imCenter
                    Request           =
                        {
                            RtLow    = align.PredictedScanTime - processParams.XicExtraction.ScanTimeWindow
                            RtHigh   = align.PredictedScanTime + processParams.XicExtraction.ScanTimeWindow
                            LockMz   = theoreticalMz |> Array.map (fun mz -> mz + correctionAt mz)
                            MzWindow = mzWindow
                            ImLow    = imCenter - imWindow
                            ImHigh   = imCenter + imWindow
                            // a fixed window is summed around the ion mobility centre of the ion, which needs the profiles of the traces
                            Point =
                                match processParams.XicExtraction.IonMobilityWindow with
                                | Domain.Window.Fixed _ when isFinite imCenter -> SumWithProfiles
                                | _ -> Sum
                        }
                }
        let prepared = alignments |> Array.mapi prepare |> Array.choose id
        if performLabeled then
            logger.Info (sprintf "%s: labeled partner not quantified for %i of %i ions (overlap ions with a fitted PSM based peak)" runName (prepared |> Array.filter (fun ion -> ion.PartnerSkipped) |> Array.length) prepared.Length)
        let withIonMobility = prepared |> Array.filter (fun ion -> isFinite ion.Request.ImLow) |> Array.length
        logger.Info (
            sprintf "%s: %i ions to extract (XIC +-%.2f min); candidate window median +-%.3f min (%.3f - %.3f); %i with an ion mobility window (median +-%.4f)"
                runName prepared.Length processParams.XicExtraction.ScanTimeWindow (prepared |> Array.map (fun ion -> ion.ScanTimeWindow) |> median)
                (if prepared.Length > 0 then prepared |> Array.map (fun ion -> ion.ScanTimeWindow) |> Array.min else nan)
                (if prepared.Length > 0 then prepared |> Array.map (fun ion -> ion.ScanTimeWindow) |> Array.max else nan)
                withIonMobility (prepared |> Array.choose (fun ion -> if isFinite ion.IonMobilityWindow then Some ion.IonMobilityWindow else None) |> median)
            )

        // decoy transfers: DecoyFraction of the new ions (no PSM based row in the run), chosen with Seed, extracted again with the
        // predicted scan time shifted to where their own peak is not. They take the path of the other ions, are not used to train the
        // selection model and are written to <run>.aligndecoy only.
        let decoys =
            let pool = prepared |> Array.filter (fun ion -> not (psm.PsmIons.Contains ((ion.Align.StringSequence, ion.Align.GlobalMod, ion.Align.Charge))))
            let ownScanTimes = psm.OwnQuant.Values |> Seq.map getTargetPeakApex |> Array.ofSeq |> Array.sort
            if processParams.DecoyFraction <= 0. || pool.Length = 0 then
                [||]
            elif ownScanTimes.Length = 0 then
                logger.Warn (sprintf "%s: no ion with a fitted peak, so the gradient is unknown; no decoy transfers." runName)
                [||]
            else
            let low = quantileSorted ownScanTimes decoyGradientQuantile
            let high = quantileSorted ownScanTimes (1. - decoyGradientQuantile)
            let rnd = Random(processParams.Seed)
            let order = Array.init pool.Length id
            for i in order.Length - 1 .. -1 .. 1 do
                let j = rnd.Next(i + 1)
                let tmp = order.[i]
                order.[i] <- order.[j]
                order.[j] <- tmp
            let requested = int (Math.Round (processParams.DecoyFraction * float pool.Length))
            let rows =
                order.[.. requested - 1]
                |> Array.sort
                |> Array.choose (fun i ->
                    let ion = pool.[i]
                    let shift = decoyShiftMin + (decoyShiftMax - decoyShiftMin) * rnd.NextDouble()
                    let sign = if rnd.NextDouble() < 0.5 then 1. else -1.
                    let predicted = ion.Align.PredictedScanTime
                    match [sign; -sign] |> List.map (fun s -> predicted + s * shift) |> List.tryFind (fun t -> t >= low && t <= high) with
                    | None ->
                        decoySkips.Add "outside the gradient"
                        None
                    | Some t when abs (t - predicted) <= processParams.XicExtraction.ScanTimeWindow + decoyMargin ->
                        decoySkips.Add "own peak inside the decoy XIC"
                        None
                    | Some t ->
                        decoyAlign.[ion.Row] <- Some { ion.Align with PredictedScanTime = t }
                        Some (n + ion.Row)
                    )
            logger.Info (sprintf "%s: decoy transfers of %i of %i new ions requested (gradient %.2f - %.2f min), %i created" runName requested pool.Length low high rows.Length)
            rows |> Array.choose (fun row -> prepare row (alignOf row))

        logger.Info (sprintf "%s: data base look ups %i (%i hits), isotope patterns %i (%i hits)" runName lookUpCalls.Value (lookUpCalls.Value - lookUps.Count) isotopeCalls.Value (isotopeCalls.Value - isotopes.Count))
        logger.Info (sprintf "%s: isotope patterns read from the peptide data base %i, computed %i" runName isotopePatterns.Read isotopePatterns.Computed)

        // XIC extraction and candidates, as soon as an XIC window is passed. Two sweeps: first the ions the run quantified itself
        // (they train the selection model and are scored by the model of the other folds, which needs the peptides of all
        // ions), then the other ions, which are scored by the model of all training ions. If the second sweep cannot change the model
        // (eager, below) they are fitted as soon as their XIC is analysed, so only their rows, not their traces, are held until the end;
        // otherwise their traces are held until the model is trained.
        let toExtract = Array.append prepared decoys
        let analysed : AnalysedIon option [] = Array.create toExtract.Length None
        let analyse (i: int) (raw: RawXic) (ionMobilityCentre: float) =
            let ion = toExtract.[i]
            analysed.[i] <-
                guarded ion.Row "peak detection" <| fun () ->
                    let align = ion.Align
                    let nTarget = ion.TargetIsotopes.Length
                    let targetTraces = createTraceSet processParams.BaseLineCorrection true raw [|0 .. nTarget - 1|]
                    let partnerTraces =
                        if ion.Partner.IsSome then Some (createTraceSet processParams.BaseLineCorrection false raw [|nTarget .. nTarget + ion.PartnerIsotopes.Length - 1|])
                        else None
                    if targetTraces.Rt.Length < 5 then
                        dropIon ion.Row "XIC with fewer than 5 points"
                        None
                    elif targetTraces.Corrected |> Array.forall (fun v -> v <= 0.) then
                        dropIon ion.Row "empty XIC"
                        None
                    else
                        let peaks =
                            try Ok (identifyPeaks targetTraces.Rt targetTraces.Corrected)
                            with ex -> Error ex.Message
                        match peaks with
                        | Error message ->
                            logger.Trace (sprintf "Peak detection failed (%s): %s" message (ionName align))
                            dropIon ion.Row "peak detection failed"
                            None
                        | Ok peaks ->
                            let context =
                                createIonContext
                                    align.PredictedScanTime ion.ScanTimeWindow (max minScanTimeScale align.PredictedScanTimeSd)
                                    align.RtTrace_SourceFile align.IntensityTrace_SourceFile ion.LockMz
                                    (if isFinite ion.Request.ImLow && isFinite align.PredictedIonMobility && isFinite align.PredictedIonMobilitySd then align.PredictedIonMobility else nan)
                                    (max minIonMobilityScale align.PredictedIonMobilitySd)
                            let candidates = createCandidates context targetTraces peaks
                            if candidates.Length = 0 then
                                dropIon ion.Row "no peak detected in the candidate window"
                                None
                            else
                                Some { Ion = ion; Context = { context with SourceRt = [||]; SourceIntensity = [||] }; Target = targetTraces; Partner = partnerTraces; Candidates = candidates; IonMobilityCentre = ionMobilityCentre }
        let ownOfIon (ion: PreparedIon) =
            match psm.OwnQuant.TryGetValue ((ion.Align.StringSequence, ion.Align.GlobalMod, ion.Align.Charge)) with
            | true, qp -> Some qp
            | _ -> None
        let ownOf (a: AnalysedIon) = ownOfIon a.Ion
        let isCorrect (qp: QuantificationResult) (c: Candidate) = abs (c.Center - getTargetPeakApex qp) <= getTargetPeakWidth qp
        let ownIdx = [| for i in 0 .. toExtract.Length - 1 do if (ownOfIon toExtract.[i]).IsSome then yield i |]
        let newIdx = [| for i in 0 .. toExtract.Length - 1 do if (ownOfIon toExtract.[i]).IsNone then yield i |]
        // The XIC of an ion with a fixed ion mobility window is summed around the centre of its signal before anything else reads it.
        // For the log, the ions whose centre comes from the signal and the ions that keep the expected ion mobility are counted;
        // decoys are left out. Returns the XIC and the ion mobility it is summed around: the centre from the signal, else the
        // expected ion mobility.
        let centreTally = Xic.MobilityCentreTally()
        let sumAroundCentre (i: int) (raw: RawXic) (profiles: Xic.MobilityProfile [] option) =
            let ion = toExtract.[i]
            match profiles with
            | None -> raw, ion.ExpectedIonMobility
            | Some profiles ->
                let summed, centre = sumAroundMobilityCentre ion.ExpectedIonMobility ion.IonMobilityWindow ion.Align.PredictedScanTime raw profiles
                if ion.Row < n then centreTally.Add(ion.ExpectedIonMobility, centre)
                summed, defaultArg centre ion.ExpectedIonMobility
        let sweepOver (indices: int []) (callback: int -> RawXic -> float -> unit) =
            sweep readSpectrum ms1 (indices |> Array.map (fun i -> toExtract.[i].Request)) (fun j raw profiles ->
                let summed, centre = sumAroundCentre indices.[j] raw profiles
                callback indices.[j] summed centre)
        let ownSpectraRead = sweepOver ownIdx analyse
        let mutable analysedOwn = ownIdx |> Array.choose (fun i -> analysed.[i])
        let baseFeatures = [| RtZ; RtZ2; ShapeCorr; LogRelHeight; AbsLog2FwhmRatio; Log2NCandidates; IsoCorr; IsoRatioLog2; AbsPpmMono |]
        let withIonMobilityFeature (exists: bool) = if useIonMobility && exists then Array.append baseFeatures [|ImZ|] else baseFeatures
        // the ion mobility feature is used if any ion with candidates has a predicted ion mobility; the ions of the second sweep
        // can only add one if they have the data for it
        let ownHasIm = analysedOwn |> Array.exists (fun a -> isFinite a.Context.PredictedIonMobility)
        let newMayHaveIm =
            newIdx |> Array.exists (fun i ->
                let ion = toExtract.[i]
                ion.Row < n && isFinite ion.Request.ImLow && isFinite ion.Align.PredictedIonMobility && isFinite ion.Align.PredictedIonMobilitySd)
        let mutable training =
            analysedOwn
            |> Array.choose (fun a ->
                match ownOf a with
                | Some qp when abs (getTargetApexDifference qp) / getTargetPeakWidth qp < 2. -> Some (a, qp)
                | _ -> None
                )
        // the new ions are fitted as soon as they are analysed if the second pass cannot change the model; with too few training ions there is none
        let eager = ownHasIm || not newMayHaveIm || training.Length < minTrainingIons
        let fitOn (featureSet: int []) (ions: (AnalysedIon*QuantificationResult) []) =
            let x = ions |> Array.collect (fun (a, _) -> a.Candidates |> Array.map (fun c -> c.Features))
            let y = ions |> Array.collect (fun (a, qp) -> a.Candidates |> Array.map (fun c -> if isCorrect qp c then 1. else 0.))
            let positives = y |> Array.sum
            if positives = 0. || positives = float y.Length then None
            else Some (fitLogistic logisticLambda featureSet x y, y.Length, positives)
        let trainModel (featureSet: int []) =
            if training.Length < minTrainingIons then
                logger.Warn (sprintf "%s: %i clean overlap ions with candidates, at least %i are needed to train the selection model; the fixed score is used." runName training.Length minTrainingIons)
                None
            else
                match fitOn featureSet training with
                | None ->
                    logger.Warn (sprintf "%s: the training candidates are all of one class; the fixed score is used." runName)
                    None
                | Some (m, nCandidates, positives) ->
                    logger.Info (
                        sprintf "%s: selection model trained on %i clean overlap ions (%i candidates, %.1f %% correct, %.1f %% of the ions have a correct candidate); standardised weights: intercept %.3f, %s"
                            runName training.Length nCandidates (100. * positives / float nCandidates)
                            (100. * (training |> Array.filter (fun (a, qp) -> a.Candidates |> Array.exists (isCorrect qp)) |> Array.length |> float) / float training.Length)
                            m.Weights.[0]
                            (featureSet |> Array.mapi (fun j f -> sprintf "%s %.3f" featureNames.[f] m.Weights.[j + 1]) |> String.concat ", ")
                        )
                    Some m
        // set after the first sweep (then all ions can be scored by it) or, if the second sweep can change the features, after it
        let mutable featureSet = withIonMobilityFeature ownHasIm
        let mutable model : LogisticModel option = None
        let mutable foldOf : Dictionary<int,int> = Dictionary<int,int>()
        let mutable foldModels : LogisticModel option [] = [||]
        // the model that scores the candidates of the ion, None for the fixed score
        let modelOf (a: AnalysedIon) =
            match model with
            | None -> None
            | Some full -> if (ownOf a).IsSome then foldModels.[foldOf.[a.Ion.Align.PepSequenceID]] else Some full
        let selectionModelOf (m: LogisticModel option) = if m.IsSome then "Logistic" else "FixedScore"
        let score (a: AnalysedIon) (c: Candidate) =
            match modelOf a with
            | Some m -> probability m c.Features
            | None -> fixedScore c.Features
        let rank (a: AnalysedIon) = a.Candidates |> Array.map (fun c -> c, score a c) |> Array.sortByDescending snd
        let mutable partnerFailures = 0
        // quantification of the chosen candidate
        let fit (a: AnalysedIon) =
            guarded a.Ion.Row "peak fit" <| fun () ->
                let ion = a.Ion
                let align = ion.Align
                let ionModel = modelOf a
                let ranked = rank a
                let chosen, best = ranked.[0]
                let margin =
                    if ranked.Length > 1 then best - snd ranked.[1]
                    elif ionModel.IsSome then best
                    else nan
                let f = chosen.Features
                let featureRow status =
                    {
                        StringSequence       = align.StringSequence
                        GlobalMod            = align.GlobalMod
                        Charge               = align.Charge
                        PepSequenceID        = align.PepSequenceID
                        ModSequenceID        = align.ModSequenceID
                        CandidateProbability = if ionModel.IsSome then best else nan
                        CandidateMargin      = margin
                        NCandidates          = a.Candidates.Length
                        AbsPpmMono           = f.[AbsPpmMono]
                        IsoCorr              = f.[IsoCorr]
                        IsoRatioLog2         = f.[IsoRatioLog2]
                        ShapeCorr            = f.[ShapeCorr]
                        RtZ                  = f.[RtZ]
                        ImZ                  = f.[ImZ]
                        LogRelHeight         = f.[LogRelHeight]
                        AbsLog2FwhmRatio     = f.[AbsLog2FwhmRatio]
                        CandidateScore       = best
                        SelectionModel       = selectionModelOf ionModel
                        CandidateScanTime    = chosen.Center
                        ScanTimeWindow       = ion.ScanTimeWindow
                        MzWindow             = mzWindow
                        IonMobilityWindow    = ion.IonMobilityWindow
                        Status               = status
                    }
                let quantP = try Some (HULQ.quantifyPeak chosen.Peak) with _ -> None
                match quantP with
                | Some q when isValidFit q ->
                    let fittedParam1 = q.EstimatedParams.[1]
                    let fittedScanTime = PeakShape.apex q.EstimatedParams
                    let targetSide =
                        {
                            QuantMz      = ion.LockMz
                            Quant        = q.Area
                            MeasuredApex = q.MeasuredApexIntensity
                            Seo          = q.StandardErrorOfPrediction
                            Params       = q.EstimatedParams
                            Difference   = align.PredictedScanTime - fittedParam1
                            Cluster      = compareIsotopicCluster a.Target ion.TargetIsotopes fittedParam1
                            Traces       = Some a.Target
                        }
                    let topPoints = topHalf a.Target chosen
                    let measuredMz =
                        let w = topPoints |> Array.filter (fun i -> isFinite a.Target.ObservedMz.[0].[i])
                        let sum = w |> Array.sumBy (fun i -> a.Target.Intensity.[0].[i])
                        if sum > 0. then (w |> Array.sumBy (fun i -> a.Target.Intensity.[0].[i] * a.Target.ObservedMz.[0].[i])) / sum else nan
                    // ion mobility of the chosen peak: intensity weighted over the points within 1 fitted width of its apex
                    let ionMobility =
                        if not useIonMobility then nan
                        else
                            let sd = PeakShape.width q.EstimatedParams
                            let around = [| for i in 0 .. a.Target.Rt.Length - 1 do if abs (a.Target.Rt.[i] - fittedScanTime) <= sd then yield i |]
                            let im = weightedIonMobility a.Target around
                            if isFinite im then im else weightedIonMobility a.Target topPoints
                    let partnerSide, correlation =
                        match ion.Partner, a.Partner with
                        | Some _, Some traces ->
                            let failed reason =
                                logger.Debug (sprintf "Labeled partner not quantified, %s: %s" reason (ionName align))
                                if ion.Row < n then partnerFailures <- partnerFailures + 1
                                emptySide ion.PartnerLockMz (Some traces), nan
                            if traces.Rt.Length < 5 || traces.Corrected |> Array.forall (fun v -> v <= 0.) then
                                failed "empty XIC"
                            else
                            let peaks = try identifyPeaks traces.Rt traces.Corrected with _ -> [||]
                            if peaks.Length = 0 then
                                failed "no peak detected"
                            else
                            // as AlignmentBasedQuantification: the partner peak at the fitted scan time of the ion
                            let anchor = if abs (align.PredictedScanTime - fittedScanTime) <= processParams.XicExtraction.ScanTimeWindow then fittedScanTime else align.PredictedScanTime
                            let partnerPeak = HULQ.getPeakBy peaks anchor
                            match (try Some (HULQ.quantifyPeak partnerPeak) with _ -> None) with
                            | Some qp when isValidFit qp ->
                                {
                                    QuantMz      = ion.PartnerLockMz
                                    Quant        = qp.Area
                                    MeasuredApex = qp.MeasuredApexIntensity
                                    Seo          = qp.StandardErrorOfPrediction
                                    Params       = qp.EstimatedParams
                                    Difference   = anchor - qp.EstimatedParams.[1]
                                    Cluster      = compareIsotopicCluster traces ion.PartnerIsotopes fittedParam1
                                    Traces       = Some traces
                                }, calcCorrelation a.Target.Rt q qp
                            | _ -> failed "peak fit failed"
                        | Some _, None -> emptySide ion.PartnerLockMz None, nan
                        | None, _ ->
                            if performLabeled && not ion.PartnerSkipped && ion.Row < n then
                                partnerFailures <- partnerFailures + 1
                            emptySide nan None, nan
                    features.[ion.Row] <- Some (featureRow "quantified")
                    Some (createRow ion measuredMz targetSide partnerSide correlation ionMobility a.IonMobilityCentre)
                | _ ->
                    features.[ion.Row] <- Some (featureRow "peak fit failed")
                    dropIon ion.Row "peak fit failed"
                    None
        // fit of the ion at toExtract index i, then its traces are released
        let results : QuantificationResult option [] = Array.create toExtract.Length None
        let chosenCenter = Dictionary<int,float>()
        let fitAt (i: int) (a: AnalysedIon) =
            if a.Ion.Row < n && decoyAlign.[a.Ion.Row].IsSome then chosenCenter.[a.Ion.Row] <- (fst (rank a).[0]).Center
            results.[i] <- fit a
            analysed.[i] <- None
        let newPeptides = HashSet<int>()
        let mutable newCandidates = 0
        let mutable newIons = 0
        let mutable newDecoys = 0
        if eager then
            model <- trainModel featureSet
        let analyseAndFit (i: int) (raw: RawXic) (ionMobilityCentre: float) =
            analyse i raw ionMobilityCentre
            match analysed.[i] with
            | Some a ->
                if a.Ion.Row < n then
                    newPeptides.Add a.Ion.Align.PepSequenceID |> ignore
                    newIons <- newIons + 1
                    newCandidates <- newCandidates + a.Candidates.Length
                else
                    newDecoys <- newDecoys + 1
                if eager then fitAt i a
            | None -> ()
        let newSpectraRead = sweepOver newIdx analyseAndFit
        let nIons = analysedOwn.Length + newIons
        logger.Info (sprintf "%s: XICs extracted from %i + %i MS1 spectrum reads (own ions + new ions), %i ions with candidates, %.1f candidates per ion, %i decoys with candidates (%.1f s)" runName ownSpectraRead newSpectraRead nIons (if nIons = 0 then nan else float ((analysedOwn |> Array.sumBy (fun a -> a.Candidates.Length)) + newCandidates) / float nIons) newDecoys (elapsed ()))
        if not eager then
            featureSet <-
                withIonMobilityFeature (
                    ownHasIm
                    || (newIdx |> Array.exists (fun i -> toExtract.[i].Row < n && (match analysed.[i] with Some a -> isFinite a.Context.PredictedIonMobility | None -> false)))
                    )
            model <- trainModel featureSet
        let selectionModel = selectionModelOf model
        // cross-fitting: an ion the run quantified itself is scored by the model trained on the other folds (folds of peptides,
        // PepSequenceID), so its selection and probability are out-of-fold like those of the new ions, which get the model of all
        // training ions
        foldOf <- assignFolds processParams.Seed selectionFolds (Seq.append (analysedOwn |> Seq.map (fun a -> a.Ion.Align.PepSequenceID)) newPeptides)
        foldModels <-
            match model with
            | None -> [||]
            | Some _ ->
                Array.init selectionFolds (fun f ->
                    match fitOn featureSet (training |> Array.filter (fun (a, _) -> foldOf.[a.Ion.Align.PepSequenceID] <> f)) with
                    | Some (m, _, _) -> Some m
                    | None ->
                        logger.Warn (sprintf "%s: the training candidates outside fold %i are all of one class; the ions of fold %i the run quantified itself are chosen with the fixed score." runName f f)
                        None
                    )
        let mutable overlap = analysedOwn |> Array.choose (fun a -> ownOf a |> Option.map (fun qp -> a, qp))
        let chosenCorrect (sel: (AnalysedIon*QuantificationResult) []) =
            if sel.Length = 0 then nan
            else 100. * (sel |> Array.filter (fun (a, qp) -> isCorrect qp (a.Candidates |> Array.maxBy (score a))) |> Array.length |> float) / float sel.Length
        logger.Info (sprintf "%s: chosen candidate within 1 sd of the own PSM peak: %.1f %% of %i clean overlap ions, %.1f %% of %i overlap ions (%s, out-of-fold)" runName (chosenCorrect training) training.Length (chosenCorrect overlap) overlap.Length selectionModel)
        // fits of the ions still waiting: the overlap ions (all ions if the second sweep could change the features); their traces
        // are released one by one
        analysedOwn <- [||]
        training <- [||]
        overlap <- [||]
        for i in 0 .. toExtract.Length - 1 do
            match analysed.[i] with
            | Some a -> fitAt i a
            | None -> ()
        let quantResults = [| for i in 0 .. prepared.Length - 1 do match results.[i] with Some q -> yield q | None -> () |]
        logger.Info (
            sprintf "%s: %i of %i .align rows quantified%s; not quantified: %s (%.1f s)"
                runName quantResults.Length alignments.Length
                (if performLabeled then sprintf ", %i without quantified labeled partner" partnerFailures else "")
                drops.Summary (elapsed ())
            )
        match processParams.XicExtraction.IonMobilityWindow with
        | Domain.Window.Fixed halfWidth when centreTally.Count > 0 ->
            logger.Info (sprintf "%s: ion mobility: %s" runName (centreTally.Summary halfWidth))
        | _ -> ()

        // output: the PSM based rows with a fitted peak, then the alignment rows; of several rows of an ion and source
        // the one with the smallest |Difference_SearchRT_FittedRT|
        let outputRows =
            Array.append psm.PsmRowsHeld quantResults
            |> Array.groupBy (fun qp -> qp.StringSequence, qp.ModSequenceID, qp.Charge, qp.GlobalMod, qp.QuantificationSource)
            |> Array.map (fun (_, rows) ->
                rows |> Array.minBy (fun qp -> let d = abs (getTargetApexDifference qp) in if isFinite d then d else infinity)
                )
        let repeated = psm.PsmFitted + quantResults.Length - (psm.PsmWritten + outputRows.Length)
        if psm.PsmFitted < psm.PsmCount || repeated > 0 then
            logger.Info (sprintf "%s: not written: %i PSM based rows without fitted peak, %i repeated rows of an ion" runName (psm.PsmCount - psm.PsmFitted) repeated)
        if psm.PsmWrittenEarly then
            File.AppendAllLines(outTmpPath, outputRows |> SeqIO'.csv "\t" false false)
            File.Move(outTmpPath, outFilePath, true)
        else
            File.WriteAllLines(outFilePath, outputRows |> SeqIO'.csv "\t" true false)
        let featuresOf (row: int) =
            match features.[row] with
            | Some f -> f
            | None ->
                let align = alignOf row
                {
                    StringSequence       = align.StringSequence
                    GlobalMod            = align.GlobalMod
                    Charge               = align.Charge
                    PepSequenceID        = align.PepSequenceID
                    ModSequenceID        = align.ModSequenceID
                    CandidateProbability = nan
                    CandidateMargin      = nan
                    NCandidates          = 0
                    AbsPpmMono           = nan
                    IsoCorr              = nan
                    IsoRatioLog2         = nan
                    ShapeCorr            = nan
                    RtZ                  = nan
                    ImZ                  = nan
                    LogRelHeight         = nan
                    AbsLog2FwhmRatio     = nan
                    CandidateScore       = nan
                    SelectionModel       = selectionModel
                    CandidateScanTime    = nan
                    ScanTimeWindow       = nan
                    MzWindow             = mzWindow
                    IonMobilityWindow    = nan
                    Status               = status.[row]
                }
        let featureRows = Array.init n featuresOf
        File.WriteAllLines(featureFilePath, featureRows |> SeqIO'.csv "\t" true false)
        logger.Info (sprintf "%s: results written to %s and %s (%.1f s)" runName outFilePath featureFilePath (elapsed ()))

        // decoy transfers: left out if the ion's own scan time is inside the decoy XIC; the others are fitted and written to <run>.aligndecoy,
        // also those without fit: the .quant columns, PredictedScanTime of the decoy and the .alignfeatures columns without the ion key
        if processParams.DecoyFraction > 0. then
            let ownPeakInside (ion: PreparedIon) =
                let row = ion.Row - n
                let scanTime = match chosenCenter.TryGetValue row with | true, t -> t | _ -> alignments.[row].PredictedScanTime
                abs (scanTime - ion.Align.PredictedScanTime) <= processParams.XicExtraction.ScanTimeWindow + decoyMargin
            let excluded = decoys |> Array.filter ownPeakInside |> Array.map (fun ion -> ion.Row) |> HashSet
            for _ in excluded do decoySkips.Add "own peak inside the decoy XIC"
            let decoyResults = dict [ for i in prepared.Length .. toExtract.Length - 1 do if not (excluded.Contains toExtract.[i].Row) then match results.[i] with Some q -> yield toExtract.[i].Row, q | None -> () ]
            let decoyRows =
                decoys
                |> Array.filter (fun ion -> not (excluded.Contains ion.Row))
                |> Array.map (fun ion ->
                    let row =
                        match decoyResults.TryGetValue ion.Row with
                        | true, q -> q
                        | _ -> createRow ion nan (emptySide ion.LockMz None) (emptySide ion.PartnerLockMz None) nan nan nan
                    row, ion.Align.PredictedScanTime, featuresOf ion.Row
                    )
            let quantLines = decoyRows |> Array.map (fun (q, _, _) -> q) |> SeqIO'.csv "\t" true false |> Array.ofSeq
            let scanTimes = Array.append [|"PredictedScanTime"|] (decoyRows |> Array.map (fun (_, t, _) -> t.ToString("R", Globalization.CultureInfo.InvariantCulture)))
            let featureLines = decoyRows |> Array.map (fun (_, _, f) -> f) |> SeqIO'.csv "\t" true false |> Array.ofSeq
            let lines = Array.map3 (fun (q: string) t (f: string) -> String.concat "\t" [| q; t; String.concat "\t" (f.Split('\t').[5 ..]) |]) quantLines scanTimes featureLines
            File.WriteAllLines(decoyFilePath, lines)
            logger.Info (
                sprintf "%s: %i decoy transfers written to %s (%s); skipped: %s (%.1f s)"
                    runName decoyRows.Length decoyFilePath
                    (decoyRows |> Array.countBy (fun (_, _, f) -> f.Status) |> Array.map (fun (s, c) -> sprintf "%s: %i" s c) |> String.concat ", ")
                    decoySkips.Summary (elapsed ())
                )
        // AlignmentBasedQuantStatistics reads <run>.aligndecoy whenever it exists, so the decoys of an earlier run would be used
        elif File.Exists decoyFilePath then
            File.Delete decoyFilePath
            logger.Info (sprintf "%s: DecoyFraction is 0; %s of an earlier run deleted" runName decoyFilePath)
