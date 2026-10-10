namespace ProteomIQon

open System
open System.IO
open System.Globalization
open System.Collections.Generic
open FSharpAux.IO.SchemaReader
open Microsoft.ML
open Microsoft.ML.Data
open Microsoft.ML.Trainers.FastTree
open Microsoft.ML.Trainers.LightGbm
open Plotly.NET
open Dto


/// Features of a transferred ion, computed from its AlignmentBasedQuantification row and its .align row only.
module AlignmentFeatures =

    let isFinite (v: float) = not (Double.IsNaN v || Double.IsInfinity v)

    /// Distance of neighbouring isotope peaks (Da).
    let c13 = 1.0033548378

    /// Median of the values, nan for an empty array.
    let median (values: float []) =
        if values.Length = 0 then
            nan
        else
            let sorted = Array.sort values
            let n = sorted.Length
            if n % 2 = 1 then sorted.[n / 2] else 0.5 * (sorted.[n / 2 - 1] + sorted.[n / 2])

    /// Index of the first value >= v in the sorted array.
    let lowerBound (sorted: float []) (v: float) =
        let mutable lo = 0
        let mutable hi = sorted.Length
        while lo < hi do
            let mid = (lo + hi) / 2
            if sorted.[mid] < v then lo <- mid + 1 else hi <- mid
        lo

    /// Index of the first value > v in the sorted array.
    let upperBound (sorted: float []) (v: float) =
        let mutable lo = 0
        let mutable hi = sorted.Length
        while lo < hi do
            let mid = (lo + hi) / 2
            if sorted.[mid] <= v then lo <- mid + 1 else hi <- mid
        lo

    /// Population standard deviation.
    let private sd (values: float []) =
        let m = Array.average values
        values |> Array.averageBy (fun v -> (v - m) * (v - m)) |> sqrt

    let private zNorm (values: float []) =
        let m = Array.average values
        let s = sd values
        if s > 0. then values |> Array.map (fun v -> (v - m) / s) else values |> Array.map (fun v -> v - m)

    /// Pearson correlation, nan if one of the series is constant.
    let private pearson (x: float []) (y: float []) =
        let sx = sd x
        let sy = sd y
        if not (sx > 0. && sy > 0.) then
            nan
        else
            let mx = Array.average x
            let my = Array.average y
            (Array.fold2 (fun acc a b -> acc + (a - mx) * (b - my)) 0. x y / float x.Length) / (sx * sy)

    /// Square root of the minimal cumulative squared distance of a dynamic time warping without constraints,
    /// as DTW'.distance with all options None.
    let dtwDistance (x: float []) (y: float []) =
        let c = y.Length
        let mutable previous = Array.create (c + 1) infinity
        let mutable current = Array.create (c + 1) infinity
        previous.[0] <- 0.
        for i = 0 to x.Length - 1 do
            current.[0] <- infinity
            for j = 0 to c - 1 do
                let d = (x.[i] - y.[j]) * (x.[i] - y.[j])
                current.[j + 1] <- d + min previous.[j] (min previous.[j + 1] current.[j])
            let tmp = previous
            previous <- current
            current <- tmp
        sqrt previous.[c]

    /// The source trace sorted by scan time.
    let private sortTrace (rt: float []) (intensity: float []) =
        let order = Array.init rt.Length id |> Array.sortBy (fun i -> rt.[i])
        Array.map (fun i -> rt.[i]) order, Array.map (fun i -> intensity.[i]) order

    /// Linear interpolation of the sorted trace at x, 0 outside of the trace.
    let private interpolateOrZero (xs: float []) (ys: float []) (x: float) =
        let n = xs.Length
        if x < xs.[0] || x > xs.[n - 1] then 0.
        elif x = xs.[n - 1] then ys.[n - 1]
        else
            let hi = upperBound xs x
            let lo = hi - 1
            ys.[lo] + (x - xs.[lo]) / (xs.[hi] - xs.[lo]) * (ys.[hi] - ys.[lo])

    let private validTrace (rt: float []) (intensity: float []) =
        rt.Length >= 5 && rt.Length = intensity.Length

    /// Correlation of the target XIC with the source trace shifted by shift, on the target scan times.
    let traceCorrelation (rtT: float []) (itT: float []) (rtS: float []) (itS: float []) (shift: float) =
        if not (validTrace rtT itT && validTrace rtS itS && isFinite shift) then
            nan
        else
            let xs, ys = sortTrace rtS itS
            let xs = xs |> Array.map (fun x -> x + shift)
            pearson (rtT |> Array.map (interpolateOrZero xs ys)) itT

    /// Correlation of the target XIC within +-3 sd of the fitted peak with the source trace shifted apex on apex.
    let localCorrelation (rtT: float []) (itT: float []) (rtS: float []) (itS: float []) (rtFit: float) (rtSource: float) (sd: float) =
        if not (validTrace rtT itT && validTrace rtS itS && isFinite sd && sd > 0. && isFinite rtFit && isFinite rtSource) then
            nan
        else
            let inside = Array.init rtT.Length id |> Array.filter (fun i -> abs (rtT.[i] - rtFit) <= 3. * sd)
            if inside.Length < 4 then
                nan
            else
                let xs, ys = sortTrace rtS itS
                let xs = xs |> Array.map (fun x -> x + rtFit - rtSource)
                pearson (inside |> Array.map (fun i -> interpolateOrZero xs ys rtT.[i])) (inside |> Array.map (fun i -> itT.[i]))

    /// Max XIC intensity outside +-3 sd of the fitted peak divided by the max inside, and the number of local maxima
    /// of the 3 point moving average above 20 % of its max.
    let xicCompetition (rtT: float []) (itT: float []) (rtFit: float) (sd: float) =
        if not (itT.Length >= 5 && rtT.Length = itT.Length && isFinite rtFit && isFinite sd && sd > 0.) then
            nan, nan
        else
            let mutable inside = 0.
            let mutable outside = 0.
            for i = 0 to rtT.Length - 1 do
                if abs (rtT.[i] - rtFit) <= 3. * sd then inside <- max inside itT.[i] else outside <- max outside itT.[i]
            let competition = if inside > 0. then outside / inside else nan
            let n = itT.Length
            let smoothed =
                Array.init n (fun i ->
                    let left = if i > 0 then itT.[i - 1] else 0.
                    let right = if i < n - 1 then itT.[i + 1] else 0.
                    (left + itT.[i] + right) / 3.
                    )
            let threshold = 0.2 * Array.max smoothed
            let mutable peaks = 0
            for i = 1 to n - 2 do
                if smoothed.[i] > smoothed.[i - 1] && smoothed.[i] >= smoothed.[i + 1] && smoothed.[i] > threshold then
                    peaks <- peaks + 1
            competition, float peaks

    /// Intensity weighted |ppm| of the isotope peaks against the calibrated theoretical m/z, ppm of the monoisotopic
    /// peak, number of isotope peaks and cosine of the target and source isotope pattern.
    let isotopeFeatures (mzs: float []) (intensities: float []) (quantMz: float) (charge: int) (sourceMzs: float []) (sourceIntensities: float []) =
        if mzs.Length = 0 then
            nan, nan, 0., nan
        elif mzs.Length <> intensities.Length || not (isFinite quantMz) || charge <= 0 then
            nan, nan, float mzs.Length, nan
        else
            let z = float charge
            let isotope mz = Math.Round((mz - quantMz) * z / c13)
            let ks = mzs |> Array.map isotope
            let ppm = Array.map2 (fun mz k -> let theo = quantMz + k * c13 / z in (mz - theo) / theo * 1e6) mzs ks
            let total = Array.sum intensities
            let weights = if total > 0. then intensities |> Array.map (fun v -> v / total) else Array.create mzs.Length (1. / float mzs.Length)
            let ppmWeighted = Array.fold2 (fun acc w p -> acc + w * abs p) 0. weights ppm
            let ppmMono = match Array.tryFindIndex ((=) 0.) ks with Some i -> ppm.[i] | None -> nan
            let cosine =
                if sourceMzs.Length = 0 || sourceMzs.Length <> sourceIntensities.Length then
                    nan
                else
                    let sourceKs = sourceMzs |> Array.map isotope
                    let byIsotope (ks: float []) (values: float []) k =
                        Array.fold2 (fun acc ki v -> if ki = k then acc + v else acc) 0. ks values
                    let all = Array.append ks sourceKs |> Array.filter isFinite |> Array.distinct
                    let a = all |> Array.map (byIsotope ks intensities)
                    let b = all |> Array.map (byIsotope sourceKs sourceIntensities)
                    let norm v = v |> Array.sumBy (fun x -> x * x) |> sqrt
                    if norm a > 0. && norm b > 0. then Array.fold2 (fun acc x y -> acc + x * y) 0. a b / (norm a * norm b) else nan
            ppmWeighted, ppmMono, float mzs.Length, cosine

    /// The NORATIO features in this order, then the intensity ratio features.
    let baseNames =
        [|
            "RtResid"; "AbsRtResid"; "RtResidSd"; "LogSd"; "WidthRatio"; "LogSourceQuant"
            "LogSourceApex"; "RelSeo"; "KLDivObserved"; "KLDivCorrected"; "PpmWeighted"; "AbsPpmMono"; "NIsotopes"; "IsotopeCosine"
            "Dtw"; "TraceCorr"; "LocalCorr"; "ApexFracXic"; "XicLength"; "SourceScanTime"; "Mz"; "Charge"; "XicCompetition"
            "XicNPeaks"; "LocalRtScale"; "RtResidLocal"
        |]

    let ratioNames = [| "Log2QuantRatio"; "Log2ApexRatio"; "LogTargetQuant"; "LogTargetApex" |]

    let ionMobilityName = "AbsIonMobilityResidSd"

    /// Index of LocalRtScale and RtResidLocal, filled in once all ions of the run are known.
    let localRtScaleIndex = Array.IndexOf(baseNames, "LocalRtScale")
    let rtResidLocalIndex = Array.IndexOf(baseNames, "RtResidLocal")
    let absRtResidIndex = Array.IndexOf(baseNames, "AbsRtResid")

    /// The light or heavy values of a quantification row, by GlobalMod.
    type Side =
        {
            Quant           : float
            Apex            : float
            Seo             : float
            Params          : float []
            KLObserved      : float
            KLCorrected     : float
            QuantMz         : float
            IsotopeMz       : float []
            IsotopeCorrected: float []
            RtTrace         : float []
            IntensityTrace  : float []
        }

    let side (q: QuantificationResult) =
        if q.GlobalMod = 0 then
            {
                Quant = q.Quant_Light; Apex = q.MeasuredApex_Light; Seo = q.Seo_Light; Params = q.Params_Light
                KLObserved = q.KLDiv_Observed_Theoretical_Light
                KLCorrected = q.KLDiv_CorrectedObserved_Theoretical_Light; QuantMz = q.QuantMz_Light
                IsotopeMz = q.IsotopicPatternMz_Light; IsotopeCorrected = q.IsotopicPatternIntensity_Corrected_Light
                RtTrace = q.RtTrace_Light; IntensityTrace = q.IntensityTrace_Corrected_Light
            }
        else
            {
                Quant = q.Quant_Heavy; Apex = q.MeasuredApex_Heavy; Seo = q.Seo_Heavy; Params = q.Params_Heavy
                KLObserved = q.KLDiv_Observed_Theoretical_Heavy
                KLCorrected = q.KLDiv_CorrectedObserved_Theoretical_Heavy; QuantMz = q.QuantMz_Heavy
                IsotopeMz = q.IsotopicPatternMz_Heavy; IsotopeCorrected = q.IsotopicPatternIntensity_Corrected_Heavy
                RtTrace = q.RtTrace_Heavy; IntensityTrace = q.IntensityTrace_Corrected_Heavy
            }

    /// Fitted scan time and sd of the HULQ peak parameters.
    let scanTimeAndSd (peakParams: float []) =
        if peakParams.Length > 2 then peakParams.[1], peakParams.[2] else nan, nan

    /// Base and ratio features (baseNames @ ratioNames) of a transferred ion; LocalRtScale and RtResidLocal are nan here.
    let compute (target: QuantificationResult) (align: AlignmentResultV2) =
        let t = side target
        let rt, sd = scanTimeAndSd t.Params
        let qs = align.Quant_SourceFile
        let aps = align.ApexIntensity_SourceFile
        let rtS = align.RtTrace_SourceFile
        let itS = align.IntensityTrace_SourceFile
        let rtResid = rt - align.PredictedScanTime
        let dtw = if t.IntensityTrace.Length > 0 && itS.Length > 0 then dtwDistance (zNorm t.IntensityTrace) (zNorm itS) else nan
        let ppmWeighted, ppmMono, nIsotopes, isotopeCosine =
            isotopeFeatures t.IsotopeMz t.IsotopeCorrected t.QuantMz target.Charge align.IsotopicPatternMz_SourceFile align.IsotopicPatternIntensity_Observed_SourceFile
        let competition, nPeaks = xicCompetition t.RtTrace t.IntensityTrace rt sd
        let xicMax = if t.IntensityTrace.Length > 0 then Array.max t.IntensityTrace else nan
        [|
            rtResid
            abs rtResid
            abs rtResid / sd
            log10 sd
            log ((t.Quant / t.Apex) / (qs / aps)) / log 2.
            log10 qs
            log10 aps
            t.Seo / t.Apex
            t.KLObserved
            t.KLCorrected
            ppmWeighted
            abs ppmMono
            nIsotopes
            isotopeCosine
            dtw
            traceCorrelation t.RtTrace t.IntensityTrace rtS itS (align.PredictedScanTime - align.ScanTime_SourceFile)
            localCorrelation t.RtTrace t.IntensityTrace rtS itS rt align.ScanTime_SourceFile sd
            t.Apex / xicMax
            float t.IntensityTrace.Length
            align.ScanTime_SourceFile
            align.Mz
            float target.Charge
            competition
            nPeaks
            nan
            nan
            // intensity ratio features
            log (t.Quant / qs) / log 2.
            log (t.Apex / aps) / log 2.
            log10 t.Quant
            log10 t.Apex
        |]

    /// |IonMobility - PredictedIonMobility| / PredictedIonMobilitySd, nan without ion mobility.
    let ionMobilityResidual (target: QuantificationResult) (align: AlignmentResultV2) =
        let sd = align.PredictedIonMobilitySd
        if isFinite target.IonMobility && isFinite align.PredictedIonMobility && isFinite sd && sd > 0. then
            abs (target.IonMobility - align.PredictedIonMobility) / sd
        else
            nan

    /// Median |own PSM scan time - predicted scan time| of the reference ions predicted within [x - 2, x + 2) minutes,
    /// nan with fewer than 5 such ions.
    let localRtScale (referencePredicted: float []) (referenceAbsError: float []) =
        let order = Array.init referencePredicted.Length id |> Array.sortBy (fun i -> referencePredicted.[i])
        let predicted = order |> Array.map (fun i -> referencePredicted.[i])
        let absError = order |> Array.map (fun i -> referenceAbsError.[i])
        fun (x: float) ->
            let lo = lowerBound predicted (x - 2.)
            let hi = lowerBound predicted (x + 2.)
            if hi - lo >= 5 then median absError.[lo .. hi - 1] else nan


/// Gradient boosted trees with per-fold Platt calibration, and q-values from the posterior error probabilities.
module AlignmentScoring =

    open AlignmentFeatures

    [<CLIMutable>]
    type MLRow =
        {
            Label    : bool
            Features : float32 []
        }

    /// Per feature median of the finite training values (0 if there are none).
    let trainingMedians (rows: float [] []) (nFeatures: int) =
        Array.init nFeatures (fun j ->
            let values = rows |> Array.map (fun row -> row.[j]) |> Array.filter isFinite
            if values.Length = 0 then 0. else median values
            )

    /// Replaces nan and infinite values by the training median and appends a 0/1 missing flag per feature.
    let impute (medians: float []) (row: float []) =
        let n = medians.Length
        let out = Array.zeroCreate (2 * n)
        for j = 0 to n - 1 do
            if isFinite row.[j] then
                out.[j] <- float32 row.[j]
            else
                out.[j] <- float32 medians.[j]
                out.[n + j] <- 1.f
        out

    /// Fits p = 1 / (1 + exp(-(a * score + b))) by maximum likelihood (Newton steps, halved while the loss does not decrease).
    let fitPlatt (scores: float []) (labels: bool []) =
        let loss (a: float) (b: float) =
            Array.fold2 (fun acc s y ->
                let m = a * s + b
                let m = if y then m else -m
                // log(1 + exp(-m)), stable
                acc + (if m > 0. then log (1. + exp -m) else -m + log (1. + exp m))
                ) 0. scores labels
        let mutable a = 1.
        let mutable b = 0.
        let mutable current = loss a b
        let mutable iteration = 0
        let mutable converged = false
        while not converged && iteration < 100 do
            let mutable ga = 0.
            let mutable gb = 0.
            let mutable haa = 1e-12
            let mutable hab = 0.
            let mutable hbb = 1e-12
            for i = 0 to scores.Length - 1 do
                let p = 1. / (1. + exp (-(a * scores.[i] + b)))
                let r = p - (if labels.[i] then 1. else 0.)
                let w = p * (1. - p)
                ga <- ga + r * scores.[i]
                gb <- gb + r
                haa <- haa + w * scores.[i] * scores.[i]
                hab <- hab + w * scores.[i]
                hbb <- hbb + w
            let det = haa * hbb - hab * hab
            let da = (hbb * ga - hab * gb) / det
            let db = (haa * gb - hab * ga) / det
            let mutable step = 1.
            let mutable next = loss (a - da) (b - db)
            while next > current && step > 1e-6 do
                step <- step / 2.
                next <- loss (a - step * da) (b - step * db)
            if next > current || not (isFinite next) then
                converged <- true
            else
                a <- a - step * da
                b <- b - step * db
                converged <- current - next < 1e-10 * (1. + current)
                current <- next
            iteration <- iteration + 1
        a, b

    let platt (a: float, b: float) (score: float) = 1. / (1. + exp (-(a * score + b)))

    /// q-values: running mean of the posterior error probabilities 1 - p over the ions sorted by p (descending),
    /// made monotone; tied p share the q-value of the last of them.
    let qValues (p: float []) =
        let n = p.Length
        let order = Array.init n id |> Array.sortWith (fun i j -> compare p.[j] p.[i])
        let q = Array.zeroCreate n
        let mutable sum = 0.
        for k = 0 to n - 1 do
            sum <- sum + (1. - p.[order.[k]])
            q.[k] <- sum / float (k + 1)
        for k = n - 2 downto 0 do
            q.[k] <- min q.[k] q.[k + 1]
        // ties
        for k = n - 2 downto 0 do
            if p.[order.[k]] = p.[order.[k + 1]] then q.[k] <- q.[k + 1]
        let out = Array.zeroCreate n
        order |> Array.iteri (fun k i -> out.[i] <- q.[k])
        out

    /// Storey's estimate of the fraction of new ions whose analyte is absent: the fraction of new ions with p < lambda divided by the
    /// fraction of decoys with p < lambda, at most 1, and 1 if no decoy is below lambda.
    let absentFraction (lambda: float) (pNew: float []) (pDecoy: float []) =
        let below (p: float []) = float (p |> Array.filter (fun v -> v < lambda) |> Array.length) / float p.Length
        if pNew.Length = 0 || not (below pDecoy > 0.) then 1. else min 1. (below pNew / below pDecoy)

    /// q-values with the false discovery rate of absent analytes added. At threshold t = p of a new ion,
    /// FDR_absent(t) = pi * (decoys with p >= t / decoys) * newIonCount / (new ions with p >= t), and q = min(1, qWrong + FDR_absent),
    /// made monotone over the ions sorted by p. pDecoy contains the decoys without row as 0, newIonCount counts the new ions without row.
    let absentQValues (pi: float) (newIonCount: int) (pDecoy: float []) (pNew: float []) (qWrong: float []) =
        let sortedDecoy = Array.sort pDecoy
        let sortedNew = Array.sort pNew
        let atLeast (sorted: float []) t = sorted.Length - lowerBound sorted t
        let q =
            Array.init pNew.Length (fun i ->
                let t = pNew.[i]
                let fdrAbsent = pi * float (atLeast sortedDecoy t) / float sortedDecoy.Length * float newIonCount / float (atLeast sortedNew t)
                min 1. (qWrong.[i] + fdrAbsent)
                )
        let order = Array.init pNew.Length id |> Array.sortWith (fun i j -> compare pNew.[j] pNew.[i])
        for k = order.Length - 2 downto 0 do
            q.[order.[k]] <- min q.[order.[k]] q.[order.[k + 1]]
        q

    /// Area under the ROC curve (rank statistic, ties averaged).
    let auc (scores: float []) (labels: bool []) =
        let n = scores.Length
        let order = Array.init n id |> Array.sortBy (fun i -> scores.[i])
        let ranks = Array.zeroCreate n
        let mutable k = 0
        while k < n do
            let mutable l = k
            while l + 1 < n && scores.[order.[l + 1]] = scores.[order.[k]] do l <- l + 1
            let rank = float (k + l) / 2. + 1.
            for m = k to l do ranks.[order.[m]] <- rank
            k <- l + 1
        let nPos = labels |> Array.filter id |> Array.length |> float
        let nNeg = float n - nPos
        let rankSum = Array.fold2 (fun acc r y -> if y then acc + r else acc) 0. ranks labels
        (rankSum - nPos * (nPos + 1.) / 2.) / (nPos * nNeg)

    /// Assigns every peptide (PepSequenceID, so all modified forms of a sequence) to one of nFolds folds by a seeded permutation
    /// of the sorted ids.
    let assignFolds (seed: int) (nFolds: int) (peptides: seq<int>) =
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

    /// Tree model of one fold with the Platt calibration fitted on the held-out ions of the fold.
    type FoldModel =
        {
            Medians     : float []
            Model       : ITransformer
            Calibration : float * float
        }

    let private schema (nFeatures: int) =
        let schema = SchemaDefinition.Create(typeof<MLRow>)
        schema.["Features"].ColumnType <- VectorDataViewType(NumberDataViewType.Single, 2 * nFeatures)
        schema

    /// Raw scores of the rows.
    let rawScores (mlContext: MLContext) (fold: FoldModel) (rows: float [] []) =
        if rows.Length = 0 then
            [||]
        else
            let data = mlContext.Data.LoadFromEnumerable(rows |> Array.map (fun row -> { Label = false; Features = impute fold.Medians row }), schema fold.Medians.Length)
            fold.Model.Transform(data).GetColumn<float32>("Score")
            |> Seq.map float
            |> Array.ofSeq

    /// Trains the gradient boosted trees of the parameters, single threaded and seeded.
    let trainModel (mlContext: MLContext) (p: Domain.AlignmentBasedQuantStatisticsV2Params) (rows: float [] []) (labels: bool []) =
        let nFeatures = rows.[0].Length
        let medians = trainingMedians rows nFeatures
        let data = mlContext.Data.LoadFromEnumerable(Array.map2 (fun row label -> { Label = label; Features = impute medians row }) rows labels, schema nFeatures)
        let model =
            match p.Trainer with
            | Domain.AlignmentScoreTrainer.FastTree ->
                let options = FastTreeBinaryTrainer.Options()
                options.LabelColumnName <- "Label"
                options.FeatureColumnName <- "Features"
                options.NumberOfLeaves <- p.NumberOfLeaves
                options.NumberOfTrees <- p.NumberOfTrees
                options.MinimumExampleCountPerLeaf <- p.MinimumExampleCountPerLeaf
                options.LearningRate <- p.LearningRate
                // one core per file
                options.NumberOfThreads <- Nullable 1
                mlContext.BinaryClassification.Trainers.FastTree(options).Fit(data) :> ITransformer
            | Domain.AlignmentScoreTrainer.LightGbm ->
                let options = LightGbmBinaryTrainer.Options()
                options.LabelColumnName <- "Label"
                options.FeatureColumnName <- "Features"
                options.NumberOfLeaves <- Nullable p.NumberOfLeaves
                options.NumberOfIterations <- p.NumberOfTrees
                options.MinimumExampleCountPerLeaf <- Nullable p.MinimumExampleCountPerLeaf
                options.LearningRate <- Nullable p.LearningRate
                options.NumberOfThreads <- Nullable 1
                options.Seed <- Nullable p.Seed
                options.Silent <- true
                mlContext.BinaryClassification.Trainers.LightGbm(options).Fit(data) :> ITransformer
        medians, model


module AlignmentBasedQuantStatistics_v2 =

    open AlignmentFeatures
    open AlignmentScoring

    type Params = Domain.AlignmentBasedQuantStatisticsV2Params

    type IonKey =
        {
            StringSequence : string
            GlobalMod      : int
            Charge         : int
            PepSequenceID  : int
            ModSequenceID  : int
        }

    let keyOfQuant (q: QuantificationResult) =
        { StringSequence = q.StringSequence; GlobalMod = q.GlobalMod; Charge = q.Charge; PepSequenceID = q.PepSequenceID; ModSequenceID = q.ModSequenceID }

    let keyOfAlign (a: AlignmentResultV2) =
        { StringSequence = a.StringSequence; GlobalMod = a.GlobalMod; Charge = a.Charge; PepSequenceID = a.PepSequenceID; ModSequenceID = a.ModSequenceID }

    let keyName (k: IonKey) = sprintf "%s, GlobalMod %i, charge %i" k.StringSequence k.GlobalMod k.Charge

    /// Ground truth of an overlap ion: the alignment peak against the own PSM peak of the run.
    type Truth =
        | Positive
        | Ambiguous
        | Negative
        /// The own PSM quantification has no finite scan time, sd or quant; the ion is scored like a new ion.
        | NoTruth

    ///
    type Ion =
        {
            Key       : IonKey
            /// Row of AlignmentBasedQuantification (QuantificationSource Alignment).
            Row       : QuantificationResult
            Align     : AlignmentResultV2
            /// None for new ions (no own PSM quantification in the run).
            Truth     : Truth option
            /// Fitted scan time of the own PSM peak, nan without a finite one.
            OwnScanTime : float
            /// |own PSM scan time - PredictedScanTime|, nan without a finite own scan time.
            OwnRtError : float
            /// baseNames @ ratioNames
            Values    : float []
            IonMobility : float
            /// Values of the .alignfeatures columns, by column name.
            Sidecar   : Dictionary<string,float>
        }

    ///
    type Run =
        {
            Name       : string
            QuantPath  : string
            /// All rows of the PSM based quantification, written to the output unchanged.
            PsmRows    : QuantificationResult []
            Ions       : Ion []
            /// The values of the SelectionModel column of the .alignfeatures file.
            SelectionModels : string []
            /// New ions of the .align file without alignment row.
            NewWithoutRow  : int
            /// The run has a .aligndecoy file (read unless AbsentAnalyteFdr is Off).
            DecoyFile      : bool
            /// Decoy transfers with a fitted peak.
            Decoys         : Ion []
            /// Decoy transfers without fitted peak.
            UnscoredDecoys : int
        }

    let private keyColumns = [| "StringSequence"; "GlobalMod"; "Charge"; "PepSequenceID"; "ModSequenceID" |]

    let private parseFloat (s: string) =
        let s = s.Trim('"', ' ')
        match Double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture) with
        | true, v -> Some v
        | _ ->
            if s = "" then Some nan else None

    /// Reads the numeric columns of the .alignfeatures file of AlignmentBasedQuantification_v2 and the distinct values of its
    /// SelectionModel column.
    let readSidecar (logger: NLog.Logger) (path: string) =
        let lines = File.ReadAllLines path
        let header = lines.[0].Split('\t') |> Array.map (fun s -> s.Trim('"', ' '))
        let keyIndices = keyColumns |> Array.map (fun c -> Array.IndexOf(header, c))
        let rows =
            lines.[1 ..]
            |> Array.filter (fun line -> line.Trim() <> "")
            |> Array.map (fun line -> line.Split('\t'))
        let candidates = [| 0 .. header.Length - 1 |] |> Array.filter (fun j -> keyIndices |> Array.contains j |> not)
        let parsed = candidates |> Array.map (fun j -> rows |> Array.map (fun row -> parseFloat row.[j]))
        let numeric = Array.zip candidates parsed |> Array.filter (fun (_, values) -> values |> Array.forall Option.isSome)
        let skipped = Array.zip candidates parsed |> Array.filter (fun (_, values) -> values |> Array.exists Option.isNone) |> Array.map (fun (j, _) -> header.[j])
        if skipped.Length > 0 then
            logger.Info (sprintf "%s: columns that are not numeric and not used: %s" path (String.concat ", " skipped))
        let byIon = Dictionary<IonKey,Dictionary<string,float>>(HashIdentity.Structural)
        rows
        |> Array.iteri (fun r row ->
            let key = { StringSequence = row.[keyIndices.[0]].Trim('"', ' '); GlobalMod = int row.[keyIndices.[1]]; Charge = int row.[keyIndices.[2]]; PepSequenceID = int row.[keyIndices.[3]]; ModSequenceID = int row.[keyIndices.[4]] }
            let values = Dictionary<string,float>()
            numeric |> Array.iter (fun (j, column) -> values.[header.[j]] <- (column.[r] |> Option.get))
            byIon.[key] <- values
            )
        let selectionIndex = Array.IndexOf(header, "SelectionModel")
        let selectionModels = rows |> Array.map (fun row -> row.[selectionIndex].Trim('"', ' ')) |> Array.distinct |> Array.sort
        let names = numeric |> Array.map (fun (j, _) -> header.[j])
        logger.Info (sprintf "%s: %i ions, %i numeric feature columns: %s; selection models: %s" path byIon.Count names.Length (String.concat ", " names) (String.concat ", " selectionModels))
        byIon, selectionModels

    let private readQuant (path: string) =
        Csv.CsvReader<QuantificationResult>().ReadFile(path,'\t',false,1)
        |> Array.ofSeq

    /// The .aligndecoy file has the .quant columns followed by further columns.
    let private readDecoyQuant (path: string) =
        Csv.CsvReader<QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(path,'\t',false,1)
        |> Array.ofSeq

    let private count (counts: SortedDictionary<string,int>) (reason: string) =
        counts.[reason] <- (match counts.TryGetValue reason with | true, c -> c + 1 | _ -> 1)

    let private formatCounts (counts: SortedDictionary<string,int>) =
        if counts.Count = 0 then "none" else counts |> Seq.map (fun kv -> sprintf "%s: %i" kv.Key kv.Value) |> String.concat ", "

    /// Truth of an overlap ion from its alignment row and its own PSM row.
    let truth (p: Params) (row: QuantificationResult) (own: QuantificationResult) =
        let rtOwn, sdOwn = scanTimeAndSd (side own).Params
        let quantOwn = (side own).Quant
        if not (isFinite rtOwn && isFinite sdOwn && sdOwn > 0. && isFinite quantOwn && quantOwn > 0.) then
            NoTruth
        else
            let rt, _ = scanTimeAndSd (side row).Params
            let z = abs (rt - rtOwn) / sdOwn
            let absLog2Ratio = abs (log ((side row).Quant / quantOwn) / log 2.)
            // nan compares false: alignment rows without fitted peak or quant are negative
            if z <= p.PositiveRtSdFactor && absLog2Ratio <= p.MaxAbsLog2QuantRatio then Positive
            elif z > p.PositiveRtSdFactor && z <= p.NegativeRtSdFactor && absLog2Ratio <= p.MaxAbsLog2QuantRatio then Ambiguous
            else Negative

    /// Reads the PSM quantification, the .align file and the AlignmentBasedQuantification output of a run, joins the
    /// alignment rows with their .align row and computes the features.
    let readRun (p: Params) (quantPath: string, alignPath: string, alignQuantPath: string) =
        let name = Path.GetFileNameWithoutExtension quantPath
        let logger = Logging.createLogger name
        let psmRows = readQuant quantPath
        let own = Dictionary<IonKey,QuantificationResult>(HashIdentity.Structural)
        let mutable psmDuplicates = 0
        for q in psmRows do
            let key = keyOfQuant q
            if own.ContainsKey key then
                psmDuplicates <- psmDuplicates + 1
                logger.Debug (sprintf "PSM quantification has several rows, the last one is the truth: %s" (keyName key))
            own.[key] <- q
        let alignRows =
            Csv.CsvReader<AlignmentResultV2>().ReadFile(alignPath,'\t',false,1)
            |> Array.ofSeq
        let align = Dictionary<IonKey,AlignmentResultV2>(HashIdentity.Structural)
        for a in alignRows do
            align.[keyOfAlign a] <- a
        let alignmentRows =
            readQuant alignQuantPath
            |> Array.filter (fun q -> q.QuantificationSource = QuantificationSource.Alignment)
        let sidecarPath = Path.Combine(Path.GetDirectoryName(Path.GetFullPath alignQuantPath), Path.GetFileNameWithoutExtension alignQuantPath + ".alignfeatures")
        let sidecarByIon, selectionModels = readSidecar logger sidecarPath
        let drops = SortedDictionary<string,int>()
        let seen = HashSet<IonKey>(HashIdentity.Structural)
        let ions =
            alignmentRows
            |> Array.choose (fun row ->
                let key = keyOfQuant row
                match align.TryGetValue key with
                | false, _ ->
                    logger.Debug (sprintf "Alignment row not used, no .align row: %s" (keyName key))
                    count drops "alignment rows without .align row"
                    None
                | true, _ when seen.Contains key ->
                    logger.Debug (sprintf "Alignment row not used, the ion has an earlier alignment row: %s" (keyName key))
                    count drops "repeated alignment rows"
                    None
                | true, a ->
                    seen.Add key |> ignore
                    let ownScanTime = match own.TryGetValue key with | true, o -> fst (scanTimeAndSd (side o).Params) | _ -> nan
                    Some {
                        Key         = key
                        Row         = row
                        Align       = a
                        Truth       = (match own.TryGetValue key with | true, o -> Some (truth p row o) | _ -> None)
                        OwnScanTime = ownScanTime
                        OwnRtError  = abs (ownScanTime - a.PredictedScanTime)
                        Values      = compute row a
                        IonMobility = ionMobilityResidual row a
                        Sidecar     = sidecarByIon.[key]
                    }
                )
        let notQuantified = align.Keys |> Seq.filter (fun key -> seen.Contains key |> not) |> Array.ofSeq
        for key in notQuantified do
            logger.Debug (sprintf ".align ion without alignment row of AlignmentBasedQuantification: %s" (keyName key))
        // decoy transfers of AlignmentBasedQuantification_v2: .quant columns, the shifted PredictedScanTime and the .alignfeatures columns
        let decoyPath = Path.Combine(Path.GetDirectoryName(Path.GetFullPath alignQuantPath), Path.GetFileNameWithoutExtension alignQuantPath + ".aligndecoy")
        let decoyFile = p.AbsentAnalyteFdr <> Domain.AbsentAnalyteFdr.Off && File.Exists decoyPath
        let decoys, unscoredDecoys =
            if decoyFile then
                let byIon, _ = readSidecar logger decoyPath
                let decoys =
                    readDecoyQuant decoyPath
                    |> Array.choose (fun row ->
                        let key = keyOfQuant row
                        match align.TryGetValue key with
                        | true, a ->
                            let values = byIon.[key]
                            let a = { a with PredictedScanTime = values.["PredictedScanTime"] }
                            Some { Key = key; Row = row; Align = a; Truth = None; OwnScanTime = nan; OwnRtError = nan; Values = compute row a; IonMobility = ionMobilityResidual row a; Sidecar = values }
                        | _ ->
                            logger.Debug (sprintf "Decoy not used, no .align row: %s" (keyName key))
                            None
                        )
                let scored = decoys |> Array.filter (fun ion -> (side ion.Row).Params.Length > 2)
                scored, decoys.Length - scored.Length
            else [||], 0
        let truthCounts = ions |> Array.countBy (fun ion -> match ion.Truth with Some t -> sprintf "%A" t | None -> "New")
        logger.Info (
            sprintf "%s: %i PSM rows (%i repeated), %i .align rows, %i alignment rows; %i ions scored or used (%s); not used: %s; .align ions without alignment row: %i"
                name psmRows.Length psmDuplicates alignRows.Length alignmentRows.Length ions.Length
                (truthCounts |> Array.map (fun (t, c) -> sprintf "%s %i" t c) |> String.concat ", ")
                (formatCounts drops) notQuantified.Length
            )
        if decoyFile then
            logger.Info (sprintf "%s: %i decoy transfers, %i with fitted peak" name (decoys.Length + unscoredDecoys) decoys.Length)
        {
            Name           = name
            QuantPath      = quantPath
            PsmRows        = psmRows
            Ions           = ions
            SelectionModels = selectionModels
            NewWithoutRow  = notQuantified |> Array.filter (fun key -> own.ContainsKey key |> not) |> Array.length
            DecoyFile      = decoyFile
            Decoys         = decoys
            UnscoredDecoys = unscoredDecoys
        }

    /// Fills in LocalRtScale and RtResidLocal: the scale of an overlap ion leaves out the overlap ions of its
    /// cross-validation fold (folds of peptides), new ions use all overlap ions of the run.
    let fillLocalRtScale (foldOf: Dictionary<int,int>) (nFolds: int) (run: Run) =
        let reference = run.Ions |> Array.filter (fun ion -> isFinite ion.OwnRtError)
        let scaleOf (ions: Ion []) = localRtScale (ions |> Array.map (fun ion -> ion.Align.PredictedScanTime)) (ions |> Array.map (fun ion -> ion.OwnRtError))
        let scaleAll = scaleOf reference
        let scaleWithout = Array.init nFolds (fun f -> scaleOf (reference |> Array.filter (fun ion -> foldOf.[ion.Key.PepSequenceID] <> f)))
        for ion in Array.append run.Ions run.Decoys do
            let scale = if isFinite ion.OwnRtError then scaleWithout.[foldOf.[ion.Key.PepSequenceID]] else scaleAll
            let s = scale ion.Align.PredictedScanTime
            ion.Values.[localRtScaleIndex] <- s
            ion.Values.[rtResidLocalIndex] <- ion.Values.[absRtResidIndex] / (if isFinite s then max s 0.01 else nan)

    /// Names of the features the model uses.
    type FeatureSet =
        {
            Names       : string []
            UseRatio    : bool
            UseIonMobility : bool
            /// Columns of the .alignfeatures file, without the prefix of their feature name.
            Sidecar     : string []
        }

    /// The columns of the .alignfeatures file that are used as features, named Sel_<column>. CandidateProbability and CandidateMargin
    /// only when all runs chose their candidates with the same selection model, because their unit depends on it. AbsPpmMono is a base
    /// feature already; the windows are constant per run and CandidateScanTime is the scan time.
    let sidecarColumns = [| "CandidateProbability"; "CandidateMargin"; "NCandidates"; "IsoCorr"; "IsoRatioLog2"; "ShapeCorr"; "RtZ"; "ImZ"; "LogRelHeight"; "AbsLog2FwhmRatio" |]

    let sidecarPrefix = "Sel_"

    /// AbsentAnalyteFdr Estimated: lambda of the absent fraction estimate, and the decoys a run needs to be calibrated on its own decoys
    /// (else on the decoys of all runs).
    let absentLambda = 0.5
    let minRunDecoys = 200

    let featureVector (fs: FeatureSet) (ion: Ion) =
        [|
            yield! ion.Values.[0 .. baseNames.Length - 1]
            if fs.UseRatio then yield! ion.Values.[baseNames.Length ..]
            if fs.UseIonMobility then yield ion.IonMobility
            for name in fs.Sidecar do
                yield ion.Sidecar.[name]
        |]

    let private formatFloat (v: float) = v.ToString("R", CultureInfo.InvariantCulture)

    /// Out-of-fold result of an overlap ion of a learning run.
    type OutOfFold =
        {
            Ion   : Ion
            Run   : string
            Fold  : int
            Score : float
            P     : float
        }

    let private calibrationTable (p: float []) (positive: bool []) =
        let edges = [| 0.; 0.5; 0.8; 0.9; 0.95; 0.98; 0.99; 0.995; 1.0001 |]
        [| 0 .. edges.Length - 2 |]
        |> Array.choose (fun b ->
            let inBin = Array.zip p positive |> Array.filter (fun (v, _) -> v >= edges.[b] && v < edges.[b + 1])
            if inBin.Length = 0 then None
            else
                Some (sprintf "[%.3g, %.3g) n %i predicted %.3f observed %.3f" edges.[b] (min 1. edges.[b + 1]) inBin.Length (inBin |> Array.averageBy fst) (inBin |> Array.averageBy (fun (_, y) -> if y then 1. else 0.)))
            )

    /// Ions at q-value <= cut and their fraction of negatives, with the ambiguous ions counted as negative and left out.
    let private fdpAt (q: float []) (negative: bool []) (ambiguous: bool []) (cut: float) =
        let selected = Array.init q.Length id |> Array.filter (fun i -> q.[i] <= cut)
        let fraction (indices: int []) = if indices.Length = 0 then nan else float (indices |> Array.filter (fun i -> negative.[i]) |> Array.length) / float indices.Length
        sprintf "q <= %g: %i ions, true FDP %.4f, without the ambiguous ions %.4f" cut selected.Length (fraction selected) (fraction (selected |> Array.filter (fun i -> not ambiguous.[i])))

    let private writeTestSet (outputDir: string) (fs: FeatureSet) (run: string) (oof: OutOfFold []) (q: Dictionary<IonKey * string,float>) =
        let header = Array.append [| "StringSequence"; "GlobalMod"; "Charge"; "PepSequenceID"; "ModSequenceID"; "Truth"; "Fold"; "RawScore"; "AlignmentScore"; "AlignmentQValue" |] fs.Names
        let lines =
            oof
            |> Array.map (fun o ->
                let k = o.Ion.Key
                Array.append
                    [| k.StringSequence; string k.GlobalMod; string k.Charge; string k.PepSequenceID; string k.ModSequenceID; sprintf "%A" o.Ion.Truth.Value; string o.Fold; formatFloat o.Score; formatFloat o.P; formatFloat q.[(k, run)] |]
                    (featureVector fs o.Ion |> Array.map formatFloat)
                |> String.concat "\t"
                )
        File.WriteAllLines(Path.Combine(outputDir, run + ".quanttestset"), Array.append [| String.concat "\t" header |] lines)

    /// Trains the alignment score on the overlap ions of the learning runs with cross-validation by peptide, scores the
    /// new ions of the runs to score, computes their q-values and writes one .quant file per run to score to outputDir.
    /// The q-values include the false discovery rate of absent analytes when AbsentAnalyteFdr is not Off, every run
    /// (learning runs too) has its .aligndecoy file, at least one decoy has a fitted peak and the runs to score have
    /// decoys with a fitted peak. testSet also writes the overlap ions with truth of each learning run and their
    /// out-of-fold scores (.quanttestset), diagCharts the AlignmentScore histogram. Each run is a
    /// (PSM .quant, .align, aligned .quant) triple; c is the number of files read and written in parallel.
    let scoreRuns (diagCharts: bool) (testSet: bool) (outputDir: string) (c: int) (p: Params) (scoringRuns: (string * string * string) []) (learningRuns: (string * string * string) []) =
        let logger = Logging.createLogger "AlignmentBasedQuantStatistics_v2"
        let runs =
            Array.append scoringRuns learningRuns
            |> Array.distinct
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> FSharpAux.PSeq.map (fun files -> files, readRun p files)
            |> dict
        let learning = learningRuns |> Array.map (fun files -> runs.[files])
        let scoring = scoringRuns |> Array.map (fun files -> runs.[files])
        // decoy transfers: negatives in training and the null distribution of absent analytes, if every run has them
        let useDecoys =
            if p.AbsentAnalyteFdr = Domain.AbsentAnalyteFdr.Off then
                false
            elif runs.Values |> Seq.forall (fun run -> run.DecoyFile) && runs.Values |> Seq.exists (fun run -> run.Decoys.Length > 0) then
                true
            else
                let missing = runs.Values |> Seq.filter (fun run -> not run.DecoyFile) |> Seq.map (fun run -> run.Name) |> Array.ofSeq
                logger.Warn (
                    sprintf "AbsentAnalyteFdr %A needs the decoy transfers of AlignmentBasedQuantification_v2 (DecoyFraction > 0, <run>.aligndecoy next to the aligned .quant); %s. AlignmentQValue does not include the false discovery rate of absent analytes."
                        p.AbsentAnalyteFdr (if missing.Length > 0 then "missing for " + String.concat ", " missing else "no decoy has a fitted peak")
                    )
                false
        let decoysOf (run: Run) = if useDecoys then run.Decoys else [||]
        // cross-validation folds by peptide (PepSequenceID), also used for the local scan time error scale
        let foldOf = assignFolds p.Seed p.NumberOfFolds (runs.Values |> Seq.collect (fun run -> Seq.append run.Ions (decoysOf run) |> Seq.map (fun ion -> ion.Key.PepSequenceID)))
        for run in runs.Values do
            fillLocalRtScale foldOf p.NumberOfFolds run

        // features
        let sidecar =
            let selectionModels = Seq.append learning scoring |> Seq.collect (fun run -> run.SelectionModels) |> Seq.distinct |> Array.ofSeq
            if selectionModels.Length <= 1 then sidecarColumns
            else
                logger.Info (sprintf "The runs chose their candidates with different selection models (%s); CandidateProbability and CandidateMargin are not used." (String.concat ", " selectionModels))
                sidecarColumns |> Array.filter (fun name -> name <> "CandidateProbability" && name <> "CandidateMargin")
        let trainingIons = learning |> Array.collect (fun run -> run.Ions |> Array.filter (fun ion -> ion.Truth.IsSome) |> Array.map (fun ion -> run.Name, ion))
        let useIonMobility = p.UseIonMobility && (trainingIons |> Array.exists (fun (_, ion) -> isFinite ion.IonMobility))
        let fs =
            {
                Names          = Array.concat [| baseNames; (if p.UseIntensityRatioFeatures then ratioNames else [||]); (if useIonMobility then [| ionMobilityName |] else [||]); sidecar |> Array.map (fun name -> sidecarPrefix + name) |]
                UseRatio       = p.UseIntensityRatioFeatures
                UseIonMobility = useIonMobility
                Sidecar        = sidecar
            }
        logger.Info (sprintf "%i features: %s" fs.Names.Length (String.concat ", " fs.Names))

        // cross-validated training on the overlap ions of the learning runs
        let truthCounts = trainingIons |> Array.countBy (fun (_, ion) -> ion.Truth.Value)
        logger.Info (sprintf "Overlap ions of the learning runs: %s" (truthCounts |> Array.map (fun (t, n) -> sprintf "%A %i" t n) |> String.concat ", "))
        let evaluated = trainingIons |> Array.filter (fun (_, ion) -> ion.Truth.Value <> NoTruth)
        let trainingDecoys = learning |> Array.collect decoysOf
        let mlContext = MLContext(seed = p.Seed)
        let stopwatch = Diagnostics.Stopwatch.StartNew()
        let calibrationOf (heldOut: ((string * Ion) * float) []) =
            let calibrationSet = heldOut |> Array.filter (fun ((_, ion), _) -> ion.Truth.Value <> Ambiguous)
            fitPlatt (calibrationSet |> Array.map snd) (calibrationSet |> Array.map (fun ((_, ion), _) -> ion.Truth.Value = Positive))
        let folds, heldOut =
            [| 0 .. p.NumberOfFolds - 1 |]
            |> Array.map (fun f ->
                let train = evaluated |> Array.filter (fun (_, ion) -> foldOf.[ion.Key.PepSequenceID] <> f && (ion.Truth.Value = Positive || ion.Truth.Value = Negative))
                let heldOut = evaluated |> Array.filter (fun (_, ion) -> foldOf.[ion.Key.PepSequenceID] = f)
                let labels = train |> Array.map (fun (_, ion) -> ion.Truth.Value = Positive)
                if labels |> Array.forall id || labels |> Array.forall not then
                    failwithf "Fold %i has only positive or only negative training ions; more overlap ions are needed." f
                // decoys are negatives in training, but not in the calibration, which is fitted on overlap ions only
                let decoys = trainingDecoys |> Array.filter (fun ion -> foldOf.[ion.Key.PepSequenceID] <> f)
                let medians, model =
                    trainModel mlContext p
                        (Array.append (train |> Array.map (snd >> featureVector fs)) (decoys |> Array.map (featureVector fs)))
                        (Array.append labels (Array.create decoys.Length false))
                let fold0 = { Medians = medians; Model = model; Calibration = (1., 0.) }
                let heldOutScores = rawScores mlContext fold0 (heldOut |> Array.map (snd >> featureVector fs))
                let calibration = calibrationOf (Array.zip heldOut heldOutScores)
                logger.Info (sprintf "Fold %i: %i training ions (%i positive), %i decoys, %i held out, Platt a %.4f b %.4f" f train.Length (labels |> Array.filter id |> Array.length) decoys.Length heldOut.Length (fst calibration) (snd calibration))
                { fold0 with Calibration = calibration }, Array.zip heldOut heldOutScores
                )
            |> Array.unzip
        // out-of-fold probabilities of the check: the scores of a fold are calibrated on the held-out scores of the other
        // folds, so that the check does not use a calibration fitted on the ions it checks
        let oof =
            heldOut
            |> Array.mapi (fun f items ->
                let crossCalibration = calibrationOf (heldOut |> Array.indexed |> Array.filter (fun (g, _) -> g <> f) |> Array.collect snd)
                logger.Info (sprintf "Fold %i: Platt on the held-out ions of the other folds a %.4f b %.4f" f (fst crossCalibration) (snd crossCalibration))
                items |> Array.map (fun ((run, ion), s) -> { Ion = ion; Run = run; Fold = f; Score = s; P = platt crossCalibration s })
                )
            |> Array.concat
        logger.Info (sprintf "Trained %i folds in %.1f s" p.NumberOfFolds stopwatch.Elapsed.TotalSeconds)

        // out-of-fold check
        let oofP = oof |> Array.map (fun o -> o.P)
        let oofPositive = oof |> Array.map (fun o -> o.Ion.Truth.Value = Positive)
        let trainable = oof |> Array.filter (fun o -> o.Ion.Truth.Value <> Ambiguous)
        let oofQ = qValues oofP
        logger.Info (sprintf "Out-of-fold AUC (positive vs negative) %.4f" (auc (trainable |> Array.map (fun o -> o.P)) (trainable |> Array.map (fun o -> o.Ion.Truth.Value = Positive))))
        for cut in [| 0.01; 0.05; 0.1 |] do
            logger.Info (sprintf "Out-of-fold check (ambiguous counted false): %s" (fdpAt oofQ (oofPositive |> Array.map not) (oof |> Array.map (fun o -> o.Ion.Truth.Value = Ambiguous)) cut))
        for line in calibrationTable oofP oofPositive do
            logger.Info (sprintf "Out-of-fold calibration %s" line)

        // scoring of the new ions and of the overlap ions without truth whose PSM row has no fitted scan time: their
        // alignment row replaces the PSM row. The other overlap ions without truth keep their PSM row and their alignment
        // row is not written, so that the output has one row per ion.
        let score (rows: float [] []) =
            let perFold = folds |> Array.map (fun fold -> rawScores mlContext fold rows |> Array.map (platt fold.Calibration))
            Array.init rows.Length (fun r -> perFold |> Array.averageBy (fun scores -> scores.[r]))
        // a decoy is scored by the model of its fold, which was trained without it; the decoys without fitted peak have p = 0
        let decoyScores (run: Run) =
            let decoys = decoysOf run
            let p = Array.zeroCreate (decoys.Length + (if useDecoys then run.UnscoredDecoys else 0))
            decoys
            |> Array.indexed
            |> Array.groupBy (fun (_, ion) -> foldOf.[ion.Key.PepSequenceID])
            |> Array.iter (fun (f, items) ->
                let scores = rawScores mlContext folds.[f] (items |> Array.map (snd >> featureVector fs))
                Array.iter2 (fun (i, _) s -> p.[i] <- platt folds.[f].Calibration s) items scores
                )
            p
        let runDecoyScores = scoring |> Array.map decoyScores
        let pooledDecoys = Array.concat runDecoyScores
        if useDecoys then
            logger.Info (sprintf "Decoy transfers of the runs to score: %i, median AlignmentScore %.4f (out-of-fold, without fitted peak 0)" pooledDecoys.Length (median pooledDecoys))
        let scored =
            Array.zip scoring runDecoyScores
            |> Array.map (fun (run, own) ->
                let newIons = run.Ions |> Array.filter (fun ion -> ion.Truth.IsNone || (ion.Truth = Some NoTruth && not (isFinite ion.OwnScanTime)))
                let keepPsmRow = run.Ions |> Array.filter (fun ion -> ion.Truth = Some NoTruth && isFinite ion.OwnScanTime) |> Array.length
                let pNew = score (newIons |> Array.map (featureVector fs))
                let qWrong = qValues pNew
                let runLogger = Logging.createLogger run.Name
                // false discovery rate of absent analytes from the decoys of the run (of all runs if the run has fewer than minRunDecoys)
                let qNew =
                    if not useDecoys then
                        qWrong
                    else
                        let pDecoy = if own.Length >= minRunDecoys then own else pooledDecoys
                        let newIonCount = newIons.Length + run.NewWithoutRow
                        let pi =
                            match p.AbsentAnalyteFdr with
                            | Domain.AbsentAnalyteFdr.Conservative -> 1.
                            | _ -> absentFraction absentLambda (Array.append pNew (Array.zeroCreate run.NewWithoutRow)) pDecoy
                        let qNew = absentQValues pi newIonCount pDecoy pNew qWrong
                        let passing (q: float []) cut = q |> Array.filter (fun v -> v <= cut) |> Array.length
                        let decoyPassing cut =
                            let threshold = Array.fold2 (fun acc v q -> if q <= cut then min acc v else acc) infinity pNew qNew
                            100. * float (pDecoy |> Array.filter (fun v -> v >= threshold) |> Array.length) / float pDecoy.Length
                        runLogger.Info (
                            sprintf "%s: absent analytes: %i decoys (%i with fitted peak)%s, %i new ions (%i without alignment row), absent fraction %.4f (%s); decoys passing q <= 0.01: %.2f %%, q <= 0.05: %.2f %%; new ions q <= 0.01: %i -> %i, q <= 0.05: %i -> %i"
                                run.Name own.Length (decoysOf run).Length (if own.Length >= minRunDecoys then "" else sprintf ", fewer than %i, the %i decoys of all runs are used" minRunDecoys pooledDecoys.Length)
                                newIonCount run.NewWithoutRow pi (sprintf "%A" p.AbsentAnalyteFdr)
                                (decoyPassing 0.01) (decoyPassing 0.05) (passing qWrong 0.01) (passing qNew 0.01) (passing qWrong 0.05) (passing qNew 0.05)
                            )
                        qNew
                runLogger.Info (
                    sprintf "%s: %i new ions (%i of them overlap ions without truth whose PSM row has no fitted scan time) scored, mean AlignmentScore %.4f, q <= 0.01: %i, q <= 0.05: %i, q <= 0.1: %i; %i overlap ions without truth keep their PSM row"
                        run.Name newIons.Length (newIons |> Array.filter (fun ion -> ion.Truth.IsSome) |> Array.length) (if newIons.Length = 0 then nan else Array.average pNew)
                        (qNew |> Array.filter (fun q -> q <= 0.01) |> Array.length) (qNew |> Array.filter (fun q -> q <= 0.05) |> Array.length) (qNew |> Array.filter (fun q -> q <= 0.1) |> Array.length)
                        keepPsmRow
                    )
                run, newIons, pNew, qNew
                )

        // calibration per source intensity quartile (quartiles of the new ions)
        let logSourceQuant (ion: Ion) = log10 ion.Align.Quant_SourceFile
        let newLogSource = scored |> Array.collect (fun (_, ions, _, _) -> ions |> Array.map logSourceQuant) |> Array.filter isFinite |> Array.sort
        if newLogSource.Length >= 4 then
            let edges = [| 0.25; 0.5; 0.75 |] |> Array.map (fun f -> newLogSource.[int (f * float (newLogSource.Length - 1))])
            let quartile v = if not (isFinite v) then -1 else edges |> Array.filter (fun e -> v >= e) |> Array.length
            let newP = scored |> Array.collect (fun (_, ions, pNew, _) -> Array.zip ions pNew)
            for b in 0 .. 3 do
                let ov = oof |> Array.filter (fun o -> quartile (logSourceQuant o.Ion) = b)
                let nw = newP |> Array.filter (fun (ion, _) -> quartile (logSourceQuant ion) = b)
                if ov.Length > 0 then
                    logger.Info (
                        sprintf "Source intensity quartile %i (log10 source quant of the new ions, edges %s): overlap n %i observed positive %.3f mean out-of-fold p %.3f; new n %i mean p %.3f"
                            (b + 1) (edges |> Array.map (sprintf "%.2f") |> String.concat "/") ov.Length
                            (ov |> Array.averageBy (fun o -> if o.Ion.Truth.Value = Positive then 1. else 0.)) (ov |> Array.averageBy (fun o -> o.P))
                            nw.Length (if nw.Length = 0 then nan else nw |> Array.averageBy snd)
                        )

        if testSet then
            let q = Dictionary<IonKey * string,float>(HashIdentity.Structural)
            Array.iter2 (fun o v -> q.[(o.Ion.Key, o.Run)] <- v) oof oofQ
            oof |> Array.groupBy (fun o -> o.Run) |> Array.iter (fun (run, items) -> writeTestSet outputDir fs run items q)

        if diagCharts then
            let histogram (values: float []) (name: string) = values |> Chart.Histogram |> Chart.withTraceInfo name
            [
                for truth in [ Positive; Ambiguous; Negative ] do
                    histogram (oof |> Array.filter (fun o -> o.Ion.Truth.Value = truth) |> Array.map (fun o -> o.P)) (sprintf "Overlap %A (out-of-fold)" truth)
                histogram (scored |> Array.collect (fun (_, _, pNew, _) -> pNew)) "New"
                if useDecoys then histogram pooledDecoys "Decoy (out-of-fold)"
            ]
            |> Chart.combine
            |> Chart.withXAxisStyle "AlignmentScore"
            |> Chart.withYAxisStyle "Count"
            |> Chart.saveHtml (Path.Combine(outputDir, "AlignmentScoreHistogram"))

        // output: the PSM rows, without those the alignment rows replace, and the scored alignment rows of the new ions
        scored
        |> FSharpAux.PSeq.withDegreeOfParallelism c
        |> FSharpAux.PSeq.iter (fun (run, newIons, pNew, qNew) ->
            let replaced = HashSet<IonKey>(newIons |> Seq.filter (fun ion -> ion.Truth.IsSome) |> Seq.map (fun ion -> ion.Key), HashIdentity.Structural)
            let psmRows = run.PsmRows |> Array.filter (fun q -> replaced.Contains (keyOfQuant q) |> not)
            let rows =
                Array.append
                    psmRows
                    (Array.map3 (fun (ion: Ion) score q -> { ion.Row with AlignmentScore = score; AlignmentQValue = q }) newIons pNew qNew)
            let outPath = Path.Combine(outputDir, Path.GetFileName run.QuantPath)
            File.WriteAllLines(outPath, rows |> SeqIO'.csv "\t" true false)
            (Logging.createLogger run.Name).Trace (sprintf "%i PSM and %i alignment rows written to %s; %i PSM rows are replaced by alignment rows" psmRows.Length newIons.Length outPath (run.PsmRows.Length - psmRows.Length))
            )
