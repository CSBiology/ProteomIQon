namespace ProteomIQon

open System
open System.Collections.Generic
open System.Diagnostics
open Plotly.NET

/// m/z correction of a run in ppm, over scan time and m/z.
module MzCalibration =

    open Numerics

    /// A measured precursor m/z of an ion at a scan time.
    type Observation =
        {
            ModSequenceID : int
            Charge        : int
            GlobalMod     : int
            /// Scan time (min).
            ScanTime      : float
            TheoreticalMz : float
            MeasuredMz    : float
        }

    /// An offset the correction is fitted on.
    type Point =
        {
            /// Scan time (min).
            Time      : float
            /// Theoretical m/z.
            Mz        : float
            /// (measured - theoretical) / theoretical x 1e6.
            OffsetPpm : float
        }

    /// m/z correction of the run.
    type Calibration =
        {
            /// Offset (Da, measured - theoretical) at a scan time (min) and a theoretical m/z.
            CorrectionAt         : float -> float -> float
            /// mzWindowSdFactor times the standard deviation (Da) of the precursor m/z offsets, the same for all m/z.
            HalfWindow           : float
            /// Knots of the chosen model on the candidate grid (scan time, m/z); (1, 1) is a constant, (0, 0) no correction.
            Spec                 : int * int
            /// Number of points after merging repeated measurements of one ion.
            PointsUnique         : int
            /// Number of points after removing isotope assignment errors.
            PointsNoIsotopeError : int
            /// Number of points after the fences on the offsets.
            PointsFence1         : int
            /// Number of points after the fences on the residuals; the points of the final fit.
            PointsFence2         : int
            /// The points of the final fit (the points without isotope assignment error if there are too few for a fit).
            Points               : Point []
            /// Time of the fit (ms).
            Milliseconds         : float
        }

    /// Factor of the half windows of XicExtraction.MzWindow_Da Estimate. HalfWindow is this factor times the standard deviation (Da)
    /// of the precursor m/z offsets, residualHalfWindow this factor times the standard deviation (Da) of the residuals around the correction.
    let mzWindowSdFactor = 4.

    /// Window half width (Da) when the standard deviation cannot be estimated.
    let fallbackMzWindow = 0.05

    /// Isotope assignment errors: points with |offset in Da| x charge of at least this many Da are left out (half the 13C-12C spacing).
    let private isotopeLimit = 0.5

    /// Tukey fences of 3 IQR.
    let private fenceFactor = 3.

    /// Fewer points than this after the first fence: the median.
    let private minPointsForModel = 20

    let private folds = 10

    /// Floors of points per node (Cox et al. 2011).
    let private pointsPerRtNode = 50.
    let private pointsPerMzNode = 80.

    /// Candidate numbers of knots: 1, 2, 3, then 2^k and 1.5 x 2^k.
    let private rtGrid = Array.append [| 1; 2; 3 |] [| for k in 2 .. 19 do yield 1 <<< k; yield 3 <<< (k - 1) |]
    let private mzGrid = [| 1; 2; 3; 4; 6; 8; 12; 16 |]

    let private splitmix64 (x: uint64) =
        let mutable z = x + 0x9E3779B97F4A7C15UL
        z <- (z ^^^ (z >>> 30)) * 0xBF58476D1CE4E5B9UL
        z <- (z ^^^ (z >>> 27)) * 0x94D049BB133111EBUL
        z ^^^ (z >>> 31)

    /// Fold of an ion: all points of the ion share it. GlobalMod above 1 collides with the charge in the key, which only groups
    /// the ions more coarsely.
    let private foldOf (modSequenceID: int) (charge: int) (globalMod: int) =
        let key = uint64 modSequenceID * 64UL + uint64 charge * 2UL + uint64 globalMod + (2UL <<< 40)
        int (splitmix64 key % uint64 folds)

    /// Indices of the values inside the Tukey fences (explicit quantiles, linear between the order statistics).
    let private insideFences (values: float []) =
        let sorted = Array.sort values
        let q1 = quantileSorted sorted 0.25
        let q3 = quantileSorted sorted 0.75
        let iqr = q3 - q1
        [| for i in 0 .. values.Length - 1 do
            if values.[i] >= q1 - fenceFactor * iqr && values.[i] <= q3 + fenceFactor * iqr then yield i |]

    /// Points sorted by scan time (and m/z) with their ion fold.
    type private Sample =
        {
            T    : float []
            Z    : float []
            Y    : float []
            Fold : int []
        }

    let private pick (p: Sample) (index: int []) =
        {
            T    = index |> Array.map (fun i -> p.T.[i])
            Z    = index |> Array.map (fun i -> p.Z.[i])
            Y    = index |> Array.map (fun i -> p.Y.[i])
            Fold = index |> Array.map (fun i -> p.Fold.[i])
        }

    /// Knots at the quantiles of the sorted values, equal values merge into one knot. The positions are those of numpy.linspace.
    let private quantileKnots (sorted: float []) (count: int) =
        let step = 1. / float (count - 1)
        let last = float (sorted.Length - 1)
        Array.init count (fun i ->
            let fraction = if count = 1 then 0. elif i = count - 1 then 1. else float i * step
            sorted.[int (floor (fraction * last + 0.5))]
            )
        |> Array.distinct

    /// Left knot of x and the weight of the right knot; the curve is constant beyond the outer knots.
    let private hatAt (knots: float []) (x: float) =
        if knots.Length = 1 then
            struct (0, 0.)
        else
            let j = max 0 (min (SortedSearch.upperBound knots x - 1) (knots.Length - 2))
            struct (j, max 0. (min 1. ((x - knots.[j]) / (knots.[j + 1] - knots.[j]))))

    /// Knots and the hat position of every point.
    type private Hats =
        {
            Knots : float []
            J     : int []
            W     : float []
        }

    let private hatsOf (knots: float []) (x: float []) =
        let j = Array.zeroCreate x.Length
        let w = Array.zeroCreate x.Length
        for i in 0 .. x.Length - 1 do
            let struct (left, weight) = hatAt knots x.[i]
            j.[i] <- left
            w.[i] <- weight
        { Knots = knots; J = j; W = w }

    /// No m/z term: one knot and no column.
    let private noMz = { Knots = [| 0. |]; J = [||]; W = [||] }

    /// The basis of a candidate: hats over scan time, hats over m/z without the first one (so the model is identifiable).
    type private Design =
        {
            RtHats : Hats
            MzHats : Hats
        }

    let private rtColumns (d: Design) = d.RtHats.Knots.Length
    let private mzColumns (d: Design) = d.MzHats.Knots.Length - 1

    /// Writes the non-zero entries of the row of point i (columns, values) and returns their number.
    let private fillRow (d: Design) (i: int) (cols: int []) (vals: float []) =
        let rt = rtColumns d
        let mutable c = 0
        cols.[c] <- d.RtHats.J.[i]
        vals.[c] <- 1. - d.RtHats.W.[i]
        c <- c + 1
        if rt > 1 then
            cols.[c] <- d.RtHats.J.[i] + 1
            vals.[c] <- d.RtHats.W.[i]
            c <- c + 1
        if mzColumns d > 0 then
            let j = d.MzHats.J.[i]
            if j > 0 then
                cols.[c] <- rt + j - 1
                vals.[c] <- 1. - d.MzHats.W.[i]
                c <- c + 1
            cols.[c] <- rt + j
            vals.[c] <- d.MzHats.W.[i]
            c <- c + 1
        c

    /// Normal equations of the additive model. The scan time block is tridiagonal, the cross block and the m/z block are dense
    /// (at most 15 columns). The buffers have a capacity and are reused for the candidates.
    type private Normal(capacityRt: int, capacityMz: int) =
        let diag = Array.zeroCreate<float> capacityRt
        let off = Array.zeroCreate<float> capacityRt
        let cross = Array.zeroCreate<float> (capacityRt * capacityMz)
        let gram = Array.zeroCreate<float> (capacityMz * capacityMz)
        let rhs = Array.zeroCreate<float> (capacityRt + capacityMz)
        let pivot = Array.zeroCreate<float> capacityRt
        let factor = Array.zeroCreate<float> capacityRt
        let work = Array.zeroCreate<float> (capacityRt * (capacityMz + 1))
        let mutable rt = 0
        let mutable mz = 0

        member _.Rt : int = rt
        member _.Mz : int = mz
        member _.Diag : float [] = diag
        member _.Off : float [] = off
        member _.Cross : float [] = cross
        member _.Gram : float [] = gram
        member _.Rhs : float [] = rhs

        member _.Reset (rtCount: int, mzCount: int) =
            rt <- rtCount
            mz <- mzCount
            Array.Clear (diag, 0, rt)
            Array.Clear (off, 0, rt)
            Array.Clear (cross, 0, rt * mz)
            Array.Clear (gram, 0, mz * mz)
            Array.Clear (rhs, 0, rt + mz)

        member _.Add (cols: int [], vals: float [], count: int, y: float) =
            for a in 0 .. count - 1 do
                let ca = cols.[a]
                let va = vals.[a]
                rhs.[ca] <- rhs.[ca] + va * y
                for e in 0 .. count - 1 do
                    let ce = cols.[e]
                    let product = va * vals.[e]
                    if ca < rt && ce < rt then
                        if ca = ce then diag.[ca] <- diag.[ca] + product
                        elif ce = ca + 1 then off.[ca] <- off.[ca] + product
                    elif ca < rt then
                        cross.[ca * mz + ce - rt] <- cross.[ca * mz + ce - rt] + product
                    elif ce >= rt then
                        gram.[(ca - rt) * mz + ce - rt] <- gram.[(ca - rt) * mz + ce - rt] + product

        member _.AddInto (o: Normal) =
            for i in 0 .. rt - 1 do
                diag.[i] <- diag.[i] + o.Diag.[i]
                off.[i] <- off.[i] + o.Off.[i]
            for i in 0 .. rt * mz - 1 do cross.[i] <- cross.[i] + o.Cross.[i]
            for i in 0 .. mz * mz - 1 do gram.[i] <- gram.[i] + o.Gram.[i]
            for i in 0 .. rt + mz - 1 do rhs.[i] <- rhs.[i] + o.Rhs.[i]

        /// This is total - o, the training set of a fold.
        member this.SetDifference (total: Normal, o: Normal) =
            this.Reset (total.Rt, total.Mz)
            for i in 0 .. rt - 1 do
                diag.[i] <- total.Diag.[i] - o.Diag.[i]
                off.[i] <- total.Off.[i] - o.Off.[i]
            for i in 0 .. rt * mz - 1 do cross.[i] <- total.Cross.[i] - o.Cross.[i]
            for i in 0 .. mz * mz - 1 do gram.[i] <- total.Gram.[i] - o.Gram.[i]
            for i in 0 .. rt + mz - 1 do rhs.[i] <- total.Rhs.[i] - o.Rhs.[i]

        /// Least squares solution with a ridge of 1e-10 x the mean diagonal: LDL' of the tridiagonal block, Schur complement on the
        /// m/z block, Cholesky of that. No solution if a pivot is not positive.
        member _.Solve () : float [] voption =
            let mutable sumDiagonal = 0.
            for i in 0 .. rt - 1 do sumDiagonal <- sumDiagonal + diag.[i]
            for a in 0 .. mz - 1 do sumDiagonal <- sumDiagonal + gram.[a * mz + a]
            let ridge = 1e-10 * max (sumDiagonal / float (rt + mz)) 1e-30
            let mutable valid = true
            pivot.[0] <- diag.[0] + ridge
            for i in 1 .. rt - 1 do
                factor.[i - 1] <- off.[i - 1] / pivot.[i - 1]
                pivot.[i] <- diag.[i] + ridge - factor.[i - 1] * off.[i - 1]
            for i in 0 .. rt - 1 do
                if not (pivot.[i] > 0.) then valid <- false
            if not valid then
                ValueNone
            else
            // work holds the solution of the scan time block for the right side (column 0) and for each m/z column of the cross block
            let solveTridiagonal (column: int) =
                let start = column * rt
                for i in 1 .. rt - 1 do work.[start + i] <- work.[start + i] - factor.[i - 1] * work.[start + i - 1]
                for i in 0 .. rt - 1 do work.[start + i] <- work.[start + i] / pivot.[i]
                for i in rt - 2 .. -1 .. 0 do work.[start + i] <- work.[start + i] - factor.[i] * work.[start + i + 1]
            Array.blit rhs 0 work 0 rt
            solveTridiagonal 0
            if mz = 0 then
                ValueSome (Array.sub work 0 rt)
            else
            for j in 0 .. mz - 1 do
                for i in 0 .. rt - 1 do work.[(j + 1) * rt + i] <- cross.[i * mz + j]
                solveTridiagonal (j + 1)
            let schur = Array2D.zeroCreate<float> mz mz
            let right = Array.zeroCreate<float> mz
            for a in 0 .. mz - 1 do
                for e in 0 .. mz - 1 do
                    let mutable s = gram.[a * mz + e] + (if a = e then ridge else 0.)
                    for i in 0 .. rt - 1 do s <- s - cross.[i * mz + a] * work.[(e + 1) * rt + i]
                    schur.[a, e] <- s
                let mutable s = rhs.[rt + a]
                for i in 0 .. rt - 1 do s <- s - cross.[i * mz + a] * work.[i]
                right.[a] <- s
            // Cholesky of the Schur complement
            let l = Array2D.zeroCreate<float> mz mz
            for j in 0 .. mz - 1 do
                let mutable s = schur.[j, j]
                for k in 0 .. j - 1 do s <- s - l.[j, k] * l.[j, k]
                if not (s > 0.) then valid <- false
                let d = sqrt s
                l.[j, j] <- d
                for i in j + 1 .. mz - 1 do
                    let mutable s = schur.[i, j]
                    for k in 0 .. j - 1 do s <- s - l.[i, k] * l.[j, k]
                    l.[i, j] <- s / d
            if not valid then
                ValueNone
            else
            let y = Array.zeroCreate<float> mz
            for i in 0 .. mz - 1 do
                let mutable s = right.[i]
                for k in 0 .. i - 1 do s <- s - l.[i, k] * y.[k]
                y.[i] <- s / l.[i, i]
            let xm = Array.zeroCreate<float> mz
            for i in mz - 1 .. -1 .. 0 do
                let mutable s = y.[i]
                for k in i + 1 .. mz - 1 do s <- s - l.[k, i] * xm.[k]
                xm.[i] <- s / l.[i, i]
            let beta = Array.zeroCreate<float> (rt + mz)
            for i in 0 .. rt - 1 do
                let mutable s = work.[i]
                for j in 0 .. mz - 1 do s <- s - work.[(j + 1) * rt + i] * xm.[j]
                beta.[i] <- s
            Array.blit xm 0 beta rt mz
            ValueSome beta

    /// Grouped cross validation of a design on the points: the mean squared out-of-fold residual, its standard error
    /// (sd of the fold means / sqrt of the folds) and the number of coefficients. Folds with fewer than 2 points are not evaluated. No result if fewer
    /// than 2 folds are evaluated, a fit fails or the error is not finite.
    let private crossValidate (blocks: Normal []) (total: Normal) (scratch: Normal) (design: Design) (y: float []) (fold: int []) =
        let rt = rtColumns design
        let mz = mzColumns design
        let cols = Array.zeroCreate<int> 4
        let vals = Array.zeroCreate<float> 4
        for block in blocks do block.Reset (rt, mz)
        let counts = Array.zeroCreate<int> folds
        for i in 0 .. y.Length - 1 do
            let count = fillRow design i cols vals
            blocks.[fold.[i]].Add (cols, vals, count, y.[i])
            counts.[fold.[i]] <- counts.[fold.[i]] + 1
        total.Reset (rt, mz)
        for block in blocks do total.AddInto block
        let betas = Array.create folds (ValueNone: float [] voption)
        for k in 0 .. folds - 1 do
            if counts.[k] >= 2 then
                scratch.SetDifference (total, blocks.[k])
                betas.[k] <- scratch.Solve ()
        let evaluated = [| for k in 0 .. folds - 1 do if counts.[k] >= 2 then yield k |]
        if evaluated.Length < 2 || evaluated |> Array.exists (fun k -> betas.[k].IsNone) then
            ValueNone
        else
            let squares = Array.zeroCreate<float> folds
            for i in 0 .. y.Length - 1 do
                match betas.[fold.[i]] with
                | ValueSome beta ->
                    let count = fillRow design i cols vals
                    let mutable prediction = 0.
                    for a in 0 .. count - 1 do prediction <- prediction + vals.[a] * beta.[cols.[a]]
                    let r = y.[i] - prediction
                    squares.[fold.[i]] <- squares.[fold.[i]] + r * r
                | ValueNone -> ()
            let foldMeans = evaluated |> Array.map (fun k -> squares.[k] / float counts.[k])
            let mean = Array.average foldMeans
            let sd = sqrt ((foldMeans |> Array.sumBy (fun v -> (v - mean) ** 2.)) / float (foldMeans.Length - 1))
            let mse = (evaluated |> Array.sumBy (fun k -> squares.[k])) / float (evaluated |> Array.sumBy (fun k -> counts.[k]))
            let se = sd / sqrt (float foldMeans.Length)
            if isFinite mse && isFinite se then ValueSome (mse, se, rtColumns design + mzColumns design) else ValueNone

    let private newBlocks (capacityRt: int) (capacityMz: int) =
        Array.init folds (fun _ -> Normal(capacityRt, capacityMz)), Normal(capacityRt, capacityMz), Normal(capacityRt, capacityMz)

    /// The design of a candidate from the points (sorted by scan time); the knots come from the x values of the points only.
    let private designOf (p: Sample) (rtKnots: int) (mzKnots: int) =
        {
            RtHats = hatsOf (quantileKnots p.T rtKnots) p.T
            MzHats = if mzKnots > 1 then hatsOf (quantileKnots (Array.sort p.Z) mzKnots) p.Z else noMz
        }

    /// A fitted model.
    type private Model =
        {
            RtKnots : float []
            MzKnots : float []
            Beta    : float []
        }

    /// Offset in ppm of the model at a scan time and an m/z.
    let private predict (m: Model) (t: float) (z: float) =
        let rt = m.RtKnots.Length
        let struct (jt, wt) = hatAt m.RtKnots t
        let mutable f = m.Beta.[jt] * (1. - wt)
        if rt > 1 then f <- f + m.Beta.[jt + 1] * wt
        if m.MzKnots.Length > 1 then
            let struct (jz, wz) = hatAt m.MzKnots z
            if jz > 0 then f <- f + m.Beta.[rt + jz - 1] * (1. - wz)
            f <- f + m.Beta.[rt + jz] * wz
        f

    /// Fits the candidate on the points; knots are rebuilt from them. None if the fit fails.
    let private fitModel (rtKnots: int, mzKnots: int) (p: Sample) =
        let design = designOf p rtKnots mzKnots
        let rt = rtColumns design
        let mz = mzColumns design
        let normal = Normal(rt, mz)
        normal.Reset (rt, mz)
        let cols = Array.zeroCreate<int> 4
        let vals = Array.zeroCreate<float> 4
        for i in 0 .. p.Y.Length - 1 do
            let count = fillRow design i cols vals
            normal.Add (cols, vals, count, p.Y.[i])
        match normal.Solve () with
        | ValueSome beta when beta |> Array.forall isFinite ->
            Some { RtKnots = design.RtHats.Knots; MzKnots = design.MzHats.Knots; Beta = beta }
        | _ -> None

    /// The candidate with the fewest coefficients among those whose cross validation error is within one standard error of the
    /// lowest (one-standard-error rule); the lower error and then the earlier candidate in the grid break ties.
    /// None if no candidate can be evaluated.
    let private selectSpec (p: Sample) =
        let n = p.Y.Length
        let allowedRt = rtGrid |> Array.filter (fun a -> a = 1 || float n / float a >= pointsPerRtNode)
        let allowedMz = mzGrid |> Array.filter (fun b -> b = 1 || float n / float b >= pointsPerMzNode)
        let sortedZ = Array.sort p.Z
        // the hats over m/z are the same for all scan time candidates
        let mzHats = allowedMz |> Array.map (fun b -> if b > 1 then hatsOf (quantileKnots sortedZ b) p.Z else noMz)
        let blocks, total, scratch = newBlocks (Array.max allowedRt) (Array.max allowedMz - 1)
        let candidates = ResizeArray<(int * int) * int * float * float>()
        for a in allowedRt do
            let rt = hatsOf (quantileKnots p.T a) p.T
            allowedMz
            |> Array.iteri (fun bi b ->
                match crossValidate blocks total scratch { RtHats = rt; MzHats = mzHats.[bi] } p.Y p.Fold with
                | ValueSome (mse, se, coefficients) -> candidates.Add (((a, b), coefficients, mse, se))
                | ValueNone -> ()
                )
        if candidates.Count = 0 then
            None
        else
            let _, _, bestMse, bestSe = candidates |> Seq.minBy (fun (_, _, mse, _) -> mse)
            let spec, _, _, _ =
                candidates
                |> Seq.filter (fun (_, _, mse, _) -> mse <= bestMse + bestSe)
                |> Seq.minBy (fun (_, coefficients, mse, _) -> coefficients, mse)
            Some spec

    let private compareObservations (a: Observation) (b: Observation) =
        let c = compare a.ScanTime b.ScanTime
        if c <> 0 then c else
        let c = compare a.TheoreticalMz b.TheoreticalMz
        if c <> 0 then c else
        let c = compare a.ModSequenceID b.ModSequenceID
        if c <> 0 then c else
        let c = compare a.Charge b.Charge
        if c <> 0 then c else
        let c = compare a.GlobalMod b.GlobalMod
        if c <> 0 then c else
        compare a.MeasuredMz b.MeasuredMz

    let private isValid (o: Observation) =
        isFinite o.ScanTime && isFinite o.TheoreticalMz && isFinite o.MeasuredMz && o.TheoreticalMz > 0. && o.Charge > 0 && o.GlobalMod >= 0 && o.ModSequenceID >= 0

    let private ppm (theoretical: float) (measured: float) = (measured - theoretical) / theoretical * 1e6

    /// Half window (Da): mzWindowSdFactor times the standard deviation of the offsets (Da, measured - theoretical) inside the Tukey
    /// fences of 3 IQR and below 0.1 Da. The fallback window if the standard deviation cannot be estimated.
    let private halfWindowOf (observations: Observation []) =
        let diffs = observations |> Array.map (fun o -> o.MeasuredMz - o.TheoreticalMz)
        let lowerFence, upperFence =
            match (diffs |> FSharp.Stats.Signal.Outliers.tukey fenceFactor).TryToTuple with
            | Some fences -> fences
            | None -> nan, nan
        let filtered = diffs |> Array.filter (fun d -> d <= upperFence && d >= lowerFence && d < 0.1)
        let sd =
            if filtered.Length < 2 then
                nan
            else
                let mean = Array.average filtered
                sqrt ((filtered |> Array.sumBy (fun v -> (v - mean) ** 2.)) / float (filtered.Length - 1))
        if isFinite sd && sd > 0. then mzWindowSdFactor * sd else fallbackMzWindow

    /// Calibration from a correction in ppm: the Da offset is the ppm offset at the theoretical m/z.
    let private calibrationOf (correctionPpm: float -> float -> float) halfWindow spec counts (points: Point []) (stopwatch: Stopwatch) =
        let unique, noIsotopeError, fence1, fence2 = counts
        let correctionAt (t: float) (mz: float) =
            if isFinite t && isFinite mz then
                let c = correctionPpm t mz * mz * 1e-6
                if isFinite c then c else 0.
            else
                0.
        {
            CorrectionAt         = correctionAt
            HalfWindow           = halfWindow
            Spec                 = spec
            PointsUnique         = unique
            PointsNoIsotopeError = noIsotopeError
            PointsFence1         = fence1
            PointsFence2         = fence2
            Points               = points
            Milliseconds         = stopwatch.Elapsed.TotalMilliseconds
        }

    /// Fits the offset in ppm as f(scan time) + g(m/z), both piecewise linear on quantile knots, whose numbers of knots are chosen by
    /// grouped cross validation (one-standard-error rule). Steps:
    /// repeated measurements of an ion at the same m/z become one point at their mean scan time; isotope assignment errors
    /// (|offset in Da| x charge of 0.5 or more) are left out; Tukey fences on the offsets; model selection; fit; Tukey fences on the
    /// residuals of all points without isotope assignment error; refit of the chosen model. The window does not depend on the fit,
    /// it comes from the offsets of all observations. Fewer than 20 points after the first fence, or a failing fit: the median.
    /// Fewer than 3 points without isotope assignment error: no correction.
    let fit (observations: Observation []) =
        let stopwatch = Stopwatch.StartNew ()
        let valid = observations |> Array.filter isValid
        let halfWindow = halfWindowOf valid
        let sorted = valid |> Array.sortWith compareObservations
        // one point per ion and measured m/z at the mean scan time, in the order of the scan times
        let unique =
            let groups = Dictionary<struct (int * int * int * float), int>()
            let first = ResizeArray<Observation>()
            let sums = ResizeArray<float>()
            let counts = ResizeArray<int>()
            for o in sorted do
                match groups.TryGetValue (struct (o.ModSequenceID, o.Charge, o.GlobalMod, o.MeasuredMz)) with
                | true, g ->
                    sums.[g] <- sums.[g] + o.ScanTime
                    counts.[g] <- counts.[g] + 1
                | _ ->
                    groups.[struct (o.ModSequenceID, o.Charge, o.GlobalMod, o.MeasuredMz)] <- first.Count
                    first.Add o
                    sums.Add o.ScanTime
                    counts.Add 1
            first.ToArray ()
            |> Array.mapi (fun g o -> { o with ScanTime = sums.[g] / float counts.[g] })
            |> Array.sortWith compareObservations
        let withoutIsotopeError = unique |> Array.filter (fun o -> abs ((o.MeasuredMz - o.TheoreticalMz) * float o.Charge) < isotopeLimit)
        if withoutIsotopeError.Length < 3 then
            let points = withoutIsotopeError |> Array.map (fun o -> { Time = o.ScanTime; Mz = o.TheoreticalMz; OffsetPpm = ppm o.TheoreticalMz o.MeasuredMz })
            calibrationOf (fun _ _ -> 0.) halfWindow (0, 0) (unique.Length, withoutIsotopeError.Length, 0, 0) points stopwatch
        else
        let all =
            {
                T    = withoutIsotopeError |> Array.map (fun o -> o.ScanTime)
                Z    = withoutIsotopeError |> Array.map (fun o -> o.TheoreticalMz)
                Y    = withoutIsotopeError |> Array.map (fun o -> ppm o.TheoreticalMz o.MeasuredMz)
                Fold = withoutIsotopeError |> Array.map (fun o -> foldOf o.ModSequenceID o.Charge o.GlobalMod)
            }
        let selection = pick all (insideFences all.Y)
        let counts fence2 = unique.Length, withoutIsotopeError.Length, selection.Y.Length, fence2
        let pointsOf (p: Sample) = Array.init p.Y.Length (fun i -> { Time = p.T.[i]; Mz = p.Z.[i]; OffsetPpm = p.Y.[i] })
        // the median of the points inside the first fence, fences on the residuals around it
        let medianFit () =
            let m = median selection.Y
            let final = pick all (insideFences (all.Y |> Array.map (fun y -> y - m)))
            calibrationOf (fun _ _ -> m) halfWindow (1, 1) (counts final.Y.Length) (pointsOf final) stopwatch
        if selection.Y.Length < minPointsForModel then
            medianFit ()
        else
        match selectSpec selection |> Option.bind (fun spec -> fitModel spec selection |> Option.map (fun model -> spec, model)) with
        | None -> medianFit ()
        | Some (spec, firstModel) ->
            let residuals = Array.init all.Y.Length (fun i -> all.Y.[i] - predict firstModel all.T.[i] all.Z.[i])
            let final = pick all (insideFences residuals)
            match fitModel spec final with
            | None -> medianFit ()
            | Some model ->
                calibrationOf (predict model) halfWindow spec (counts final.Y.Length) (pointsOf final) stopwatch

    /// The XIC half window (Da): the fixed window, or the estimate of the calibration.
    let mzWindow (window: Domain.Window) (calibration: Calibration) =
        match window with
        | Domain.Window.Fixed v -> v
        | Domain.Window.Estimate -> calibration.HalfWindow

    /// Half width (Da): mzWindowSdFactor times the standard deviation of the differences in Da between the points of the final fit and the correction.
    /// calibration.HalfWindow if it cannot be estimated.
    let residualHalfWindow (calibration: Calibration) =
        let residuals = calibration.Points |> Array.map (fun p -> p.OffsetPpm * p.Mz * 1e-6 - calibration.CorrectionAt p.Time p.Mz)
        let sd = if residuals.Length < 3 then nan else FSharp.Stats.Seq.stDev residuals
        if isFinite sd && sd > 0. then mzWindowSdFactor * sd else calibration.HalfWindow

    /// One line for the log: the XIC half window in Da and where it comes from.
    let describeWindow (window: Domain.Window) (halfWindow: float) (isFallback: bool) =
        match window with
        | Domain.Window.Fixed _ -> sprintf "m/z window +-%.4f Da (Fixed)" halfWindow
        | Domain.Window.Estimate when isFallback -> sprintf "m/z window +-%.4f Da (fallback)" halfWindow
        | Domain.Window.Estimate -> sprintf "m/z window +-%.4f Da (estimate)" halfWindow

    /// The offsets of the calibration points in ppm over scan time with the correction at the same points.
    let calibrationChart (calibration: Calibration) =
        let time = calibration.Points |> Array.map (fun p -> p.Time)
        [
            Chart.Point(time, calibration.Points |> Array.map (fun p -> p.OffsetPpm))
            |> Chart.withTraceInfo "Raw"
            Chart.Point(time, calibration.Points |> Array.map (fun p -> calibration.CorrectionAt p.Time p.Mz / p.Mz * 1e6))
            |> Chart.withTraceInfo "Correction"
        ]
        |> Chart.combine

    /// One line for the log: model, point counts and time.
    let describe (calibration: Calibration) =
        let rt, mz = calibration.Spec
        sprintf "model %i scan time knots x %i m/z knots, points %i unique, %i without isotope error, %i inside the fences of the offsets, %i inside the fences of the residuals, %.0f ms"
            rt mz calibration.PointsUnique calibration.PointsNoIsotopeError calibration.PointsFence1 calibration.PointsFence2 calibration.Milliseconds
