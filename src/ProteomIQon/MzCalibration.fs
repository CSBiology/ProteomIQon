namespace ProteomIQon

open System

/// m/z correction over scan time.
module MzCalibration =

    open Numerics

    /// Running median of y over the k neighbouring points of the points sorted by x, linearly interpolated in between and constant beyond the ends.
    let runningMedianCurve (k: int) (x: float []) (y: float []) =
        let n = x.Length
        let k = max 1 (min k n)
        let fitted =
            Array.init n (fun i ->
                let lo = max 0 (min (i - k / 2) (n - k))
                y.[lo .. lo + k - 1] |> median
                )
        let knots =
            Array.zip x fitted
            |> Array.groupBy fst
            |> Array.map (fun (xv, group) -> xv, group |> Array.averageBy snd)
        let kx = Array.map fst knots
        let ky = Array.map snd knots
        fun (v: float) ->
            let m = kx.Length
            if m = 0 || Double.IsNaN v then nan
            elif v <= kx.[0] then ky.[0]
            elif v >= kx.[m - 1] then ky.[m - 1]
            else
                let hi = SortedSearch.upperBound kx v
                let lo = hi - 1
                let t = (v - kx.[lo]) / (kx.[hi] - kx.[lo])
                ky.[lo] + t * (ky.[hi] - ky.[lo])

    /// m/z correction of the run.
    type MzCalibration =
        {
            /// Scan time -> m/z offset (measured - theoretical).
            Correction : float -> float
            /// Standard deviation of the offsets (inside the fences) around the correction; nan if it cannot be estimated.
            Sd         : float
            /// (scan time, offset) of the ions used.
            Points     : (float*float) []
        }

    /// Fits the m/z offset over scan time: offsets outside the Tukey fences (3 IQR) or above 0.1 Da are left out, as in the PSM based
    /// quantification; the correction is a running median over scan time for runs longer than
    /// 20 minutes with more than 500 points, else the median.
    let fitMzCalibration (points: (float*float) []) =
        let points = points |> Array.filter (fun (t, d) -> isFinite t && isFinite d) |> Array.sortBy fst
        if points.Length < 3 then
            { Correction = (fun _ -> 0.); Sd = nan; Points = points }
        else
        let lower, upper =
            match (points |> Array.map snd |> FSharp.Stats.Signal.Outliers.tukey 3.).TryToTuple with
            | Some borders -> borders
            | None -> nan, nan
        let filtered =
            if Double.IsNaN upper || Double.IsNaN lower then points
            else points |> Array.filter (fun (_, d) -> d <= upper && d >= lower && d < 0.1)
        let runTime = points |> Array.last |> fst
        let correction =
            if runTime > 20. && filtered.Length > 500 then
                runningMedianCurve 101 (Array.map fst filtered) (Array.map snd filtered)
            else
                let m = filtered |> Array.map snd |> median
                fun _ -> m
        let residuals = filtered |> Array.map (fun (t, d) -> d - correction t)
        let sd = if residuals.Length < 3 then nan else FSharp.Stats.Seq.stDev residuals
        { Correction = correction; Sd = sd; Points = points }
