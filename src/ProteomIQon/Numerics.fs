namespace ProteomIQon

open System

/// Small numeric helpers on float arrays.
module Numerics =

    let isFinite (v: float) = not (Double.IsNaN v || Double.IsInfinity v)

    /// Median of the values, nan for an empty array.
    let median (values: float []) =
        if values.Length = 0 then
            nan
        else
            let sorted = Array.sort values
            let n = sorted.Length
            if n % 2 = 1 then sorted.[n / 2] else 0.5 * (sorted.[n / 2 - 1] + sorted.[n / 2])

    /// Pearson correlation, nan for fewer than 4 points or a constant input.
    let pearson (a: float []) (b: float []) =
        let n = a.Length
        if n < 4 || b.Length <> n then
            nan
        else
            let ma = Array.average a
            let mb = Array.average b
            let mutable sab = 0.
            let mutable saa = 0.
            let mutable sbb = 0.
            for i in 0 .. n - 1 do
                let da = a.[i] - ma
                let db = b.[i] - mb
                sab <- sab + da * db
                saa <- saa + da * da
                sbb <- sbb + db * db
            if saa <= 0. || sbb <= 0. then nan else sab / sqrt (saa * sbb)

    /// Linear interpolation in (xs, ys) sorted by xs, 0 outside.
    let interpolateOrZero (xs: float []) (ys: float []) (x: float) =
        let n = xs.Length
        if n = 0 || x < xs.[0] || x > xs.[n - 1] then
            0.
        else
            let hi = SortedSearch.upperBound xs x
            if hi >= n then ys.[n - 1]
            else
                let lo = hi - 1
                if lo < 0 then ys.[0]
                else
                    let t = (x - xs.[lo]) / (xs.[hi] - xs.[lo])
                    ys.[lo] + t * (ys.[hi] - ys.[lo])

    /// Quantile of a sorted array, linear between neighbouring values (numpy's default); nan for an empty array.
    let quantileSorted (sorted: float []) (q: float) =
        if sorted.Length = 0 then
            nan
        else
            let h = q * float (sorted.Length - 1)
            let lo = int (floor h)
            let hi = min (lo + 1) (sorted.Length - 1)
            sorted.[lo] + (h - float lo) * (sorted.[hi] - sorted.[lo])

    /// Median distance between neighbouring values of a sorted array, nan for fewer than two values.
    let medianSpacing (sorted: float []) =
        if sorted.Length < 2 then nan else Array.init (sorted.Length - 1) (fun i -> sorted.[i + 1] - sorted.[i]) |> median
