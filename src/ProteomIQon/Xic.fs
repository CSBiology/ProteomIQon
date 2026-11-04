namespace ProteomIQon

open System

/// Processing of the raw XIC of a precursor.
module Xic =

    /// Half width of the ion mobility window of an XIC (1/K0) when it is not estimated.
    let ionMobilityHalfWindow = 0.05

    /// A zero point that follows a point more than 100 times higher is a drop-out; the mask is false there.
    let keepMask (y: float []) =
        Array.init y.Length (fun i ->
            if i = 0 || i = y.Length - 1 || y.[i] > 0. then
                true
            else
                let previous = y.[i - 1]
                previous = 0. || previous <= 100. * (y.[i] + 1.)
            )

    /// ALS baseline subtraction for XICs of 3 to 500 points; negative values are set to 0.
    let subtractBaseline (baseline: Domain.BaseLineCorrection option) (y: float []) =
        match baseline with
        | Some p when y.Length >= 3 && y.Length <= 500 ->
            let b = BaselineAls.baselineAls p.MaxIterations p.Lambda p.P y
            Array.map2 (fun v bv -> max 0. (v - bv)) y b
        | _ -> Array.copy y

    /// The XIC points of a trace without drop-outs and the baseline corrected intensities: scan times, corrected and raw intensities.
    let processXic (baseline: Domain.BaseLineCorrection option) (rt: float []) (intensity: float []) =
        let keep = keepMask intensity
        let pick (a: float []) = [| for i in 0 .. a.Length - 1 do if keep.[i] then yield a.[i] |]
        let intensity' = pick intensity
        pick rt, subtractBaseline baseline intensity', intensity'

    /// One XIC point: the peak closest to lockMz within +-halfWindow; intensity 0 if there is none.
    /// With ion mobility the peaks inside the ion mobility window (no window if imLow is nan) are first summed per m/z,
    /// because one m/z is spread over the mobility scans.
    /// Returns intensity, observed m/z and the intensity weighted ion mobility of the chosen m/z (nan without ion mobility).
    let extractPoint (spectrum: TimsPeaks.MzSortedPeaks) (lockMz: float) (halfWindow: float) (imLow: float) (imHigh: float) =
        let mz = spectrum.Mz
        let intensity = spectrum.Intensity
        let high = lockMz + halfWindow
        let mutable i = SortedSearch.lowerBound mz (lockMz - halfWindow)
        let mutable bestDistance = infinity
        let mutable bestIntensity = 0.
        let mutable bestMz = nan
        let mutable bestIonMobility = nan
        if spectrum.IonMobility.Length = 0 then
            while i < mz.Length && mz.[i] <= high do
                let distance = abs (mz.[i] - lockMz)
                if distance < bestDistance then
                    bestDistance <- distance
                    bestIntensity <- intensity.[i]
                    bestMz <- mz.[i]
                i <- i + 1
        else
            let im = spectrum.IonMobility
            let filterIonMobility = not (Double.IsNaN imLow)
            while i < mz.Length && mz.[i] <= high do
                let groupMz = mz.[i]
                let mutable n = 0
                let mutable sum = 0.
                let mutable sumIm = 0.
                while i < mz.Length && mz.[i] = groupMz do
                    if not filterIonMobility || (im.[i] > imLow && im.[i] < imHigh) then
                        n <- n + 1
                        sum <- sum + intensity.[i]
                        sumIm <- sumIm + intensity.[i] * im.[i]
                    i <- i + 1
                if n > 0 then
                    let distance = abs (groupMz - lockMz)
                    if distance < bestDistance then
                        bestDistance <- distance
                        bestIntensity <- sum
                        bestMz <- groupMz
                        bestIonMobility <- if sum > 0. then sumIm / sum else nan
        bestIntensity, bestMz, bestIonMobility
