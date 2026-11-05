namespace ProteomIQon

open System
open System.Collections.Generic
open BioFSharp.Mz
open MzIO.Commons.Arrays
open MzIO.Processing
open MzIO.Processing.MzIOLinq

/// Processing of the raw XIC of a precursor.
module Xic =

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
    let extractPoint (spectrum: SortedPeaks.MzSortedPeaks) (lockMz: float) (halfWindow: float) (imLow: float) (imHigh: float) =
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

    /// Reads the peaks of a spectrum sorted by m/z. The ion mobilities are empty when the spectrum has none.
    let readSortedPeaks (reader: MzIO.MzSQL.MzSQL) (spectrumID: string) =
        let struct (mzs, intensities, ionMobilities) = reader.ReadSpectrumPeaksAsArrays spectrumID
        SortedPeaks.sortPeakArraysByMz mzs intensities (if isNull ionMobilities then [||] else ionMobilities)

    /// The peaks of a spectrum with lowMz < m/z < highMz, in file order.
    let getPeaksInMzRange (peaks: SortedPeaks.MzSortedPeaks) lowMz highMz =
        let n = peaks.Mz.Length
        // first peak with Mz > lowMz
        let lo = SortedSearch.upperBound peaks.Mz lowMz
        let mutable last = lo
        while last < n && peaks.Mz.[last] < highMz do
            last <- last + 1
        Array.init (last - lo) (fun k -> lo + k)
        |> Array.sortBy (fun k -> peaks.Index.[k])
        |> Array.map (fun k -> Peak(peaks.Mz.[k], peaks.Intensity.[k]))

    /// Cache of the sorted peaks of MS1 spectra that keeps at most maxEntries spectra. When it is full, the tenth of
    /// the spectra that was used least recently is dropped.
    type PeakCache(read: string -> SortedPeaks.MzSortedPeaks, maxEntries: int) =
        let cache = Dictionary<string, struct (SortedPeaks.MzSortedPeaks * int64 ref)>()
        let mutable clock = 0L
        let mutable reads = 0
        /// The number of spectra that were read, including the ones that were read again after they were dropped.
        member _.Reads = reads
        member _.Capacity = maxEntries
        member _.Get(spectrumID: string) =
            clock <- clock + 1L
            match cache.TryGetValue spectrumID with
            | true, struct (peaks, lastUse) ->
                lastUse.Value <- clock
                peaks
            | _ ->
                let peaks = read spectrumID
                reads <- reads + 1
                if cache.Count >= maxEntries then
                    cache
                    |> Seq.map (fun kv -> let struct (_, lastUse) = kv.Value in kv.Key, lastUse.Value)
                    |> Seq.sortBy snd
                    |> Seq.truncate (max 1 (maxEntries / 10))
                    |> Seq.toArray
                    |> Array.iter (fun (key, _) -> cache.Remove key |> ignore)
                cache.[spectrumID] <- struct (peaks, ref clock)
                peaks

    /// The size of a PeakCache for ions quantified in the order of their scan time, whose XICs read the MS1 spectra of a window of
    /// 2 x scanTimeWindow. Returns the cache size (twice the most MS1 spectra in such a window, at least 16) and that most number of spectra.
    let peakCacheSize (ms1ScanTimes: float []) (scanTimeWindow: float) =
        let windowWidth = 2. * scanTimeWindow
        let mutable windowEnd = 0
        let mutable widestWindow = 1
        for windowStart in 0 .. ms1ScanTimes.Length - 1 do
            while windowEnd < ms1ScanTimes.Length && ms1ScanTimes.[windowEnd] < ms1ScanTimes.[windowStart] + windowWidth do
                windowEnd <- windowEnd + 1
            widestWindow <- max widestWindow (windowEnd - windowStart)
        max 16 (2 * widestWindow), widestWindow

    /// Maps f over xs in the order of key (ascending) and returns the results in the order of xs. f gets the position of the
    /// element in the order of key.
    let mapiOrderedBy (key: 'a -> float) (f: int -> 'a -> 'b) (xs: 'a []) =
        let order = Array.init xs.Length id |> Array.sortBy (fun i -> key xs.[i], i)
        let results = Array.zeroCreate xs.Length
        order |> Array.iteri (fun position i -> results.[i] <- f position xs.[i])
        results

    /// The XIC of an ion: the scan times of the MS1 spectra inside rtRange and the intensity of each (Xic.extractPoint at lockMz
    /// +- mzHalfWindow). With an ion mobility window (imCenter +- imHalfWindow) the peaks inside it are summed per m/z; imHalfWindow
    /// nan means no ion mobility window, imCenter nan with a window means no peak is inside it.
    let rtProfile (readSortedPeaks: string -> SortedPeaks.MzSortedPeaks) (rtIndex: IMzIOArray<RtIndexEntry>) (rtRange: RangeQuery) (lockMz: float) (mzHalfWindow: float) (imCenter: float) (imHalfWindow: float) =
        RtIndexEntry.Search(rtIndex, rtRange)
        |> Seq.map (fun entry ->
            let intensity =
                if Double.IsNaN imHalfWindow then
                    let intensity, _, _ = extractPoint (readSortedPeaks entry.SpectrumID) lockMz mzHalfWindow nan nan
                    intensity
                elif Double.IsNaN imCenter then
                    0.
                else
                    let intensity, _, _ = extractPoint (readSortedPeaks entry.SpectrumID) lockMz mzHalfWindow (imCenter - imHalfWindow) (imCenter + imHalfWindow)
                    intensity
            entry.Rt, intensity
            )
        |> Array.ofSeq
        |> Array.unzip
