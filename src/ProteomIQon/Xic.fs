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

    /// The peak closest to lockMz within +-halfWindow; intensity 0 if there is none. It is the XIC point of a spectrum without ion
    /// mobility and the measurement of the m/z calibration of AlignmentBasedQuantification.
    /// With ion mobility the peaks inside the ion mobility window (no window if imLow is nan) are first summed per m/z,
    /// because one m/z is spread over the mobility scans, and the sum closest to lockMz is taken.
    /// Returns intensity, observed m/z and the intensity weighted ion mobility of the chosen m/z (nan without ion mobility).
    let closestPoint (spectrum: SortedPeaks.MzSortedPeaks) (lockMz: float) (halfWindow: float) (imLow: float) (imHigh: float) =
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

    /// The (ion mobility, intensity) pairs of the peaks of a series of XIC points that lie inside the m/z window of the point and
    /// inside a recorded ion mobility range; peaks without ion mobility are not recorded. The pairs of all points are stored one
    /// after the other. It is filled by Xic.extractPointWithProfile, one point per call, and emptied with Clear to be filled again.
    type MobilityProfile() =
        let mobility = ResizeArray<float32>()
        let intensity = ResizeArray<float32>()
        let starts = ResizeArray<int>()
        /// The number of points.
        member _.Points = starts.Count
        /// Starts a point: the pairs added next belong to it.
        member _.StartPoint () = starts.Add mobility.Count
        /// Adds a pair to the point that was started last.
        member _.Add (im: float, y: float) =
            mobility.Add (float32 im)
            intensity.Add (float32 y)
        /// The pairs of a point are the positions from (inclusive) to (exclusive).
        member _.Range (point: int) =
            starts.[point], (if point + 1 < starts.Count then starts.[point + 1] else mobility.Count)
        /// The ion mobility of the pair at a position.
        member _.Mobility (pair: int) = float mobility.[pair]
        /// The intensity of the pair at a position.
        member _.Intensity (pair: int) = float intensity.[pair]
        /// Removes all points; the buffers keep their capacity.
        member _.Clear () =
            mobility.Clear()
            intensity.Clear()
            starts.Clear()

    /// The loop of Xic.extractPoint. record is called for every peak inside the m/z window with its ion mobility and intensity.
    let inline private extractPointWith
            ([<InlineIfLambda>] record: float -> float -> unit)
            (spectrum: SortedPeaks.MzSortedPeaks) (lockMz: float) (halfWindow: float) (imLow: float) (imHigh: float) =
        if spectrum.IonMobility.Length = 0 then
            closestPoint spectrum lockMz halfWindow imLow imHigh
        else
            let mz = spectrum.Mz
            let intensity = spectrum.Intensity
            let im = spectrum.IonMobility
            let high = lockMz + halfWindow
            let filterIonMobility = not (Double.IsNaN imLow)
            let mutable i = SortedSearch.lowerBound mz (lockMz - halfWindow)
            let mutable n = 0
            let mutable sum = 0.
            let mutable sumMz = 0.
            let mutable sumIm = 0.
            let mutable sumImWeight = 0.
            while i < mz.Length && mz.[i] <= high do
                if not filterIonMobility || (im.[i] > imLow && im.[i] < imHigh) then
                    n <- n + 1
                    sum <- sum + intensity.[i]
                    sumMz <- sumMz + intensity.[i] * mz.[i]
                    if not (Double.IsNaN im.[i]) then
                        sumIm <- sumIm + intensity.[i] * im.[i]
                        sumImWeight <- sumImWeight + intensity.[i]
                record im.[i] intensity.[i]
                i <- i + 1
            if n > 0 && sum > 0. then
                sum, sumMz / sum, (if sumImWeight > 0. then sumIm / sumImWeight else nan)
            else
                sum, nan, nan

    /// One XIC point. Without ion mobility it is the peak closest to lockMz within +-halfWindow (Xic.closestPoint). With ion
    /// mobility it is the sum of all peaks within +-halfWindow of lockMz inside the ion mobility window (no window if imLow is
    /// nan), because the signal of an ion is spread over neighbouring m/z values and mobility scans.
    /// Returns intensity, intensity weighted m/z and intensity weighted ion mobility of these peaks; intensity 0 and nan for both
    /// if no peak is inside (nan without ion mobility).
    let extractPoint (spectrum: SortedPeaks.MzSortedPeaks) (lockMz: float) (halfWindow: float) (imLow: float) (imHigh: float) =
        extractPointWith (fun _ _ -> ()) spectrum lockMz halfWindow imLow imHigh

    /// Xic.extractPoint that also appends the point to the mobility profile, in the same pass over the peaks. The point sums the
    /// ion mobility window (imLow, imHigh) only; the profile records the pairs up to the half-width of the window beyond it on both
    /// sides (edges exclusive). The centre of an ion lies within that half-width of its expected ion mobility, so the window of that
    /// half-width around the centre lies inside the recorded range. imLow and imHigh must be finite.
    let extractPointWithProfile (profile: MobilityProfile) (spectrum: SortedPeaks.MzSortedPeaks) (lockMz: float) (halfWindow: float) (imLow: float) (imHigh: float) =
        let halfWidth = (imHigh - imLow) / 2.
        let recordLow = imLow - halfWidth
        let recordHigh = imHigh + halfWidth
        profile.StartPoint ()
        extractPointWith
            (fun im intensity -> if im > recordLow && im < recordHigh then profile.Add (im, intensity))
            spectrum lockMz halfWindow imLow imHigh

    /// The summed trace of a profile: for every point the sum of its pairs with centre - halfWidth < ion mobility < centre + halfWidth
    /// and their intensity weighted mean ion mobility (nan where the sum is 0).
    let summedTrace (profile: MobilityProfile) (centre: float) (halfWidth: float) =
        let low = centre - halfWidth
        let high = centre + halfWidth
        let intensity = Array.zeroCreate profile.Points
        let ionMobility = Array.create profile.Points nan
        for point in 0 .. profile.Points - 1 do
            let first, last = profile.Range point
            let mutable sum = 0.
            let mutable sumIm = 0.
            for k in first .. last - 1 do
                let im = profile.Mobility k
                if im > low && im < high then
                    let y = profile.Intensity k
                    sum <- sum + y
                    sumIm <- sumIm + y * im
            intensity.[point] <- sum
            if sum > 0. then ionMobility.[point] <- sumIm / sum
        intensity, ionMobility

    /// The points of a raw XIC trace whose ion mobility pairs give the mobility centre of an ion. From the point closest to the scan
    /// time the climb goes to the neighbour with the higher intensity while it is higher than the current point, which ends at the
    /// apex; the points are the contiguous run around the apex with at least half its intensity. No points for an empty trace or an
    /// apex that is not positive.
    let centrePoints (rt: float []) (intensity: float []) (scanTime: float) =
        if rt.Length = 0 then
            [||]
        else
            let at i = if i < 0 || i >= intensity.Length then -infinity else intensity.[i]
            let mutable apex = 0
            for i in 1 .. rt.Length - 1 do
                if abs (rt.[i] - scanTime) < abs (rt.[apex] - scanTime) then apex <- i
            let mutable climbing = true
            while climbing do
                let left = at (apex - 1)
                let right = at (apex + 1)
                if max left right > intensity.[apex] then
                    apex <- if left >= right then apex - 1 else apex + 1
                else
                    climbing <- false
            if not (intensity.[apex] > 0.) then
                [||]
            else
                let half = 0.5 * intensity.[apex]
                let mutable first = apex
                while at (first - 1) >= half do first <- first - 1
                let mutable last = apex
                while at (last + 1) >= half do last <- last + 1
                [| first .. last |]

    /// The intensity weighted mean ion mobility of the pairs of the given points of the profile that lie within halfWidth of the
    /// expected ion mobility (|ion mobility - expected| < halfWidth); nan if there is none.
    let mobilityCentre (profile: MobilityProfile) (points: int []) (expected: float) (halfWidth: float) =
        let mutable sum = 0.
        let mutable sumIm = 0.
        for point in points do
            let first, last = profile.Range point
            for k in first .. last - 1 do
                let im = profile.Mobility k
                if abs (im - expected) < halfWidth then
                    let intensity = profile.Intensity k
                    sum <- sum + intensity
                    sumIm <- sumIm + intensity * im
        if sum > 0. then sumIm / sum else nan

    /// The ion mobility centre of an ion from the signal around its expected ion mobility: Xic.mobilityCentre over the Xic.centrePoints
    /// of the raw XIC trace (scan times and intensities of the points of the profile) at the scan time; None if there is no pair.
    let tryMobilityCentre (profile: MobilityProfile) (rt: float []) (intensity: float []) (scanTime: float) (expected: float) (halfWidth: float) =
        let centre = mobilityCentre profile (centrePoints rt intensity scanTime) expected halfWidth
        if Numerics.isFinite centre then Some centre else None

    /// Counts the ions whose XIC is summed around the centre from the signal, with their centre minus the expected ion mobility,
    /// and the ions that keep the expected ion mobility because they have no centre.
    type MobilityCentreTally() =
        let shifts = ResizeArray<float>()
        let mutable kept = 0
        /// The number of ions counted.
        member _.Count = shifts.Count + kept
        /// Counts an ion with its expected ion mobility and its centre from the signal, if it has one.
        member _.Add(expected: float, centre: float option) =
            match centre with
            | Some centre -> shifts.Add (centre - expected)
            | None -> kept <- kept + 1
        /// The counts for a summed half-width, and the median, 5th and 95th percentile of centre minus expected ion mobility.
        member _.Summary(halfWidth: float) =
            let sorted = Array.sort (shifts.ToArray())
            sprintf "%i ions summed over +-%.4f around the centre of their signal, %i keep the expected ion mobility; centre minus expected ion mobility: median %.5f, 5th percentile %.5f, 95th percentile %.5f"
                sorted.Length halfWidth kept (Numerics.median sorted) (Numerics.quantileSorted sorted 0.05) (Numerics.quantileSorted sorted 0.95)

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
    /// +- mzHalfWindow). With an ion mobility window (imCenter +- imHalfWindow) only the peaks inside it count; imHalfWindow
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

    /// Xic.rtProfile with an ion mobility window (imCenter +- imHalfWindow) that also fills the mobility profile of its points
    /// (emptied first, Xic.extractPointWithProfile). imCenter must be finite.
    let rtProfileWithMobilityProfile (readSortedPeaks: string -> SortedPeaks.MzSortedPeaks) (rtIndex: IMzIOArray<RtIndexEntry>) (rtRange: RangeQuery) (lockMz: float) (mzHalfWindow: float) (imCenter: float) (imHalfWindow: float) (profile: MobilityProfile) =
        profile.Clear()
        RtIndexEntry.Search(rtIndex, rtRange)
        |> Seq.map (fun entry ->
            let intensity, _, _ = extractPointWithProfile profile (readSortedPeaks entry.SpectrumID) lockMz mzHalfWindow (imCenter - imHalfWindow) (imCenter + imHalfWindow)
            entry.Rt, intensity
            )
        |> Array.ofSeq
        |> Array.unzip
