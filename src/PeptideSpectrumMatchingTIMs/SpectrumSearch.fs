namespace ProteomIQon

open System
open BioFSharp
open BioFSharp.Mz
open FSharp.Stats
open PeptideIndex

/// Spectrum preprocessing and the fragment index search for one MS2 spectrum, with the
/// hyperscore and the expectation value model built from all candidates of that spectrum.
module SpectrumSearch =

    let protonMass = Mass.Table.PMassInU

    /// Mass difference between consecutive isotope peaks.
    let isotopeSpacing = Isotopes.Table.C13.Mass - Isotopes.Table.C12.Mass

    type SearchSettings =
        {
            /// Precursor tolerance in ppm of the precursor neutral mass.
            PrecursorTolerancePpm : float
            /// Fragment tolerance in ppm of the fragment m/z.
            FragmentTolerancePpm  : float
            /// Isotope errors of the precursor to consider, for example [0; 1; 2].
            IsotopeErrors         : int[]
            /// Highest fragment charge state to match, at most the precursor charge minus one and at least 1.
            MaxFragmentCharge     : int
            /// Only the most intense peaks are used for matching.
            TopNPeaks             : int
            /// Peaks below this fraction of the base peak are discarded.
            MinimumRatio          : float
            /// Peaks within +/- this m/z of the precursor m/z are discarded.
            RemovePrecursorRange  : float
            /// Removes peaks that sit one or more isotope spacings above a kept peak.
            Deisotope             : bool
            /// Spectra with fewer processed peaks are skipped.
            MinimumPeaks          : int
            /// Candidates need at least this many matched fragments to be reported. The index
            /// pass applies it to the binned matches, the exact rescoring applies it again.
            MinMatchedFragments   : int
            /// Candidates need at least this many matched fragments to enter the expectation model.
            MinFragmentsModelling : int
            /// Number of best targets and best decoys kept per spectrum.
            ReportedHitsPerLabel  : int
        }

    /// A processed peak list sorted by m/z with intensities scaled to a base peak of 100.
    type ProcessedSpectrum =
        {
            Mz        : float[]
            Intensity : float[]
        }

    /// One candidate of the index search.
    [<Struct>]
    type Candidate =
        {
            Rank       : int
            IsDecoy    : bool
            MatchedB   : int
            MatchedY   : int
            Hyperscore : float
        }

    /// A precursor mass window as a rank range [Lo, Hi) together with its offset in the
    /// concatenated accumulator of all windows of a spectrum.
    [<Struct>]
    type Window =
        {
            Lo     : int
            Hi     : int
            Offset : int
        }

    /// X!Tandem hyperscore as BioFSharp.Mz computes it: ln(sum of matched intensities) + ln(Nb!)
    /// + ln(Ny!), intensities relative to a base peak of 100.
    let hyperscore (matchedB: int) (matchedY: int) (sumB: float) (sumY: float) =
        XScoring.createCountedMatches 0 0 0 0 matchedB matchedY 0 (sumB + sumY)
        |> XScoring.calcHyperScore

    /// Index of the first peak at or above the given m/z, or the array length.
    let private firstAtOrAbove (mz: float[]) (value: float) =
        let rec loop l h =
            if l >= h then l
            else
                let m = l + ((h - l) >>> 1)
                if mz.[m] < value then loop (m + 1) h else loop l m
        loop 0 mz.Length

    /// Marks the peaks that sit one or more isotope spacings above a kept peak, walking each
    /// envelope upwards from its lowest member so that the monoisotopic peak survives. Peaks
    /// must be sorted by m/z. Runs for every spectrum, so it works on a flag array.
    let private isotopePartners (settings: SearchSettings) (precursorCharge: int) (mz: float[]) =
        let removed = Array.zeroCreate<bool> mz.Length
        let partnerOf (position: int) (charge: int) (step: int) =
            let expected = mz.[position] + float step * isotopeSpacing / float charge
            let tol = expected * settings.FragmentTolerancePpm * 1e-6
            let i = firstAtOrAbove mz (expected - tol)
            if i > position && i < mz.Length && mz.[i] <= expected + tol then i else -1
        let rec envelope (position: int) (charge: int) (step: int) =
            let partner = partnerOf position charge step
            if partner >= 0 then
                removed.[partner] <- true
                envelope position charge (step + 1)
        for i = 0 to mz.Length - 1 do
            if not removed.[i] then
                for charge = 1 to max 1 precursorCharge do
                    envelope i charge 1
        removed

    /// Removes precursor peaks, weak peaks and isotope partners, keeps the top N peaks and
    /// returns them sorted by m/z with intensities scaled to a base peak of 100.
    let preprocess (settings: SearchSettings) (precursorMz: float) (precursorCharge: int) (mz: float[]) (intensity: float[]) =
        let peaks =
            Array.zip mz intensity
            |> Array.filter (fun (m, i) -> abs (m - precursorMz) > settings.RemovePrecursorRange && i > 0. && not (Double.IsNaN i))
            |> Array.sortBy fst
        let basePeak = if peaks.Length = 0 then 0. else peaks |> Array.maxBy snd |> snd
        let strong = peaks |> Array.filter (fun (_, i) -> i >= basePeak * settings.MinimumRatio)
        let deisotoped =
            if settings.Deisotope then
                let partners = isotopePartners settings precursorCharge (Array.map fst strong)
                strong |> Array.indexed |> Array.filter (fun (i, _) -> not partners.[i]) |> Array.map snd
            else strong
        let selected =
            deisotoped
            |> Array.sortByDescending snd
            |> Array.truncate settings.TopNPeaks
            |> Array.sortBy fst
        let scale = if basePeak > 0. then 100. / basePeak else 1.
        {
            Mz = selected |> Array.map fst
            Intensity = selected |> Array.map (fun (_, i) -> i * scale)
        }

    /// Accumulators of matched fragments per candidate, reused across spectra so that the
    /// scatter loop allocates nothing. Slot index: (window offset + rank - window start) * 2 + decoy.
    type Scratch() =
        let mutable countB : int[] = Array.zeroCreate 1024
        let mutable countY : int[] = Array.zeroCreate 1024
        let mutable sumB : float[] = Array.zeroCreate 1024
        let mutable sumY : float[] = Array.zeroCreate 1024
        let mutable fragments : float[] = Array.zeroCreate 512
        member _.CountB = countB
        member _.CountY = countY
        member _.SumB = sumB
        member _.SumY = sumY
        member _.Fragments = fragments
        /// Clears the slots of all windows and grows the arrays when needed.
        member _.Reset (slots: int) (maxFragments: int) =
            if slots > countB.Length then
                let capacity = max slots (countB.Length * 2)
                countB <- Array.zeroCreate capacity
                countY <- Array.zeroCreate capacity
                sumB <- Array.zeroCreate capacity
                sumY <- Array.zeroCreate capacity
            else
                Array.Clear(countB, 0, slots)
                Array.Clear(countY, 0, slots)
                Array.Clear(sumB, 0, slots)
                Array.Clear(sumY, 0, slots)
            if maxFragments > fragments.Length then fragments <- Array.zeroCreate maxFragments

    /// Precursor rank windows for every isotope error, merged where they overlap and laid out
    /// one after the other in the accumulator.
    let precursorWindows (settings: SearchSettings) (table: PeptideTable) (precursorMass: float) =
        let tol = precursorMass * settings.PrecursorTolerancePpm * 1e-6
        let merged =
            settings.IsotopeErrors
            |> Array.map (fun k ->
                let center = precursorMass - float k * isotopeSpacing
                rankLowerBound table (center - tol), rankLowerBound table (center + tol))
            |> Array.filter (fun (lo, hi) -> hi > lo)
            |> Array.sortBy fst
            |> Array.fold (fun acc (lo, hi) ->
                match acc with
                | (plo, phi) :: rest when lo <= phi -> (plo, max phi hi) :: rest
                | _ -> (lo, hi) :: acc) []
            |> List.rev
        let offsets = merged |> List.scan (fun offset (lo, hi) -> offset + (hi - lo)) 0
        List.zip merged (List.take merged.Length offsets)
        |> List.map (fun ((lo, hi), offset) -> { Lo = lo; Hi = hi; Offset = offset })
        |> Array.ofList

    /// Number of accumulator slots the windows need.
    let slotCount (windows: Window[]) =
        if windows.Length = 0 then 0
        else
            let last = windows.[windows.Length - 1]
            2 * (last.Offset + last.Hi - last.Lo)

    /// Neutral mass intervals a spectrum matches: every peak at every fragment charge gives
    /// the interval of its tolerance, and intervals closer than one index bin are merged with
    /// the larger intensity. The merged intervals cover disjoint bins even after a shift by
    /// the water mass, so every index entry is counted at most once per pass, the way the
    /// exact rescoring counts every theoretical fragment once with its closest peak.
    let massIntervals (settings: SearchSettings) (binWidth: float) (maxFragmentCharge: int) (spectrum: ProcessedSpectrum) =
        Array.init (spectrum.Mz.Length * maxFragmentCharge) (fun i ->
            let mz = spectrum.Mz.[i / maxFragmentCharge]
            let z = float (i % maxFragmentCharge + 1)
            let mass = (mz - protonMass) * z
            let tol = mz * settings.FragmentTolerancePpm * 1e-6 * z
            struct (mass - tol, mass + tol, spectrum.Intensity.[i / maxFragmentCharge]))
        |> Array.sortBy (fun struct (lo, _, _) -> lo)
        |> Array.fold (fun acc (struct (lo, hi, intensity)) ->
            match acc with
            | struct (plo, phi, pIntensity) :: rest when lo < phi + binWidth ->
                struct (plo, max phi hi, max pIntensity intensity) :: rest
            | _ -> struct (lo, hi, intensity) :: acc) []
        |> List.rev
        |> Array.ofList

    /// Scatters the mass intervals of a spectrum over the fragment index and accumulates the
    /// matched b and y ions per candidate. Decoy candidates are served from the target index:
    /// a decoy b ion of mass m is a target y ion of mass m + water, a decoy y ion of mass m is
    /// a target b ion of mass m - water. This loop runs for every interval of every spectrum,
    /// so it uses index arithmetic and mutable cursors.
    let scatter (settings: SearchSettings) (index: FragmentIndex) (scratch: Scratch) (windows: Window[]) (precursorCharge: int) (spectrum: ProcessedSpectrum) =
        let maxFragmentCharge = max 1 (min settings.MaxFragmentCharge (precursorCharge - 1))
        let countB = scratch.CountB
        let countY = scratch.CountY
        let sumB = scratch.SumB
        let sumY = scratch.SumY
        // Adds the index entries of the neutral mass interval to the accumulators of the slot
        // parity decoy. takeB and takeY select the entry types that count. For decoys the ion
        // type flips: a target y entry counts as a decoy b ion and a target b entry as a decoy y ion.
        let accumulate (lo: float) (hi: float) (intensity: float) (decoy: int) (takeB: bool) (takeY: bool) =
            let b0 = max 0 (int (lo / index.BinWidth))
            let b1 = min (index.BinCount - 1) (int (hi / index.BinWidth))
            for b = b0 to b1 do
                let len = index.BinLength.[b]
                if len > 0 then
                    let chunk = index.Chunks.[index.BinChunk.[b]]
                    let s = index.BinLocal.[b]
                    let e = s + len
                    for w = 0 to windows.Length - 1 do
                        let window = windows.[w]
                        let mutable i = lowerBound chunk s e window.Lo
                        while i < e && (chunk.[i] &&& RankMask) < window.Hi do
                            let entry = chunk.[i]
                            let isY = entry &&& YFlag <> 0
                            if (isY && takeY) || (not isY && takeB) then
                                let slot = ((window.Offset + (entry &&& RankMask) - window.Lo) <<< 1) ||| decoy
                                if isY = (decoy = 0) then
                                    countY.[slot] <- countY.[slot] + 1
                                    sumY.[slot] <- sumY.[slot] + intensity
                                else
                                    countB.[slot] <- countB.[slot] + 1
                                    sumB.[slot] <- sumB.[slot] + intensity
                            i <- i + 1
        massIntervals settings index.BinWidth maxFragmentCharge spectrum
        |> Array.iter (fun struct (lo, hi, intensity) ->
            accumulate lo hi intensity 0 true true
            accumulate (lo + waterMass) (hi + waterMass) intensity 1 false true
            accumulate (lo - waterMass) (hi - waterMass) intensity 1 true false)

    /// Collects the candidates of the searched windows. Returns the hyperscores of all
    /// candidates that enter the expectation model and the best targets and decoys. The two
    /// thresholds are independent, a candidate can be reported without entering the model
    /// and enter the model without being reported.
    let collect (settings: SearchSettings) (scratch: Scratch) (windows: Window[]) =
        let minMatched = min settings.MinFragmentsModelling settings.MinMatchedFragments
        let candidates =
            [|
                for window in windows do
                    for rank in window.Lo .. window.Hi - 1 do
                        for decoy in 0 .. 1 do
                            let slot = ((window.Offset + rank - window.Lo) <<< 1) ||| decoy
                            let nb = scratch.CountB.[slot]
                            let ny = scratch.CountY.[slot]
                            if nb + ny >= minMatched then
                                yield
                                    {
                                        Rank = rank
                                        IsDecoy = decoy = 1
                                        MatchedB = nb
                                        MatchedY = ny
                                        Hyperscore = hyperscore nb ny scratch.SumB.[slot] scratch.SumY.[slot]
                                    }
            |]
        let scores =
            candidates
            |> Array.filter (fun c -> c.MatchedB + c.MatchedY >= settings.MinFragmentsModelling)
            |> Array.map (fun c -> c.Hyperscore)
        let top isDecoy =
            candidates
            |> Array.filter (fun c -> c.IsDecoy = isDecoy && c.MatchedB + c.MatchedY >= settings.MinMatchedFragments)
            |> Array.sortByDescending (fun c -> c.Hyperscore)
            |> Array.truncate settings.ReportedHitsPerLabel
        scores, top false, top true

    /// Hyperscores go into the survival histogram on the X!Tandem scale (4 log10 of the
    /// product form, which is 4 / ln 10 of the natural log form used here).
    let private histogramScale = 4. / log 10.

    /// Expectation value model of X!Tandem (mhistogram::survival and mhistogram::model), which
    /// MSFragger follows. The candidate hyperscores form a histogram with unit bins on the
    /// X!Tandem scale, rounded to the nearest bin. Its survival function is cleaned of
    /// candidates that sit above a gap in the upper part of the distribution (the potentially
    /// valid matches), then a line is fitted to log10 of the survival function between half of
    /// the remaining count and ten survivors. With fewer than 200 remaining candidates
    /// X!Tandem's default line is used. Returns a function from a hyperscore to the expected
    /// number of random candidates scoring at least as high, floored at 1e-15. The scores put
    /// into the model and the score of a hit must come from the same matching, so callers
    /// evaluate it with the index pass hyperscore of the hit.
    let expectationModel (scores: float[]) =
        let expect (a0: float) (a1: float) (score: float) =
            max 1e-15 (Math.Pow(10., a0 + a1 * score * histogramScale))
        let bins = scores |> Array.map (fun s -> max 0 (int (s * histogramScale + 0.5)))
        let histogram = Array.zeroCreate<int> (if bins.Length = 0 then 1 else Array.max bins + 2)
        bins |> Array.iter (fun b -> histogram.[b] <- histogram.[b] + 1)
        let raw = Array.scanBack (+) histogram 0 |> Array.take histogram.Length
        let n = raw.[0]
        if n = 0 then expect 3.5 -0.18
        else
            // Walk down from the top. A plateau (equal neighbours) above the point where the
            // survival function drops below a fifth of the count marks candidates that are
            // separated from the random distribution by an empty bin. Their count is removed
            // from every lower value.
            let mid = raw |> Array.findIndex (fun v -> v <= n / 5)
            let top = raw |> Array.findIndexBack (fun v -> v > 0)
            let rec clean a (removed: int) (acc: int list) =
                if a < 0 then acc
                elif a > 0 && raw.[a] = raw.[a - 1] && raw.[a] <> raw.[0] && a > mid then
                    let plateau = raw.[a]
                    let rec skip b (acc: int list) =
                        if b >= 0 && raw.[b] = plateau then skip (b - 1) ((raw.[b] - plateau) :: acc)
                        else b, acc
                    let b, acc = skip a acc
                    clean b plateau acc
                else clean (a - 1) removed ((raw.[a] - removed) :: acc)
            let survival = clean top 0 [] |> Array.ofList
            let survival = Array.append survival (Array.zeroCreate (raw.Length - survival.Length))
            let maxLimit = int (0.5 + float survival.[0] / 2.)
            let points =
                [| 0 .. survival.Length - 2 |]
                |> Array.skipWhile (fun a -> survival.[a] > maxLimit)
                |> Array.takeWhile (fun a -> survival.[a] > 10)
            if survival.[0] < 200 || points.Length = 0 then expect 3.5 -0.18
            else
                // The survival function does not increase, so the fit starts at the first point.
                let xs = points |> Array.map float |> vector
                let ys = points |> Array.map (fun a -> log10 (float survival.[a])) |> vector
                if xs.Length < 2 then expect 3.5 -0.18
                else
                    let coefficients = Fitting.LinearRegression.OLS.Linear.Univariable.fit xs ys
                    if Double.IsNaN coefficients.Linear then expect 3.5 -0.18
                    else expect coefficients.Constant coefficients.Linear

    /// Exact rescoring of one candidate against the processed spectrum with the ppm tolerance,
    /// independent of the index bin width. Returns matched b, matched y, the matched intensity
    /// sums of both series and the number of theoretical fragments. Runs for every reported
    /// hit, so it scans the peak array with index arithmetic.
    let rescore (settings: SearchSettings) (table: PeptideTable) (scratch: Scratch) (precursorCharge: int)
                (spectrum: ProcessedSpectrum) (rank: int) (isDecoy: bool) =
        let maxFragmentCharge = max 1 (min settings.MaxFragmentCharge (precursorCharge - 1))
        let buffer = scratch.Fragments
        let count = fragmentMasses table rank buffer
        let half = count / 2
        let mz = spectrum.Mz
        // Intensity of the closest peak within tolerance over the fragment charge states, or -1.
        let matchedIntensity (mass: float) =
            let rec overCharges z (bestDelta: float) (best: float) =
                if z > maxFragmentCharge then best
                else
                    let target = mass / float z + protonMass
                    let tol = target * settings.FragmentTolerancePpm * 1e-6
                    let rec scan i (bestDelta: float) (best: float) =
                        if i < mz.Length && mz.[i] <= target + tol then
                            let delta = abs (mz.[i] - target)
                            if delta < bestDelta then scan (i + 1) delta spectrum.Intensity.[i]
                            else scan (i + 1) bestDelta best
                        else bestDelta, best
                    let bestDelta, best = scan (firstAtOrAbove mz (target - tol)) bestDelta best
                    overCharges (z + 1) bestDelta best
            overCharges 1 Double.MaxValue -1.
        let series (fragmentMass: int -> float) =
            Seq.init half fragmentMass
            |> Seq.map matchedIntensity
            |> Seq.filter (fun intensity -> intensity >= 0.)
            |> Seq.fold (fun (n, sum) intensity -> n + 1, sum + intensity) (0, 0.)
        let nb, sumB = series (fun i -> if isDecoy then buffer.[half + i] - waterMass else buffer.[i])
        let ny, sumY = series (fun i -> if isDecoy then buffer.[i] + waterMass else buffer.[half + i])
        nb, ny, sumB, sumY, count
