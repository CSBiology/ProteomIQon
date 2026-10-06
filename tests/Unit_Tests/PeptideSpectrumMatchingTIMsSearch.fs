/// Tests of the search: the candidate order and the cut, the slice merge, the ceiling of the
/// index, the bin admission shared by the slice boundaries and the index build, the histogram
/// form of the expectation model, the streaming header parser, and the two library functions
/// BioFSharp.Mz 0.2.2 rewrote against copies of their originals.
module ProteomIQon.UnitTestsPeptideSpectrumMatchingTIMsSearch

open System
open System.Collections.Generic
open Expecto
open BioFSharp
open BioFSharp.Mz
open ProteomIQon
open ProteomIQon.PeptideIndex
open ProteomIQon.SpectrumSearch
open ProteomIQon.PeptideSpectrumMatchingTIMs

let private rng = Random 7

let private candidate rank decoy score exact =
    { Rank = rank; IsDecoy = decoy; MatchedB = 0; MatchedY = 0; Hyperscore = score; Exact = exact; ExactMatched = 0; Total = 0 }

/// Candidates with many ties: scores from a handful of values, exact scores from another handful.
let private tiedCandidates (count: int) =
    Array.init count (fun i ->
        candidate i (rng.Next 2 = 0) (float (rng.Next 5)) (float (rng.Next 3)))

let private shuffled (items: 'a[]) =
    let copy = Array.copy items
    for i = copy.Length - 1 downto 1 do
        let j = rng.Next (i + 1)
        let t = copy.[i]
        copy.[i] <- copy.[j]
        copy.[j] <- t
    copy

let private sortedBy (order: IComparer<Candidate>) (items: Candidate[]) =
    let copy = Array.copy items
    Array.Sort(copy, order)
    copy

let orderTests =
    testList "candidate order" [
        testCase "sorting gives the same array from any input order" <| fun _ ->
            let items = tiedCandidates 500
            let reference = sortedBy candidateOrder items
            for _ in 1 .. 20 do
                Expect.equal (sortedBy candidateOrder (shuffled items)) reference "the order is total"
        testCase "index pass score first, then exact score, then rank, then target before decoy" <| fun _ ->
            let a = candidate 5 false 3. 1.
            let b = candidate 2 false 2. 9.
            let c = candidate 9 false 3. 2.
            let d = candidate 5 true 3. 1.
            Expect.equal (sortedBy candidateOrder [| d; b; c; a |]) [| c; a; d; b |] "c by exact, a before d by label, b last by score"
        testCase "index order ignores the exact score" <| fun _ ->
            let a = candidate 5 false 3. 1.
            let c = candidate 9 false 3. 2.
            Expect.equal (sortedBy indexOrder [| c; a |]) [| a; c |] "rank decides"
    ]

let cutTests =
    testList "cutWithTies" [
        testCase "no tie at the cut gives exactly the limit" <| fun _ ->
            // exact scores of 0 and not NaN, since NaN is unequal to itself and would fail the comparison
            let items = Array.init 10 (fun i -> candidate i false (10. - float i) 0.)
            Expect.equal (cutWithTies items items.Length 4) (Array.sub items 0 4) "four"
        testCase "ties at the cut are kept" <| fun _ ->
            let items = Array.init 10 (fun i -> candidate i false (if i < 2 then 9. else 5.) nan)
            Expect.equal (cutWithTies items items.Length 4).Length 10 "the entire cutoff tie"
        testCase "cutoff ties are not capped before exact scoring" <| fun _ ->
            let items = Array.init 30 (fun i -> candidate i false 5. nan)
            Expect.equal (cutWithTies items items.Length 4).Length 30 "all ties must be considered"
        testCase "fewer than the limit come back as they are" <| fun _ ->
            let items = Array.init 3 (fun i -> candidate i false 5. 0.)
            Expect.equal (cutWithTies items items.Length 4) items "all three"
        testCase "a limit of zero keeps nothing" <| fun _ ->
            let items = Array.init 3 (fun i -> candidate i false 5. nan)
            Expect.equal (cutWithTies items items.Length 0) [||] "nothing"
        testCase "only the first count entries of the buffer are read" <| fun _ ->
            let items = Array.init 10 (fun i -> candidate i false 5. nan)
            Expect.equal (cutWithTies items 3 4).Length 3 "count bounds the buffer"
    ]

let mergeTests =
    testList "mergeTop" [
        testCase "merging slices in any partition gives the top of the whole" <| fun _ ->
            for _ in 1 .. 30 do
                let items = tiedCandidates (50 + rng.Next 200)
                let limit = 1 + rng.Next 12
                let reference = Array.truncate limit (sortedBy candidateOrder items)
                let slices = 1 + rng.Next 6
                let parts = shuffled items |> Array.mapi (fun i c -> i % slices, c) |> Array.groupBy fst |> Array.map (fun (_, xs) -> xs |> Array.map snd)
                let merged =
                    parts
                    |> Array.fold (fun kept part ->
                        // a slice delivers its own top, sorted and cut like the index pass does
                        let top = Array.truncate limit (sortedBy candidateOrder part)
                        mergeTop kept top limit) [||]
                Expect.equal merged reference "same kept set and order"
        testCase "the merge is bounded by the limit" <| fun _ ->
            let items = tiedCandidates 100
            Expect.isLessThanOrEqual (mergeTop (Array.sub items 0 50) (Array.sub items 50 50) 7).Length 7 "seven"
    ]

let private settings =
    {
        PrecursorTolerancePpm = 20.
        FragmentTolerancePpm = 20.
        IsotopeErrors = [| -1; 0; 1; 2 |]
        MaxFragmentCharge = 2
        TopNPeaks = 150
        MinimumRatio = 0.01
        RemovePrecursorRange = 1.5
        Deisotope = true
        MinimumPeaks = 5
        MinMatchedFragments = 4
        MinFragmentsModelling = 1
        ReportedHitsPerLabel = 10
    }

/// A small peptide table from residue strings over three residues, in mass order.
let private smallTable (sequences: string[]) =
    let codes =
        [|
            { Token = "G"; Masses = [| 57.02146; 57.02146 |] }
            { Token = "A"; Masses = [| 71.03711; 71.03711 |] }
            { Token = "K"; Masses = [| 128.09496; 128.09496 + 4. |] }
        |]
    let code (c: char) = codes |> Array.findIndex (fun r -> r.Token = string c) |> byte
    let residues = sequences |> Array.collect (fun s -> s.ToCharArray() |> Array.map code)
    let starts = Array.scan (fun start (s: string) -> start + s.Length) 0 sequences
    let massOf (s: string) = (s.ToCharArray() |> Array.sumBy (fun c -> codes.[int (code c)].Masses.[0])) + waterMass
    let order = Array.init sequences.Length id |> Array.sortBy (fun row -> massOf sequences.[row])
    {
        Count = sequences.Length
        Mass = order |> Array.map (fun row -> massOf sequences.[row])
        Order = order
        ModSequenceIdRow = Array.init sequences.Length id
        PepSequenceIdRow = Array.init sequences.Length id
        GlobalModRow = Array.zeroCreate sequences.Length
        SeqStartRow = starts
        Residues = residues
        Codes = codes
        MaxLength = sequences |> Array.map String.length |> Array.max
        RawSequenceRow = dict []
    }

let private table =
    smallTable [| "GAK"; "AAK"; "GGAK"; "AGAK"; "KKAG"; "GGGGK"; "AAAAK"; "KAKAK"; "GAGAGAK"; "AKAKAKAK"; "GGGGGGGGGGK"; "KKKKKKKK" |]

let ceilingTests =
    testList "fragment ceiling" [
        testCase "the reach of the windows is the widest isotope window" <| fun _ ->
            let precursor = 1500.
            let tol = precursor * 20e-6
            Expect.floatClose Accuracy.veryHigh (windowCeiling settings precursor) (precursor + isotopeSpacing + tol) "the -1 isotope error reaches highest"
        testCase "no isotope errors, no reach" <| fun _ ->
            Expect.equal (windowCeiling { settings with IsotopeErrors = [||] } 1500.) 0. "zero"
        testCase "every peptide of a precursor window is lighter than the reach" <| fun _ ->
            for _ in 1 .. 200 do
                let precursor = 200. + rng.NextDouble() * 1400.
                let reach = windowCeiling settings precursor
                for window in precursorWindows settings table precursor do
                    for rank in window.Lo .. window.Hi - 1 do
                        Expect.isLessThan table.Mass.[rank] reach "in the window means below the reach"
    ]

/// Entries the index build admits for a rank range, counted with the shared rule.
let private admitted (binWidth: float) (ceiling: float) (rankLo: int) (rankHi: int) =
    let binCount = binCountOf binWidth ceiling
    let buffer = Array.zeroCreate (2 * table.MaxLength)
    let mutable n = 0
    for rank = rankLo to rankHi - 1 do
        let count = fragmentMasses table rank buffer
        for i = 0 to count - 1 do
            if binOf binWidth binCount buffer.[i] >= 0 then n <- n + 1
    n

let sliceTests =
    testList "slices and the index" [
        testCase "the admission rule is the bin range" <| fun _ ->
            let binCount = binCountOf 0.02 100.
            Expect.equal binCount 5002 "ceiling over width plus two"
            Expect.equal (binOf 0.02 binCount 100.) 5000 "the ceiling itself is admitted"
            Expect.equal (binOf 0.02 binCount 100.03) 5001 "one bin of rounding room"
            Expect.equal (binOf 0.02 binCount 100.04) -1 "beyond the room is out"
            Expect.equal (binOf 0.02 binCount -0.5) -1 "negative masses are out"
        testCase "slice boundaries start at zero, end at the count and do not decrease" <| fun _ ->
            for slices in 1 .. 6 do
                let boundaries = sliceBoundaries table slices 0.01 5000.
                Expect.equal boundaries.[0] 0 "zero"
                Expect.equal boundaries.[slices] table.Count "count"
                Expect.isTrue (boundaries |> Array.pairwise |> Array.forall (fun (a, b) -> a <= b)) "monotone"
        testCase "slices hold about the same number of admitted entries under a ceiling" <| fun _ ->
            let ceiling = 500.
            let boundaries = sliceBoundaries table 3 0.01 ceiling
            let counts = Array.init 3 (fun s -> admitted 0.01 ceiling boundaries.[s] boundaries.[s + 1])
            let total = admitted 0.01 ceiling 0 table.Count
            for count in counts do
                Expect.isLessThanOrEqual (abs (count - total / 3)) (2 * table.MaxLength) "within one peptide of the share"
        testCase "the index of a slice holds exactly the admitted entries" <| fun _ ->
            for ceiling in [ 150.; 300.; 500.; 5000. ] do
                let boundaries = sliceBoundaries table 3 0.01 ceiling
                for s in 0 .. 2 do
                    let index = buildFragmentIndex table 0.01 ceiling 2 boundaries.[s] boundaries.[s + 1] ignore
                    Expect.equal index.EntryCount (int64 (admitted 0.01 ceiling boundaries.[s] boundaries.[s + 1])) "same rule"
        testCase "one index and the union of the slices hold the same entries" <| fun _ ->
            let entriesOf (index: FragmentIndex) =
                [ for bin in 0 .. index.BinCount - 1 do
                    for i in 0 .. index.BinLength.[bin] - 1 do
                        yield bin, index.Chunks.[index.BinChunk.[bin]].[index.BinLocal.[bin] + i] ]
            let whole = buildFragmentIndex table 0.01 400. 1 0 table.Count ignore |> entriesOf
            let boundaries = sliceBoundaries table 4 0.01 400.
            let parts = [ for s in 0 .. 3 do yield! buildFragmentIndex table 0.01 400. 3 boundaries.[s] boundaries.[s + 1] ignore |> entriesOf ]
            Expect.equal (List.sort parts) (List.sort whole) "same entries"
    ]

let histogramTests =
    testList "expectation model from a histogram" [
        testCase "the histogram of two slices gives the one pass model" <| fun _ ->
            for _ in 1 .. 20 do
                let scores = Array.init (50 + rng.Next 2000) (fun _ -> rng.NextDouble() * 30.)
                let onePass = expectationModel scores
                let half = scores.Length / 2
                let histogram = Array.zeroCreate 64 |> ref
                let add (part: float[]) =
                    for score in part do
                        let bin = histogramBin score
                        if bin + 2 > histogram.Value.Length then
                            let grown = Array.zeroCreate (max (bin + 2) (histogram.Value.Length * 2))
                            Array.blit histogram.Value 0 grown 0 histogram.Value.Length
                            histogram.Value <- grown
                        histogram.Value.[bin] <- histogram.Value.[bin] + 1
                add (Array.sub scores 0 half)
                add (Array.sub scores half (scores.Length - half))
                let last = histogram.Value |> Array.findIndexBack (fun c -> c > 0)
                let sliced = expectationModelOfHistogram (Array.sub histogram.Value 0 (last + 2))
                for score in [ 0.; 5.; 12.5; 20.; 35. ] do
                    Expect.equal (sliced score) (onePass score) "same expectation"
    ]

let headerTests =
    let sample = IO.File.ReadAllText(IO.Path.Combine(AppContext.BaseDirectory, "sample_description.json"))
    testList "streaming header parser" [
        testCase "a real description" <| fun _ ->
            let h = TimSpectrumHeader.tryParse sample
            Expect.isNull h.Error "no error"
            Expect.equal h.ID "merged=253418 frame=35072 scanStart=763 scanEnd=787" "id"
            Expect.equal h.MsLevel (Some 2) "ms level"
            Expect.equal h.ChargeState (Some 3) "charge"
            Expect.equal h.PrecursorMz (Some 419.215935113248) "selected ion m/z, MS:1000744"
            Expect.equal h.IonMobility (Some 0.737820952224) "mobility"
            Expect.equal h.ScanTime (Some 62.17991573) "scan time"
        testCase "the isolation target MS:1002234 wins over the selected ion m/z" <| fun _ ->
            let marker = "\"MS:1000041\":{\"$id\":\"1\",\"CvAccession\":\"MS:1000041\""
            Expect.equal (sample.Split(marker).Length) 2 "the sample carries the charge once"
            let json = sample.Replace(marker, "\"MS:1002234\":{\"$id\":\"1\",\"CvAccession\":\"MS:1002234\",\"Type\":\"WithCvUnitAccession\",\"Values\":[\"420.5\",\"MS:1000040\"]}," + marker)
            let h = TimSpectrumHeader.tryParse json
            Expect.isNull h.Error "no error"
            Expect.equal h.PrecursorMz (Some 420.5) "target"
            Expect.equal h.ChargeState (Some 3) "charge still read"
        testCase "missing scan, precursor and charge give None and no error" <| fun _ ->
            let h = TimSpectrumHeader.tryParse """{"$id":"1","properties":{"$id":"2","MS:1000511":{"CvAccession":"MS:1000511","Type":"CvValue","Values":["1"]}},"ID":"x"}"""
            Expect.isNull h.Error "no error"
            Expect.equal h.MsLevel (Some 1) "ms level"
            Expect.equal h.ScanTime None "scan time"
            Expect.equal h.PrecursorMz None "precursor"
            Expect.equal h.ChargeState None "charge"
            Expect.equal h.IonMobility None "mobility"
        testCase "a malformed description names the failure and carries nothing else" <| fun _ ->
            let h = TimSpectrumHeader.tryParse (sample.Substring(0, sample.Length / 2))
            Expect.isNotNull h.Error "error"
            Expect.isNull h.ID "no id"
            Expect.equal h.MsLevel None "no level"
    ]

/// The autocorrelation of the released library, copied from BioFSharp.Mz 0.2.1.
let private originalAutoCorrelation (plusMinusMaxDelay: int) (vector: FSharp.Stats.Vector<float>) =
    let shifted (vector: FSharp.Stats.Vector<float>) (tau: int) =
        vector |> FSharp.Stats.Vector.mapi (fun i _ ->
            let index = i - tau
            if index < 0 || index > vector.Length - 1 then 0. else vector.[index])
    let rec accumVector accum (state: int) (max: int) =
        if state = max then accum else accumVector (accum + shifted vector state) (state - 1) max
    let empty = FSharp.Stats.Vector.zero vector.Length
    let plus = accumVector empty plusMinusMaxDelay 1
    let minus = accumVector empty -1 (-plusMinusMaxDelay)
    (plus + minus) |> FSharp.Stats.Vector.map (fun x -> x / (float plusMinusMaxDelay * 2.))

/// The theoretical spectrum builder of the released library, copied from BioFSharp.Mz 0.2.1.
let private originalPredictOf (lowerScanLimit, upperScanLimit) chargeState (fragments: PeakFamily<TaggedMass.TaggedMass> list) =
    let predictPeak charge (taggedMass: TaggedMass.TaggedMass) =
        TaggedPeak.TaggedPeak(taggedMass.Iontype, Mass.toMZ taggedMass.Mass charge, nan)
    let computePeakFamily charge fragments =
        let mainPeak = predictPeak charge fragments.MainPeak
        let dependentPeaks =
            fragments.DependentPeaks
            |> List.fold (fun acc dependent -> if charge <= 1. then predictPeak charge dependent :: acc else acc) []
        Peaks.createPeakFamily mainPeak dependentPeaks
    let rec recloop ions charge (fragments: PeakFamily<TaggedMass.TaggedMass> list) =
        match fragments with
        | fragments :: rest ->
            if charge <= 1. then recloop (computePeakFamily 1. fragments :: ions) charge rest
            else
                let tempIons = [ for z = 1 to 2 do yield computePeakFamily (float z) fragments ]
                recloop (tempIons @ ions) charge rest
        | [] -> ions
    recloop [] chargeState fragments |> List.toArray

let libraryTests =
    testList "rewritten library functions" [
        testCase "autocorrelation equals the recursion of the released library" <| fun _ ->
            for _ in 1 .. 50 do
                let n = 5 + rng.Next 300
                let vector = FSharp.Stats.Vector.init n (fun _ -> if rng.Next 3 = 0 then 0. else rng.NextDouble() * 100.)
                let delay = 1 + rng.Next 80
                let expected = originalAutoCorrelation delay vector |> Seq.toArray
                let actual = SequestLike.autoCorrelation delay vector |> Seq.toArray
                Expect.equal actual expected "bitwise the same vector"
        testCase "predictOf gives the families of the released library in the same order" <| fun _ ->
            for peptide in [ "PEPTIDEK"; "GAGAGAK"; "MKWVTFISLLFLFSSAYSR"; "AK" ] do
                let aal = BioList.ofAminoAcidString peptide
                let series = Fragmentation.Series.bOfBioList BioItem.monoisoMass aal @ Fragmentation.Series.yOfBioList BioItem.monoisoMass aal
                for charge in [ 1.; 2.; 3. ] do
                    let expected = originalPredictOf (0., 2000.) charge series
                    let actual = XScoring.predictOf (0., 2000.) charge series
                    Expect.equal actual.Length expected.Length "count"
                    for i in 0 .. expected.Length - 1 do
                        let a : PeakFamily<TaggedPeak.TaggedPeak> = actual.[i]
                        let e : PeakFamily<TaggedPeak.TaggedPeak> = expected.[i]
                        Expect.equal a.MainPeak.Mz e.MainPeak.Mz "main peak m/z"
                        Expect.equal a.MainPeak.Iontype e.MainPeak.Iontype "ion type"
                        Expect.equal (a.DependentPeaks |> List.map (fun (p: TaggedPeak.TaggedPeak) -> p.Mz)) (e.DependentPeaks |> List.map (fun (p: TaggedPeak.TaggedPeak) -> p.Mz)) "dependent peaks"
    ]

[<Tests>]
let all =
    testList "PeptideSpectrumMatchingTIMs search" [ orderTests; cutTests; mergeTests; ceilingTests; sliceTests; histogramTests; headerTests; libraryTests ]
