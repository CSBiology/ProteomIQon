/// Tests of the ion mobility profile of XIC points, the ion mobility centre of an ion and the XICs summed around it.
module ProteomIQon.UnitTestsIonMobilityCentre

open System
open Expecto
open ProteomIQon
open MzIO.Commons.Arrays
open MzIO.Processing.MzIOLinq
open FSharp.Stats

let private spectrum (rows: (float * float * float) list) =
    SortedPeaks.sortPeakArraysByMz
        (rows |> List.map (fun (m, _, _) -> m) |> Array.ofList)
        (rows |> List.map (fun (_, i, _) -> i) |> Array.ofList)
        (rows |> List.map (fun (_, _, im) -> im) |> Array.ofList)

let private bits (x: float) = BitConverter.DoubleToInt64Bits x

/// A profile with the given points; every point is a list of (ion mobility, intensity).
let private profileOf (points: (float * float) list list) =
    let profile = Xic.MobilityProfile()
    for point in points do
        profile.StartPoint ()
        for (im, intensity) in point do
            profile.Add (im, intensity)
    profile

/// The pairs of a point of a profile, sorted.
let private pairsOf (profile: Xic.MobilityProfile) (point: int) =
    let first, last = profile.Range point
    [ for k in first .. last - 1 -> profile.Mobility k, profile.Intensity k ] |> List.sort

/// The peaks (m/z, intensity, mobility) of the spectrum i of a species at the m/z `mz` that elutes as a Gaussian (centre and width in
/// points) with its intensity, split on the mobility scans of `scans` (mobility, share).
let private speciesAt (mz: float) (i: int) (intensity: float) (centre: float) (width: float) (scans: (float * float) list) =
    let level = intensity * exp (-((float i - centre) ** 2.) / (2. * width * width))
    scans |> List.map (fun (im, share) -> mz, Math.Round (level * share), im)

let private species = speciesAt 500.0

let private spacing = 0.04

let private ms1Of (spectra: SortedPeaks.MzSortedPeaks []) = spectra |> Array.mapi (fun i _ -> float i * spacing, string i)

/// Sweeps the spectra with one request and returns its raw XIC and the profiles of its traces.
let private sweepOne (spectra: SortedPeaks.MzSortedPeaks []) (request: XicSweep.XicRequest) =
    let result = ref None
    XicSweep.sweep (fun id -> spectra.[int id]) (ms1Of spectra) [| request |] (fun _ raw profiles -> result.Value <- Some (raw, profiles)) |> ignore
    result.Value.Value

/// The request of the traces at the given m/z with the ion mobility window expected +- halfWidth, recording the profiles.
let private requestFor (lockMz: float []) expected halfWidth : XicSweep.XicRequest =
    {
        RtLow = 0.
        RtHigh = 1000.
        LockMz = lockMz
        MzWindow = 0.02
        ImLow = expected - halfWidth
        ImHigh = expected + halfWidth
        Point = XicSweep.SumWithProfiles
    }

/// The values of the traces of an XIC as bits, so that nan equals nan.
let private bitsOf (traces: float [] []) = traces |> Array.map (Array.map bits)

let private tolerance = { absolute = 1e-6; relative = 0. }

let private tolerance4 = { absolute = 1e-4; relative = 0. }

[<Tests>]
let tests =
    testList "Ion mobility centre" [

        testCase "the point with a profile is the point of extractPoint, bit for bit, and the profile holds the pairs of the point up to the half-width of the window beyond it" <| fun _ ->
            let rnd = Random 11
            for _ in 1 .. 300 do
                let n = rnd.Next(0, 60)
                let mzs = Array.init n (fun _ -> 500. + float (rnd.Next(0, 400)) * 0.001 + 1e-7)
                let intensities = Array.init n (fun _ -> float (rnd.Next(1, 1000)))
                let ims = Array.init n (fun _ -> if rnd.Next 15 = 0 then nan else 0.8 + rnd.NextDouble() * 0.4)
                let peaks = SortedPeaks.sortPeakArraysByMz mzs intensities ims
                let lock = 500.2 + rnd.NextDouble() * 0.0003
                let halfWindow = 0.05 + rnd.NextDouble() * 0.05
                let halfWidth = [| 0.02; 0.05; 0.1 |].[rnd.Next 3]
                let imLow, imHigh = 1.0 - halfWidth, 1.0 + halfWidth
                let profile = Xic.MobilityProfile()
                let i1, m1, im1 = Xic.extractPointWithProfile profile peaks lock halfWindow imLow imHigh
                let i2, m2, im2 = Xic.extractPoint peaks lock halfWindow imLow imHigh
                Expect.equal (bits i1, bits m1, bits im1) (bits i2, bits m2, bits im2) "same point"
                let expected =
                    [ for k in 0 .. peaks.Mz.Length - 1 do
                        let inMz = peaks.Mz.[k] >= lock - halfWindow && peaks.Mz.[k] <= lock + halfWindow
                        let inRange = peaks.IonMobility.[k] > 1.0 - 2. * halfWidth && peaks.IonMobility.[k] < 1.0 + 2. * halfWidth
                        if inMz && inRange then yield float (float32 peaks.IonMobility.[k]), float (float32 peaks.Intensity.[k]) ]
                Expect.equal (pairsOf profile 0) (List.sort expected) "pairs of the point"

        testCase "the profile takes the m/z edges, leaves out the edges of the recorded range and the peaks without mobility, and has one start per point" <| fun _ ->
            // the window is 0.5 .. 1.5, so the recorded range is 0 .. 2
            let peaks =
                spectrum [ 500.25, 1., 1.0     // exactly at the lower m/z edge: inside
                           500.75, 2., 1.0     // exactly at the upper m/z edge: inside
                           500.76, 4., 1.0     // outside the m/z window
                           500.50, 8., 0.5     // at the edge of the window: not summed, recorded
                           500.50, 16., 1.5    // at the other edge of the window: not summed, recorded
                           500.50, 32., nan    // no mobility
                           500.50, 64., 1.25   // inside the window
                           500.50, 128., 0.25  // outside the window, inside the recorded range
                           500.50, 256., 0.    // at the edge of the recorded range: not recorded
                           500.50, 512., 2.    // at the other edge of the recorded range: not recorded
                           500.50, 1024., 2.5 ]
            let profile = Xic.MobilityProfile()
            let intensity, _, _ = Xic.extractPointWithProfile profile peaks 500.5 0.25 0.5 1.5
            Expect.equal intensity (1. + 2. + 64.) "sum of the window of the point"
            Expect.equal (pairsOf profile 0) [ 0.25, 128.; 0.5, 8.; 1.0, 1.; 1.0, 2.; 1.25, 64.; 1.5, 16. ] "the recorded pairs"
            // a second point without peaks and a third one from a spectrum without ion mobility
            Xic.extractPointWithProfile profile peaks 700. 0.05 0.5 1.5 |> ignore
            let withoutMobility = SortedPeaks.sortPeakArraysByMz [| 500.5 |] [| 5. |] [||]
            let intensityWithout, _, _ = Xic.extractPointWithProfile profile withoutMobility 500.5 0.05 0.5 1.5
            Expect.equal intensityWithout 5. "the closest peak of a spectrum without ion mobility"
            Expect.equal profile.Points 3 "one start per point"
            Expect.isEmpty (pairsOf profile 1) "no pairs"
            Expect.isEmpty (pairsOf profile 2) "a spectrum without ion mobility adds no pairs"
            profile.Clear()
            Expect.equal profile.Points 0 "emptied"

        testCase "summedTrace leaves out the edges of the window, gives the weighted ion mobility and 0 and nan for a point without pairs inside" <| fun _ ->
            let profile = profileOf [ [ 0.75, 1.; 0.8, 2.; 1.0, 4.; 1.2, 8.; 1.25, 16. ]; []; [ 0.5, 3. ] ]
            let intensity, im = Xic.summedTrace profile 1.0 0.25
            Expect.equal intensity.Length 3 "one value per point"
            Expect.floatClose Accuracy.high intensity.[0] (2. + 4. + 8.) "edges exclusive"
            Expect.floatClose tolerance im.[0] ((0.8 * 2. + 1.0 * 4. + 1.2 * 8.) / 14.) "weighted ion mobility"
            Expect.equal intensity.[1] 0. "empty point"
            Expect.isTrue (Double.IsNaN im.[1]) "no ion mobility"
            Expect.equal intensity.[2] 0. "no pair inside the window"
            Expect.isTrue (Double.IsNaN im.[2]) "no ion mobility"
            let wider, _ = Xic.summedTrace profile 1.0 0.26
            Expect.floatClose Accuracy.high wider.[0] (1. + 2. + 4. + 8. + 16.) "a wider window takes them"

        testCase "mobilityCentre of one species is its intensity weighted mean" <| fun _ ->
            let species = [ for i in 0 .. 4 -> [ 1.10 + 0.002 * float i, 100. + float i ] ]
            let profile = profileOf species
            let pairs = species |> List.collect id
            let mean = (pairs |> List.sumBy (fun (im, w) -> im * w)) / (pairs |> List.sumBy snd)
            Expect.floatClose tolerance (Xic.mobilityCentre profile [| 0 .. 4 |] 1.12 0.025) mean "the weighted mean"
            Expect.floatClose tolerance (Xic.mobilityCentre profile [| 0 .. 4 |] 1.10 0.025) mean "the weighted mean from another expected mobility"
            Expect.floatClose tolerance (Xic.mobilityCentre profile [| 2 |] 1.12 0.025) (1.104) "only the given points"

        testCase "mobilityCentre does not move to a stronger species outside the half-width around the expected mobility" <| fun _ ->
            let point = [ 1.000, 300.; 1.004, 200.; 1.045, 1000.; 1.047, 800. ]
            let profile = profileOf [ point; point ]
            let centre = Xic.mobilityCentre profile [| 0; 1 |] 1.002 0.025
            Expect.floatClose tolerance centre ((1.000 * 300. + 1.004 * 200.) / 500.) "the species at the expected mobility"
            // with the other species inside the half-width, both are averaged
            let both = Xic.mobilityCentre profile [| 0; 1 |] 1.02 0.04
            Expect.floatClose tolerance both ((1.000 * 300. + 1.004 * 200. + 1.045 * 1000. + 1.047 * 800.) / 2300.) "everything inside the half-width"

        testCase "mobilityCentre is one step: a species beyond the half-width of the expected mobility does not pull the centre further" <| fun _ ->
            // A at e (weight 1), B at e + 0.8 v (3), C at e + 1.5 v (10): only A and B are within v of e. A second step from their mean
            // 0.6 v would take C (0.9 v away) in.
            let e, v = 1.0, 0.04
            let profile = profileOf [ [ e, 1.; e + 0.8 * v, 3.; e + 1.5 * v, 10. ] ]
            Expect.floatClose { absolute = 1e-6; relative = 0. } (Xic.mobilityCentre profile [| 0 |] e v) (e + 0.6 * v) "the mean of A and B"

        testCase "mobilityCentre is nan when nothing lies inside the half-width or the points are empty" <| fun _ ->
            let profile = profileOf [ [ 1.2, 100. ]; [] ]
            Expect.isTrue (Double.IsNaN (Xic.mobilityCentre profile [| 0; 1 |] 1.0 0.025)) "nothing near the expected mobility"
            Expect.isTrue (Double.IsNaN (Xic.mobilityCentre profile [| 1 |] 1.2 0.025)) "an empty point"
            Expect.isTrue (Double.IsNaN (Xic.mobilityCentre profile [||] 1.2 0.025)) "no points"
            Expect.isTrue (Double.IsNaN (Xic.mobilityCentre (profileOf [ [ 1.0, 100. ] ]) [| 0 |] 1.25 0.25)) "a pair at the edge is outside"

        testCase "tryMobilityCentre takes the centre from the points around the scan time and is None without a pair" <| fun _ ->
            let rt = [| 0.; 1.; 2.; 3.; 4. |]
            let intensity = [| 1.; 5.; 100.; 5.; 1. |]
            let profile = profileOf [ [ 1.0, 1. ]; [ 1.0, 5. ]; [ 1.01, 60.; 1.03, 40.; 1.2, 500. ]; [ 1.0, 5. ]; [ 1.0, 1. ] ]
            match Xic.tryMobilityCentre profile rt intensity 0.5 1.0 0.05 with
            | Some centre -> Expect.floatClose tolerance centre ((1.01 * 60. + 1.03 * 40.) / 100.) "the peak the climb from the scan time ends at"
            | None -> failtest "a centre"
            Expect.isNone (Xic.tryMobilityCentre profile rt intensity 2. 2.0 0.05) "no pair near the expected mobility"
            Expect.isNone (Xic.tryMobilityCentre profile rt (Array.zeroCreate 5) 2. 1.0 0.05) "no centre points"

        testCase "centrePoints climbs from the scan time to the local apex and takes the run above half of it" <| fun _ ->
            let rt = Array.init 12 float
            let points y scanTime = Xic.centrePoints rt y scanTime
            let peak = [| 1.; 2.; 4.; 8.; 20.; 40.; 100.; 60.; 30.; 5.; 2.; 1. |]
            // start on the left flank, on the right flank and at the apex: the same run (at least 50, up to the first point below half)
            Expect.equal (points peak 1.2) [| 6; 7 |] "from the left flank"
            Expect.equal (points peak 9.) [| 6; 7 |] "from the right flank"
            Expect.equal (points peak 6.) [| 6; 7 |] "from the apex"
            // the run ends at the first point below half, a later point above half does not belong to it
            let twoTops = [| 1.; 100.; 60.; 10.; 70.; 5.; 1.; 1.; 1.; 1.; 1.; 1. |]
            Expect.equal (points twoTops 1.) [| 1; 2 |] "stops at the first point below half"
            // a drop-out at the start: the climb goes to the higher neighbour
            let dropOut = [| 1.; 5.; 30.; 0.; 80.; 100.; 40.; 5.; 1.; 1.; 1.; 1. |]
            Expect.equal (points dropOut 3.) [| 4; 5 |] "a zero start climbs to the higher neighbour"
            // equal neighbours end the climb
            let plateau = [| 0.; 0.; 10.; 10.; 10.; 0.; 0.; 0.; 0.; 0.; 0.; 0. |]
            Expect.equal (points plateau 3.) [| 2; 3; 4 |] "a plateau ends the climb at the start"
            // equal higher neighbours: the climb goes to the lower position and the other peak is not part of the run
            let equalNeighbours = [| 0.; 0.; 10.; 4.; 10.; 0.; 0.; 0.; 0.; 0.; 0.; 0. |]
            Expect.equal (points equalNeighbours 3.) [| 2 |] "equal higher neighbours: the lower position"
            // no signal
            Expect.isEmpty (points (Array.zeroCreate 12) 5.) "all zero"
            Expect.isEmpty (Xic.centrePoints [||] [||] 1.) "empty trace"
            // the first and the last point of the trace
            Expect.equal (points [| 50.; 20.; 5.; 1.; 1.; 1.; 1.; 1.; 1.; 1.; 1.; 1. |] 0.) [| 0 |] "start of the trace"
            Expect.equal (points [| 1.; 1.; 1.; 1.; 1.; 1.; 1.; 1.; 1.; 5.; 20.; 50. |] 11.) [| 11 |] "end of the trace"

        testCase "the XIC of the target of the PSM based quantification is summed around the ion mobility centre and leaves out a foreign species" <| fun _ ->
            // 100 spectra, the ion (500.0) elutes at point 50 at the mobilities 1.000 and 1.004, a species at 1.045 elutes with it
            let spectra = Array.init 100 (fun i -> spectrum (species i 3000. 50. 6. [ 1.000, 0.6; 1.004, 0.4 ] @ species i 1500. 50. 6. [ 1.045, 1.0 ]))
            let index = spectra |> Array.mapi (fun i _ -> RtIndexEntry(float i * spacing, string i)) |> MzIOArray.ToMzIOArray
            let read (id: string) = spectra.[int id]
            let scanTime = 50. * spacing
            let getXicAtMobility = PSMBasedQuantificationTIMs.initGetProcessedXIC None read index 1.5 0.02 0.025
            let getTargetXic = PSMBasedQuantificationTIMs.initGetTargetXic None read index 1.5 0.02 0.025
            let ionMobility = 1.000 * 0.6 + 1.004 * 0.4
            // the expected mobility is a little off: the signal is found around it, the foreign species is not in the sum
            for expected in [ 1.000; 1.008; 1.015 ] do
                match getTargetXic scanTime 500.0 expected with
                | Some centre, _ -> Expect.floatClose tolerance4 centre ionMobility (sprintf "centre for the expected mobility %g" expected)
                | None, _ -> failtestf "no centre for the expected mobility %g" expected
            // the XIC is the one of a window centred on the ion
            let centre, (rt, corrected, raw) = getTargetXic scanTime 500.0 1.015
            let directRt, directCorrected, directRaw = getXicAtMobility scanTime 500.0 centre.Value
            Expect.equal rt directRt "scan times"
            Expect.equal raw directRaw "intensities of a window centred on the ion"
            Expect.equal corrected directCorrected "corrected intensities"
            let _, _, offRaw = getXicAtMobility scanTime 500.0 1.030
            Expect.isTrue (Array.sum offRaw < 0.8 * Array.sum raw) "a window centred on an expected mobility off the ion misses signal"
            Expect.isTrue (Array.max raw < 3001.) "the foreign species is left out"
            // a scan time next to the peak climbs to it
            match getTargetXic (scanTime + 0.5) 500.0 1.0 with
            | Some next, _ -> Expect.floatClose tolerance4 next ionMobility "the centre of the peak next to the scan time"
            | None, _ -> failtest "no centre next to the peak"
            // an expected mobility that is not a number
            let notANumber, xicAtNan = getTargetXic scanTime 500.0 nan
            Expect.isNone notANumber "no expected mobility"
            Expect.equal xicAtNan (getXicAtMobility scanTime 500.0 nan) "the XIC at the mobility nan"

        testCase "the XIC of the target holds a shoulder of the ion that lies beyond the window around the expected mobility" <| fun _ ->
            // the ion at 1.000 / 1.004 has a shoulder at 1.020; the expected mobility 0.99 gives the window ending at 1.015
            let spectra = Array.init 100 (fun i -> spectrum (species i 3000. 50. 6. [ 1.000, 0.5; 1.004, 0.3; 1.020, 0.2 ]))
            let index = spectra |> Array.mapi (fun i _ -> RtIndexEntry(float i * spacing, string i)) |> MzIOArray.ToMzIOArray
            let read (id: string) = spectra.[int id]
            let scanTime = 50. * spacing
            let getXicAtMobility = PSMBasedQuantificationTIMs.initGetProcessedXIC None read index 1.5 0.02 0.025
            let getTargetXic = PSMBasedQuantificationTIMs.initGetTargetXic None read index 1.5 0.02 0.025
            let centre, (_, _, raw) = getTargetXic scanTime 500.0 0.99
            Expect.floatClose tolerance4 centre.Value ((1.000 * 0.5 + 1.004 * 0.3) / 0.8) "the centre of the ion without the shoulder"
            let _, _, plain = getXicAtMobility scanTime 500.0 0.99
            Expect.floatClose tolerance (Array.max plain) 2400. "the window around the expected mobility cuts the shoulder"
            Expect.floatClose tolerance (Array.max raw) 3000. "the XIC holds the whole ion"
            Expect.isTrue (Array.sum raw > Array.sum plain) "more than the plain window"

        testCase "the XIC of the target without a centre is the XIC at the expected mobility, and the profile is emptied for the next ion" <| fun _ ->
            // the ion (500.0) elutes at point 50 at 1.000, a second ion (700.0) at point 30 at 1.100, a species at 500.0 at 1.040 elutes with the ion
            let spectra =
                Array.init 100 (fun i ->
                    spectrum (species i 3000. 50. 6. [ 1.000, 1.0 ] @ speciesAt 700.0 i 2000. 30. 6. [ 1.100, 1.0 ] @ speciesAt 600.0 i 1000. 50. 6. [ 1.040, 1.0 ]))
            let index = spectra |> Array.mapi (fun i _ -> RtIndexEntry(float i * spacing, string i)) |> MzIOArray.ToMzIOArray
            let read (id: string) = spectra.[int id]
            let scanTime = 50. * spacing
            let getXicAtMobility = PSMBasedQuantificationTIMs.initGetProcessedXIC None read index 1.5 0.02 0.025
            let create () = PSMBasedQuantificationTIMs.initGetTargetXic None read index 1.5 0.02 0.025
            // no signal at the m/z, and signal beyond the half-width of the expected mobility: no centre and the XIC at the expected mobility
            let getTargetXic = create ()
            for mz, expected in [ 650.0, 1.0; 600.0, 1.0 ] do
                let centre, xic = getTargetXic scanTime mz expected
                Expect.isNone centre (sprintf "no centre at the m/z %g" mz)
                Expect.equal xic (getXicAtMobility scanTime mz expected) (sprintf "the XIC at the expected mobility at the m/z %g" mz)
            // another ion first, then the ion: the same XIC as for the ion alone
            let alone = (create ()) scanTime 500.0 1.01
            let shared = create ()
            let other = shared (30. * spacing) 700.0 1.09
            Expect.isSome (fst other) "the other ion has a centre"
            Expect.equal (shared scanTime 500.0 1.01) alone "the ion after another ion"
            Expect.equal (shared scanTime 500.0 1.01) alone "the ion again"
            // signal of 0.1, which float32 does not hold exactly, at other scan times and none around the start point: no centre, and the trace
            // is the one that was extracted, not a sum of recorded float32 pairs
            let faint = Array.init 100 (fun i -> spectrum (if i >= 20 && i <= 24 then [ 500.0, 0.1, 1.0 ] else []))
            let faintIndex = faint |> Array.mapi (fun i _ -> RtIndexEntry(float i * spacing, string i)) |> MzIOArray.ToMzIOArray
            let readFaint (id: string) = faint.[int id]
            let faintCentre, (_, _, faintRaw) = PSMBasedQuantificationTIMs.initGetTargetXic None readFaint faintIndex 1.5 0.02 0.025 scanTime 500.0 1.0
            Expect.isNone faintCentre "no centre"
            Expect.equal (Array.max faintRaw) 0.1 "the intensity as it was extracted"
            Expect.equal faintRaw (let _, _, y = PSMBasedQuantificationTIMs.initGetProcessedXIC None readFaint faintIndex 1.5 0.02 0.025 scanTime 500.0 1.0 in y) "the trace at the expected mobility"

        testCase "the sweep of a request that records profiles gives the raw XIC of the request that sums, and hands over one profile per trace" <| fun _ ->
            let spectra = Array.init 60 (fun i -> spectrum (species i 3000. 30. 3. [ 1.000, 0.6; 1.004, 0.4 ] @ speciesAt 600.0 i 500. 30. 3. [ 1.010, 1.0 ]))
            let request = requestFor [| 500.0; 600.0 |] 1.0 0.025
            let raw, profiles = sweepOne spectra request
            let plain, noProfiles = sweepOne spectra { request with Point = XicSweep.Sum }
            Expect.equal raw.Rt plain.Rt "scan times"
            Expect.equal (bitsOf raw.Intensity) (bitsOf plain.Intensity) "intensities"
            Expect.equal (bitsOf raw.ObservedMz) (bitsOf plain.ObservedMz) "observed m/z"
            Expect.equal (bitsOf raw.IonMobility) (bitsOf plain.IonMobility) "ion mobilities"
            Expect.isNone noProfiles "no profile requested"
            Expect.equal (profiles |> Option.map (Array.map (fun p -> p.Points))) (Some [| 60; 60 |]) "a profile with a start per point for every trace"

        testCase "the sweep takes the closest peak for a request that asks for it" <| fun _ ->
            let spectra = [| spectrum [ 500.0, 10., 1.0; 500.001, 30., 1.0; 500.015, 5., 1.0 ] |]
            let request = { requestFor [| 500.0 |] 1.0 0.025 with Point = XicSweep.ClosestPeak }
            let raw, profiles = sweepOne spectra request
            Expect.isNone profiles "no profile requested"
            // the sum of the m/z group closest to the lock m/z, not the sum of all peaks of the window
            Expect.equal raw.Intensity.[0].[0] 10. "the closest peak"
            let summed, _ = sweepOne spectra { request with Point = XicSweep.Sum }
            Expect.equal summed.Intensity.[0].[0] 45. "the sum"

        testCase "the summed XIC has every trace summed at one centre, taken from the peak at the predicted scan time" <| fun _ ->
            // The ion elutes at point 40 at 1.000 / 1.004 with a shoulder at 1.020. A stronger species at 1.010 elutes at point 15. The
            // expected mobility 0.99 gives the window around it ending at 1.015, so the shoulder is outside it. The isotope trace 500.5
            // and the partner (600.0) are at the mobilities of the ion; the isotope trace has a species at 1.035 of its own.
            let spectra =
                Array.init 80 (fun i ->
                    spectrum (
                        species i 3000. 40. 3. [ 1.000, 0.6; 1.004, 0.4 ]
                        @ species i 600. 40. 3. [ 1.020, 1.0 ]
                        @ species i 9000. 15. 3. [ 1.010, 1.0 ]
                        @ speciesAt 500.5 i 1500. 40. 3. [ 1.002, 1.0 ]
                        @ speciesAt 500.5 i 700. 40. 3. [ 1.035, 1.0 ]
                        @ speciesAt 600.0 i 800. 40. 3. [ 1.003, 1.0 ]))
            let request = requestFor [| 500.0; 500.5; 600.0 |] 0.99 0.025
            let raw, profiles = sweepOne spectra request
            Expect.isTrue (Array.max raw.Intensity.[0] > 8000.) "the stronger species is the highest point of the trace"
            Expect.equal raw.Intensity.[0].[40] 3000. "the window around the expected mobility cuts the shoulder"
            let summed, centre = XicSweep.sumAroundMobilityCentre 0.99 0.025 (40. * spacing) raw profiles.Value
            Expect.floatClose tolerance4 centre.Value (1.000 * 0.6 + 1.004 * 0.4) "the centre of the peak at the scan time"
            Expect.equal summed.Intensity.[0].[40] 3600. "the window around the centre holds the shoulder"
            // every trace is the sum of the pairs within the half-width around the centre, with their weighted ion mobility
            for trace in 0 .. 2 do
                for point in 0 .. raw.Rt.Length - 1 do
                    let intensity, _, im = Xic.extractPoint spectra.[point] request.LockMz.[trace] 0.02 (centre.Value - 0.025) (centre.Value + 0.025)
                    Expect.equal summed.Intensity.[trace].[point] intensity (sprintf "intensity of trace %i at point %i" trace point)
                    if intensity > 0. then Expect.floatClose tolerance summed.IonMobility.[trace].[point] im "ion mobility"
                    else Expect.isTrue (Double.IsNaN summed.IonMobility.[trace].[point]) "no ion mobility without signal"
            Expect.equal summed.Rt raw.Rt "scan times"
            Expect.equal summed.Intensity.[1].[40] 1500. "the species of the isotope trace beyond the window is left out"
            // a predicted scan time at the stronger species takes the centre from it
            let _, other = XicSweep.sumAroundMobilityCentre 0.99 0.025 (15. * spacing) raw profiles.Value
            Expect.floatClose tolerance4 other.Value 1.010 "the centre of the peak at the other scan time"

        testCase "the observed m/z of a point is that of the sweep where the summed intensity is positive and nan where it is 0" <| fun _ ->
            // the target is at 1.015, so its centre is 1.015 and the window around it is 0.99 .. 1.04; the partner (600.0) is at 0.980,
            // inside the window around the expected mobility 1.0 (0.975 .. 1.025) and outside the one around the centre
            let spectra = Array.init 60 (fun i -> spectrum (species i 3000. 30. 3. [ 1.015, 1.0 ] @ speciesAt 600.0 i 800. 30. 3. [ 0.980, 1.0 ]))
            let raw, profiles = sweepOne spectra (requestFor [| 500.0; 600.0 |] 1.0 0.025)
            Expect.isTrue (raw.ObservedMz.[1] |> Array.exists (fun m -> not (Double.IsNaN m))) "the sweep has an m/z for the partner"
            let summed, centre = XicSweep.sumAroundMobilityCentre 1.0 0.025 (30. * spacing) raw profiles.Value
            Expect.floatClose tolerance centre.Value 1.015 "centre of the target"
            Expect.isTrue (summed.Intensity.[1] |> Array.forall (fun v -> v = 0.)) "nothing of the partner inside the window around the centre"
            Expect.isTrue (summed.ObservedMz.[1] |> Array.forall Double.IsNaN) "no m/z where the summed intensity is 0"
            for point in 0 .. raw.Rt.Length - 1 do
                let expected = if summed.Intensity.[0].[point] > 0. then bits raw.ObservedMz.[0].[point] else bits nan
                Expect.equal (bits summed.ObservedMz.[0].[point]) expected "the m/z of the sweep where the summed intensity is positive"

        testCase "a peak off the expected ion mobility is not cut by the summed XIC" <| fun _ ->
            // the ion at 1.020 / 1.024 with a shoulder at 1.040 lies off the expected mobility 1.0: the window around the expected mobility
            // ends at 1.025, the one around the centre takes the shoulder
            let spectra = Array.init 80 (fun i -> spectrum (species i 3000. 40. 3. [ 1.020, 0.5; 1.024, 0.3; 1.040, 0.2 ]))
            let request = requestFor [| 500.0 |] 1.0 0.025
            let raw, profiles = sweepOne spectra request
            let summed, centre = XicSweep.sumAroundMobilityCentre 1.0 0.025 (40. * spacing) raw profiles.Value
            Expect.floatClose tolerance4 centre.Value (1.020 * 5. / 8. + 1.024 * 3. / 8.) "the centre"
            Expect.floatClose tolerance (Array.max summed.Intensity.[0]) 3000. "the whole peak"
            Expect.isTrue (Array.max raw.Intensity.[0] < 2900.) "the window around the expected mobility cuts the shoulder"

        testCase "the XIC stays as the sweep gave it without a signal around the expected ion mobility" <| fun _ ->
            // nothing at the m/z, and a species that lies within the recorded range but farther than the half-width from the expected mobility
            let nothing = Array.init 60 (fun i -> spectrum (species i 3000. 30. 3. [ 1.000, 1.0 ] |> List.map (fun (_, y, im) -> 700.0, y, im)))
            let beside = Array.init 60 (fun i -> spectrum (species i 3000. 30. 3. [ 1.040, 1.0 ]))
            let noMobility = Array.init 60 (fun i -> spectrum (species i 3000. 30. 3. [ nan, 1.0 ]))
            for spectra in [ nothing; beside; noMobility ] do
                let raw, profiles = sweepOne spectra (requestFor [| 500.0 |] 1.0 0.025)
                let summed, centre = XicSweep.sumAroundMobilityCentre 1.0 0.025 (30. * spacing) raw profiles.Value
                Expect.isNone centre "no centre from the signal"
                Expect.equal summed.Rt raw.Rt "scan times"
                Expect.equal (bitsOf summed.Intensity) (bitsOf raw.Intensity) "intensities"
                Expect.equal (bitsOf summed.ObservedMz) (bitsOf raw.ObservedMz) "observed m/z"
                Expect.equal (bitsOf summed.IonMobility) (bitsOf raw.IonMobility) "ion mobilities"
            // no points
            let empty, emptyProfiles = sweepOne [||] (requestFor [| 500.0 |] 1.0 0.025)
            let summedEmpty, emptyCentre = XicSweep.sumAroundMobilityCentre 1.0 0.025 1. empty emptyProfiles.Value
            Expect.isNone emptyCentre "no points"
            Expect.isEmpty summedEmpty.Rt "no points"

        testCase "the XIC stays as extracted where only a faint signal at another scan time lies in the window" <| fun _ ->
            // 0.1 is not exactly representable as float32, and there is nothing at and around the scan time: no centre
            let spectra = Array.init 60 (fun i -> spectrum (if i >= 5 && i <= 8 then [ 500.0, 0.1, 1.0 ] else []))
            let raw, profiles = sweepOne spectra (requestFor [| 500.0 |] 1.0 0.025)
            Expect.equal (Array.max raw.Intensity.[0]) 0.1 "the intensity as it was extracted"
            let summed, centre = XicSweep.sumAroundMobilityCentre 1.0 0.025 (30. * spacing) raw profiles.Value
            Expect.isNone centre "no centre"
            Expect.equal (bitsOf summed.Intensity) (bitsOf raw.Intensity) "the intensities of the sweep, bit for bit"

        testCase "the labeled partner is summed at the centre of the target" <| fun _ ->
            // the target (500.0) is at 1.000, the partner (600.0) at 1.020 with a second species at 1.030. The expected mobility is 1.005:
            // the centre of the target is 1.000, so the window around it takes the partner at 1.020 but not the species at 1.030, which
            // the window around the centre of the partner (1.025) would take
            let spectra =
                Array.init 60 (fun i ->
                    spectrum (species i 3000. 30. 3. [ 1.000, 1.0 ] @ speciesAt 600.0 i 800. 30. 3. [ 1.020, 1.0 ] @ speciesAt 600.0 i 800. 30. 3. [ 1.030, 1.0 ]))
            let raw, profiles = sweepOne spectra (requestFor [| 500.0; 600.0 |] 1.005 0.025)
            let summed, centre = XicSweep.sumAroundMobilityCentre 1.005 0.025 (30. * spacing) raw profiles.Value
            Expect.floatClose tolerance centre.Value 1.000 "centre of the target"
            Expect.floatClose tolerance (Array.max summed.Intensity.[1]) 800. "the partner at 1.020 only"
            Expect.floatClose tolerance4 (summed.IonMobility.[1] |> Array.filter (fun v -> not (Double.IsNaN v)) |> Array.max) 1.020 "its ion mobility"
    

        testCase "the tally counts the ions with and without a centre" <| fun _ ->
            let tally = Xic.MobilityCentreTally()
            Expect.equal tally.Count 0 "nothing counted"
            tally.Add (1.00, Some 1.01)
            tally.Add (1.00, Some 0.98)
            tally.Add (1.00, None)
            Expect.equal tally.Count 3 "three ions"
            Expect.equal (tally.Summary 0.025) "2 ions summed over +-0.0250 around the centre of their signal, 1 keep the expected ion mobility; centre minus expected ion mobility: median -0.00500, 5th percentile -0.01850, 95th percentile 0.00850" "summary"
    ]
