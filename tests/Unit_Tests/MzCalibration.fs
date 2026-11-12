/// Tests of the m/z correction the quantification tools share.
///
/// The reference test fits subsamples of two runs (every k-th deduplicated point with k = ceil(n / 10000), in tests/data/MzCalibration/in)
/// and compares them with a Python implementation of the same steps run on the same files. Probes are the 1, 25, 50, 75 and 99 % quantiles of scan time and m/z of the points.
module ProteomIQon.UnitTestsMzCalibration

open System
open System.IO
open Expecto
open FSharp.Stats
open ProteomIQon

let private ppmTolerance = { absolute = 1e-6; relative = 0. }

let private observationOf id charge globalMod time theoreticalMz (offsetPpm: float) : MzCalibration.Observation =
    {
        ModSequenceID = id
        Charge        = charge
        GlobalMod     = globalMod
        ScanTime      = time
        TheoreticalMz = theoreticalMz
        MeasuredMz    = theoreticalMz * (1. + offsetPpm * 1e-6)
    }

let private observation id charge = observationOf id charge 0

/// Offset in ppm of the correction (Da) at a scan time and m/z.
let private correctionPpm (c: MzCalibration.Calibration) t mz = c.CorrectionAt t mz / mz * 1e6

/// The half window of the old definition: mzWindowSdFactor x the sample sd of the offsets (Da) of all observations inside Tukey fences
/// of 3 IQR and below 0.1 Da, the fallback window if it is not finite or not positive.
let private expectedHalfWindow (observations: MzCalibration.Observation []) =
    let diffs = observations |> Array.map (fun o -> o.MeasuredMz - o.TheoreticalMz)
    let lower, upper =
        match (diffs |> FSharp.Stats.Signal.Outliers.tukey 3.).TryToTuple with
        | Some fences -> fences
        | None -> nan, nan
    let inside = diffs |> Array.filter (fun d -> d >= lower && d <= upper && d < 0.1)
    let sd = if inside.Length < 2 then nan else inside |> Seq.stDev
    if Double.IsFinite sd && sd > 0. then MzCalibration.mzWindowSdFactor * sd else MzCalibration.fallbackMzWindow

let private readPoints (name: string) =
    File.ReadAllLines (Path.Combine (Pipeline.repoRoot, "tests", "data", "MzCalibration", "in", name))
    |> Array.skip 1
    |> Array.map (fun line ->
        let f = line.Split '\t'
        observationOf (int f.[0]) (int f.[1]) (int f.[2]) (float f.[3]) (float f.[4]) (float f.[5])
        )

type private Reference =
    {
        File   : string
        Spec   : int * int
        N1     : int
        N2     : int
        /// Scan time, theoretical m/z and the offset in ppm of the prototype.
        Probes : (float * float * float) []
    }

let private references =
    [
        {
            File   = "D1_UPS1_04_subsample.tsv"
            Spec   = (32, 6)
            N1     = 6546
            N2     = 6488
            Probes =
                [|
                    20.322601239999997, 354.64093099353926, 3.9439695846028924
                    65.200102, 589.5514619807944, 6.1498059271574395
                    103.10856, 745.8998024726969, 7.416164808458614
                    147.70943, 920.4654148965219, 8.41301646776769
                    245.0406688, 1371.1055321344809, 10.453732237722898
                |]
        }
        {
            File   = "D4_Slot1-77_subsample.tsv"
            Spec   = (4, 3)
            N1     = 8467
            N2     = 8406
            Probes =
                [|
                    6.866495012280001, 334.90142364027264, 0.26115823716699416
                    27.28706268, 519.737291214857, 4.575809095556495
                    44.562100053333324, 643.368798372772, 6.009986519632387
                    60.94048880666667, 798.0223629685787, 5.700987693857758
                    84.844268545, 1187.5012418937195, 6.510067854731114
                |]
        }
    ]

/// 3000 ions: offset 2 + 0.05 t + 0.004 (m/z - 400) ppm and Gaussian noise of 1 ppm.
let private drifting =
    let rnd = Random 42
    let gaussian () = sqrt (-2. * log (1. - rnd.NextDouble ())) * cos (2. * Math.PI * rnd.NextDouble ())
    Array.init 3000 (fun i ->
        let t = 60. * rnd.NextDouble ()
        let mz = 400. + 800. * rnd.NextDouble ()
        observation i 2 t mz (2. + 0.05 * t + 0.004 * (mz - 400.) + gaussian ())
        )

[<Tests>]
let mzCalibrationTests =
    testList "MzCalibration" [
        testList "reference values of the prototype" [
            for r in references ->
                testCase r.File <| fun () ->
                    let points = readPoints r.File
                    let c = MzCalibration.fit points
                    Expect.equal c.Spec r.Spec "model"
                    Expect.equal c.PointsFence1 r.N1 "points after the first fences"
                    Expect.equal c.PointsFence2 r.N2 "points after the second fences"
                    Expect.floatClose Accuracy.high c.HalfWindow (expectedHalfWindow points) "half window"
                    for t, mz, expected in r.Probes do
                        Expect.floatClose ppmTolerance (correctionPpm c t mz) expected (sprintf "correction at %f min, m/z %f" t mz)
        ]
        testCase "the correction follows the drift in scan time and m/z" <| fun () ->
            let c = MzCalibration.fit drifting
            for t, mz in [ 5., 450.; 30., 800.; 55., 1100. ] do
                let truth = 2. + 0.05 * t + 0.004 * (mz - 400.)
                Expect.floatClose { absolute = 0.3; relative = 0. } (correctionPpm c t mz) truth (sprintf "correction at %f min, m/z %f" t mz)
            Expect.equal c.PointsUnique 3000 "unique points"
        testCase "the window is one constant, mzWindowSdFactor x the sd of the offsets in Da" <| fun () ->
            let c = MzCalibration.fit drifting
            Expect.floatClose Accuracy.high c.HalfWindow (expectedHalfWindow drifting) "half window"
            Expect.isGreaterThan c.HalfWindow 0. "positive"
        testCase "residualHalfWindow is mzWindowSdFactor x the sd of the noise around the correction" <| fun () ->
            let c = MzCalibration.fit drifting
            let noiseSdDa = 1e-6 * sqrt (drifting |> Array.averageBy (fun o -> o.TheoreticalMz ** 2.))
            let w = MzCalibration.residualHalfWindow c
            Expect.floatClose { absolute = 0.; relative = 0.1 } w (MzCalibration.mzWindowSdFactor * noiseSdDa) "4 x noise sd"
            Expect.isLessThan w c.HalfWindow "smaller than the window of the offsets"
        testCase "residualHalfWindow with fewer than 3 points is the half window of the calibration" <| fun () ->
            let c = MzCalibration.fit [| observation 1 2 10. 500. 3.; observation 2 2 11. 600. 3. |]
            Expect.equal (MzCalibration.residualHalfWindow c) c.HalfWindow "half window"
        testCase "the fit is deterministic" <| fun () ->
            let a = MzCalibration.fit drifting
            let b = MzCalibration.fit (Array.rev drifting)
            Expect.equal a.Spec b.Spec "model"
            Expect.equal (a.CorrectionAt 20. 700.) (b.CorrectionAt 20. 700.) "correction, independent of the order of the points"
        testCase "repeated measurements of an ion at the same m/z are one point" <| fun () ->
            let repeated = Array.append drifting (drifting |> Array.take 100 |> Array.map (fun o -> { o with ScanTime = o.ScanTime + 0.01 }))
            let c = MzCalibration.fit repeated
            Expect.equal c.PointsUnique 3000 "unique points"
        testCase "the half window is computed from all observations, before the deduplication" <| fun () ->
            let repeated = Array.append drifting (drifting |> Array.take 500 |> Array.map (fun o -> { o with ScanTime = o.ScanTime + 0.01; MeasuredMz = o.MeasuredMz + 0.004 }))
            let c = MzCalibration.fit repeated
            Expect.floatClose Accuracy.high c.HalfWindow (expectedHalfWindow repeated) "half window of all observations"
            Expect.notEqual c.HalfWindow (MzCalibration.fit drifting).HalfWindow "the repeated measurements count"
        testCase "the correction is constant beyond the outer knots" <| fun () ->
            let c = MzCalibration.fit drifting
            let tMin = drifting |> Array.map (fun o -> o.ScanTime) |> Array.min
            let tMax = drifting |> Array.map (fun o -> o.ScanTime) |> Array.max
            Expect.equal (correctionPpm c (tMin - 1.) 700.) (correctionPpm c (tMin - 100.) 700.) "below the lowest scan time"
            Expect.equal (correctionPpm c (tMax + 1.) 700.) (correctionPpm c (tMax + 100.) 700.) "above the highest scan time"
            Expect.equal (correctionPpm c 30. 100.) (correctionPpm c 30. 50.) "below the lowest m/z"
            Expect.equal (correctionPpm c 30. 3000.) (correctionPpm c 30. 5000.) "above the highest m/z"
        testCase "fewer than 3 points: no correction and the fallback window" <| fun () ->
            let c = MzCalibration.fit [| observation 1 2 10. 500. 3. |]
            Expect.equal (c.CorrectionAt 10. 500.) 0. "no correction"
            Expect.equal c.HalfWindow MzCalibration.fallbackMzWindow "fallback window"
            Expect.equal c.Spec (0, 0) "no model"
        testCase "no points: no correction and the fallback window" <| fun () ->
            let c = MzCalibration.fit [||]
            Expect.equal (c.CorrectionAt 10. 500.) 0. "no correction"
            Expect.equal c.HalfWindow MzCalibration.fallbackMzWindow "fallback window"
        testCase "only isotope assignment errors: no correction and the fallback window" <| fun () ->
            // 1 Da at charge 2 is 2 x the limit
            let c = MzCalibration.fit (Array.init 100 (fun i -> observation i 2 (float i) 500. (1. / 500. * 1e6)))
            Expect.equal c.PointsUnique 100 "unique points"
            Expect.equal c.PointsNoIsotopeError 0 "no point is left"
            Expect.equal (c.CorrectionAt 5. 500.) 0. "no correction"
            Expect.equal c.HalfWindow MzCalibration.fallbackMzWindow "fallback window"
        testCase "isotope assignment errors are left out of the fit" <| fun () ->
            let errors = Array.init 300 (fun i -> observation (10000 + i) 2 (60. * float i / 300.) 600. (1. / 600. * 1e6))
            let c = MzCalibration.fit (Array.append drifting errors)
            Expect.equal c.PointsNoIsotopeError 3000 "points without isotope error"
            Expect.floatClose { absolute = 0.3; relative = 0. } (correctionPpm c 30. 800.) (2. + 1.5 + 1.6) "correction"
        testCase "fewer than 20 points inside the first fences: the median" <| fun () ->
            let points = Array.init 15 (fun i -> observation i 2 (float i) (400. + 10. * float i) (3. + 0.1 * float (i % 3)))
            let c = MzCalibration.fit points
            Expect.equal c.Spec (1, 1) "constant"
            Expect.floatClose Accuracy.high (correctionPpm c 0. 400.) 3.1 "median"
            Expect.floatClose Accuracy.high (correctionPpm c 14. 540.) 3.1 "median at any scan time and m/z"
            Expect.floatClose Accuracy.high c.HalfWindow (expectedHalfWindow points) "window"
        testCase "constant offsets: the correction is the offset" <| fun () ->
            let c = MzCalibration.fit (Array.init 100 (fun i -> observation i 2 (float i) 500. 2.))
            Expect.floatClose Accuracy.high (correctionPpm c 50. 500.) 2. "correction"
        testCase "non-finite input is ignored" <| fun () ->
            let bad =
                [|
                    { drifting.[0] with ScanTime = nan }
                    { drifting.[1] with MeasuredMz = nan }
                    { drifting.[2] with TheoreticalMz = infinity }
                    { drifting.[3] with MeasuredMz = -infinity }
                |]
            let clean = MzCalibration.fit (drifting |> Array.skip 4)
            let withBad = MzCalibration.fit (Array.append bad (drifting |> Array.skip 4))
            Expect.equal withBad.PointsUnique clean.PointsUnique "unique points"
            Expect.equal withBad.HalfWindow clean.HalfWindow "half window"
            Expect.equal (withBad.CorrectionAt 20. 700.) (clean.CorrectionAt 20. 700.) "correction"
            Expect.equal (clean.CorrectionAt nan 700.) 0. "no correction at a scan time that is nan"
            Expect.equal (clean.CorrectionAt 20. nan) 0. "no correction at an m/z that is nan"
        testCase "only non-finite points: no correction" <| fun () ->
            let c = MzCalibration.fit (Array.init 10 (fun i -> observation i 2 nan 500. 1.))
            Expect.equal c.PointsUnique 0 "unique points"
            Expect.equal (c.CorrectionAt 1. 500.) 0. "no correction"
        testCase "the window: fixed or estimated" <| fun () ->
            let c = MzCalibration.fit drifting
            Expect.equal (MzCalibration.mzWindow (Domain.Window.Fixed 0.01) c) 0.01 "fixed"
            Expect.equal (MzCalibration.mzWindow Domain.Window.Estimate c) c.HalfWindow "estimate"
    ]
