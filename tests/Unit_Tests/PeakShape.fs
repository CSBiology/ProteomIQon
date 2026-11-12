/// Tests of the apex and the width of fitted Gaussian and EMG peaks.
module ProteomIQon.UnitTestsPeakShape

open System
open Expecto
open FSharp.Stats
open ProteomIQon

/// The EMG function with mean 0 and sigma 1 as in BioFSharp.Mz (Quantification, Table.emgModel), with the erfc of FSharp.Stats.
let private emgLibrary (tau: float) (x: float) =
    (1. / tau) * sqrt (Math.PI / 2.)
    * exp (0.5 / (tau * tau) - x / tau)
    * SpecialFunctions.Errorfunction.Erfc ((1. / tau - x) / sqrt 2.)

/// The same function with exp(0.5 / tau^2 - x / tau) * erfc(t) = exp(-x^2 / 2) * erfcx(t), t = (1 / tau - x) / sqrt 2, with an accurate erfcx.
let private emgExact (tau: float) (x: float) =
    (1. / tau) * sqrt (Math.PI / 2.) * exp (-x * x / 2.) * PeakShape.erfcx ((1. / tau - x) / sqrt 2.)

/// Position of the maximum of f on [lo, hi] by golden section search.
let private argMax (f: float -> float) (lo: float) (hi: float) =
    let phi = (sqrt 5. - 1.) / 2.
    let mutable a = lo
    let mutable b = hi
    while b - a > 1e-10 do
        let c = b - phi * (b - a)
        let d = a + phi * (b - a)
        if f c > f d then b <- d else a <- c
    0.5 * (a + b)

[<Tests>]
let tests =
    testList "PeakShape" [
        testList "mode of the EMG" [
            testCase "reference values" <| fun () ->
                let g = PeakShape.emgModeOffset
                Expect.floatClose { absolute = 5e-5; relative = 0. } (g 0.1) 0.0990 "r = 0.1"
                Expect.floatClose { absolute = 5e-5; relative = 0. } (g 1.) 0.6974 "r = 1"
                Expect.floatClose { absolute = 5e-5; relative = 0. } (g 3.) 1.2154 "r = 3"
                Expect.floatClose { absolute = 5e-3; relative = 0. } (g 10.) 1.79 "r = 10"
            testCase "equals the maximum of the EMG function" <| fun () ->
                for r in [0.05; 0.3; 1.; 2.5; 8.; 30.] do
                    let numerical = argMax (emgExact r) -3. 15.
                    Expect.floatClose { absolute = 1e-6; relative = 0. } (PeakShape.emgModeOffset r) numerical (sprintf "r = %g" r)
            testCase "equals the maximum of the EMG function of BioFSharp.Mz where its erfc is accurate" <| fun () ->
                // the erfc of FSharp.Stats has an absolute error of about 1e-7, so the function is off for small r (tail at large arguments)
                for r in [0.3; 1.; 2.5; 8.; 30.] do
                    let numerical = argMax (emgLibrary r) -3. 15.
                    Expect.floatClose { absolute = 2e-3; relative = 0. } (PeakShape.emgModeOffset r) numerical (sprintf "r = %g" r)
            testCase "not usable r gives 0" <| fun () ->
                for r in [0.; -1.; nan; infinity; -infinity] do
                    Expect.equal (PeakShape.emgModeOffset r) 0. (sprintf "r = %g" r)
            testCase "continuous and increasing in r from 1e-3 to 1e3" <| fun () ->
                let rs = [| for i in 0 .. 6000 -> 10. ** (-3. + float i / 1000.) |]
                let g = rs |> Array.map PeakShape.emgModeOffset
                Expect.isTrue (g |> Array.forall (fun v -> not (Double.IsNaN v || Double.IsInfinity v))) "finite"
                for i in 1 .. g.Length - 1 do
                    Expect.isTrue (g.[i] > g.[i - 1]) (sprintf "increasing at r = %g" rs.[i])
                    // a step of 0.23 % in r moves the mode by far less than 0.01 sigma
                    Expect.isTrue (g.[i] - g.[i - 1] < 0.01) (sprintf "no jump at r = %g" rs.[i])
                Expect.isTrue (g.[0] > 0. && g.[0] < 0.001) "mode offset at r = 1e-3 is about r"
            testCase "no overflow for extreme r" <| fun () ->
                for r in [1e-9; 1e-6; 1e6; 1e9] do
                    let v = PeakShape.emgModeOffset r
                    Expect.isTrue (not (Double.IsNaN v || Double.IsInfinity v)) (sprintf "finite, r = %g" r)
                    // about r for small r (up to the rounding of 1 / r), several sigma for large r
                    Expect.isTrue ((if r < 1. then abs v < 1e-5 else v > 3.)) (sprintf "size, r = %g" r)
            testCase "erfcx against known values" <| fun () ->
                let close = Expect.floatClose { absolute = 0.; relative = 1e-10 }
                close (PeakShape.erfcx 0.) 1. "t = 0"
                close (PeakShape.erfcx 1.) 0.427583576155807 "t = 1"
                close (PeakShape.erfcx 2.) 0.255395676310506 "t = 2"
                close (PeakShape.erfcx 10.) 0.0561409927438226 "t = 10"
                close (PeakShape.erfcx (-1.)) 5.00898008076228 "t = -1"
            testCase "erfcx is continuous at the switch point and matches the asymptotic series" <| fun () ->
                Expect.floatClose { absolute = 0.; relative = 1e-12 } (PeakShape.erfcx (Math.BitDecrement 2.)) (PeakShape.erfcx 2.) "t = 2"
                // 1 / (t sqrt pi) (1 - 1 / (2 t^2) + 3 / (4 t^4) - 15 / (8 t^6)) at t = 700 and 50
                for t in [700.; 50.] do
                    let series = 1. / (t * sqrt Math.PI) * (1. - 1. / (2. * t * t) + 3. / (4. * t ** 4.) - 15. / (8. * t ** 6.))
                    Expect.floatClose { absolute = 0.; relative = 1e-10 } (PeakShape.erfcx t) series (sprintf "t = %g" t)
        ]
        testList "apex and width" [
            testCase "Gaussian parameters give mean and sigma" <| fun () ->
                let p = [| 5.; 31.2; 0.07 |]
                Expect.equal (PeakShape.apex p) 31.2 "apex"
                Expect.equal (PeakShape.width p) 0.07 "width"
            testCase "EMG: apex right of the mean, width from sigma and tau" <| fun () ->
                let p = [| 5.; 31.2; 0.1; 0.1 |]
                Expect.floatClose { absolute = 1e-4; relative = 0. } (PeakShape.apex p) (31.2 + 0.1 * 0.6974) "apex"
                Expect.floatClose { absolute = 1e-12; relative = 0. } (PeakShape.width p) (sqrt 0.02) "width"
                Expect.floatClose { absolute = 1e-5; relative = 0. } (PeakShape.apex p) (31.2 + 0.1 * argMax (emgExact 1.) -3. 15.) "apex is the maximum"
            testCase "EMG with unusable tau is read as Gaussian" <| fun () ->
                for tau in [0.; -0.2; nan; infinity] do
                    let p = [| 5.; 31.2; 0.1; tau |]
                    Expect.equal (PeakShape.apex p) 31.2 (sprintf "apex, tau = %g" tau)
                    Expect.equal (PeakShape.width p) 0.1 (sprintf "width, tau = %g" tau)
            testCase "EMG with unusable sigma" <| fun () ->
                let p = [| 5.; 31.2; 0.; 0.1 |]
                Expect.equal (PeakShape.apex p) 31.2 "apex"
                Expect.equal (PeakShape.width p) 0. "width"
            testCase "fewer than 3 parameters give nan" <| fun () ->
                for p in [ [||]; [| 1. |]; [| 1.; 2. |] ] do
                    Expect.isTrue (Double.IsNaN (PeakShape.apex p)) "apex"
                    Expect.isTrue (Double.IsNaN (PeakShape.width p)) "width"
        ]
        testList "Dto helpers" [
            let row (globalMod: int) (light: float []) (heavy: float []) (differenceLight: float) (differenceHeavy: float) : Dto.QuantificationResult =
                let values = dict [ "GlobalMod", box globalMod; "Params_Light", box light; "Params_Heavy", box heavy; "Difference_SearchRT_FittedRT_Light", box differenceLight; "Difference_SearchRT_FittedRT_Heavy", box differenceHeavy ]
                let fields = Reflection.FSharpType.GetRecordFields typeof<Dto.QuantificationResult>
                let defaultOf (t: Type) = if t.IsValueType then Activator.CreateInstance t else null
                Reflection.FSharpValue.MakeRecord (typeof<Dto.QuantificationResult>, fields |> Array.map (fun f -> match values.TryGetValue f.Name with | true, v -> v | _ -> defaultOf f.PropertyType)) :?> Dto.QuantificationResult
            testCase "target label, apex, width and apex difference" <| fun () ->
                let light = [| 5.; 31.2; 0.1; 0.1 |]
                let heavy = [| 5.; 31.2; 0.2 |]
                let apex = 31.2 + 0.1 * PeakShape.emgModeOffset 1.
                let a = row 0 light heavy 0.5 0.4
                Expect.equal (Dto.QuantificationResult.getTargetScanTime a) 31.2 "parameter 1 kept"
                Expect.equal (Dto.QuantificationResult.getTargetPeakApex a) apex "apex"
                Expect.floatClose { absolute = 1e-12; relative = 0. } (Dto.QuantificationResult.getTargetPeakWidth a) (sqrt 0.02) "width"
                Expect.floatClose { absolute = 1e-12; relative = 0. } (Dto.QuantificationResult.getTargetApexDifference a) (0.5 - (apex - 31.2)) "apex difference"
                let b = row 1 light heavy 0.5 0.1
                Expect.equal (Dto.QuantificationResult.getTargetPeakApex b) 31.2 "heavy apex"
                Expect.equal (Dto.QuantificationResult.getTargetPeakWidth b) 0.2 "heavy width"
                Expect.equal (Dto.QuantificationResult.getTargetApexDifference b) 0.1 "Gaussian: the apex difference is the stored difference, bit for bit"
            testCase "no parameters give nan" <| fun () ->
                let a = row 0 [||] [||] 0.5 0.5
                Expect.isTrue (Double.IsNaN (Dto.QuantificationResult.getTargetPeakApex a)) "apex"
                Expect.isTrue (Double.IsNaN (Dto.QuantificationResult.getTargetPeakWidth a)) "width"
                Expect.isTrue (Double.IsNaN (Dto.QuantificationResult.getTargetApexDifference a)) "difference"
        ]
    ]
