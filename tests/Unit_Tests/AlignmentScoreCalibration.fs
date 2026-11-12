/// Tests of the calibration of the alignment score of AlignmentBasedQuantStatistics: Platt scaling and isotonic regression.
module ProteomIQon.UnitTestsAlignmentScoreCalibration

open System
open Expecto
open ProteomIQon
open ProteomIQon.AlignmentScoring

let private tolerance = { absolute = 1e-12; relative = 0. }

let private isNonDecreasing (values: float []) = values |> Array.pairwise |> Array.forall (fun (a, b) -> a <= b)

[<Tests>]
let tests =
    testList "AlignmentScoreCalibration" [

        testCase "pool adjacent violators of a hand-checked example" <| fun () ->
            // y 1, 3, 2, 4 with equal weights: 3 and 2 pool to (2.5, 2.5)
            let knots = poolAdjacentViolators [| 1.; 2.; 3.; 4. |] [| 1.; 3.; 2.; 4. |] [| 1.; 1.; 1.; 1. |]
            Expect.equal knots [| (1., 1.); (2.5, 2.5); (4., 4.) |] "pooled knots"
            // weights 1 and 3: x = (0 + 3) / 4, y = (0.9 + 0.9) / 4
            let weighted = poolAdjacentViolators [| 0.; 1. |] [| 0.9; 0.3 |] [| 1.; 3. |]
            Expect.equal weighted.Length 1 "one block"
            Expect.floatClose tolerance (fst weighted.[0]) 0.75 "weighted x"
            Expect.floatClose tolerance (snd weighted.[0]) 0.45 "weighted y"
            // a cascade: 4, 3, 2, 1 pools into one block
            let cascade = poolAdjacentViolators [| 0.; 1.; 2.; 3. |] [| 4.; 3.; 2.; 1. |] [| 1.; 1.; 1.; 1. |]
            Expect.equal cascade [| (1.5, 2.5) |] "cascade"

        testCase "bin size is n / 20 clamped to 100 to 2000" <| fun () ->
            Expect.equal (isotonicBinSize 10) 100 "small"
            Expect.equal (isotonicBinSize 4000) 200 "n / 20"
            Expect.equal (isotonicBinSize 1000000) 2000 "large"

        testCase "isotonic fit with two bins" <| fun () ->
            // 200 rows, bin size 100: scores 0 .. 199, the first bin has no positive, the second bin 100 positives
            let scores = Array.init 200 float
            let labels = Array.init 200 (fun i -> i >= 100)
            let knots = fitIsotonic scores labels
            Expect.equal knots.X [| 49.5; 149.5 |] "mean scores"
            Expect.floatClose tolerance knots.Y.[0] (1. / 102.) "first y"
            Expect.floatClose tolerance knots.Y.[1] (101. / 102.) "second y"

        testCase "a small last bin joins the bin before, a larger one stays" <| fun () ->
            let fit n = fitIsotonic (Array.init n float) (Array.init n (fun i -> i % 2 = 0 && i >= 100))
            // 249 rows: bins 0..99 and 100..248 (49 < 50 rows left over)
            Expect.equal (fit 249).X.Length 2 "merged"
            Expect.floatClose tolerance (fit 249).X.[1] 174. "mean of 100 .. 248"
            // 250 rows: bins 0..99, 100..199 and 200..249 (50 rows are not fewer than 50)
            Expect.equal (fit 250).X.Length 3 "kept"

        testCase "bins are cut by score, ties by index" <| fun () ->
            // scores in descending order of the index: the lowest scores are the last rows
            let n = 200
            let scores = Array.init n (fun i -> float (n - 1 - i))
            let labels = Array.init n (fun i -> i < 100)
            let knots = fitIsotonic scores labels
            Expect.equal knots.X [| 49.5; 149.5 |] "mean scores"
            // rows 0 .. 99 have the highest scores and are positive
            Expect.floatClose tolerance knots.Y.[1] (101. / 102.) "high scores are positive"

        testCase "isotonic knots are monotone and probabilities within 0 and 1 on random data" <| fun () ->
            let rnd = Random(1)
            for n in [| 50; 300; 5000; 60000 |] do
                let scores = Array.init n (fun _ -> rnd.NextDouble() * 10. - 5.)
                // the labels follow a logistic probability with noise, so that the bin shares are not monotone by chance
                let labels = scores |> Array.map (fun s -> rnd.NextDouble() < 1. / (1. + exp (-s)) * 0.8 + 0.1)
                let knots = fitIsotonic scores labels
                Expect.isTrue (isNonDecreasing knots.X) (sprintf "x ascending, n %i" n)
                Expect.isTrue (isNonDecreasing knots.Y) (sprintf "y non-decreasing, n %i" n)
                Expect.isTrue (knots.Y |> Array.forall (fun y -> y > 0. && y < 1.)) (sprintf "y in (0, 1), n %i" n)
                let predicted = Array.init 2001 (fun i -> predictIsotonic knots (-6. + float i * 0.006))
                Expect.isTrue (isNonDecreasing predicted) (sprintf "prediction non-decreasing, n %i" n)
                Expect.isTrue (predicted |> Array.forall (fun y -> y >= 0. && y <= knots.Y.[knots.Y.Length - 1])) (sprintf "prediction between 0 and the last y, n %i" n)

        testCase "prediction at, between and outside the knots" <| fun () ->
            // a flat tail (a = 0): the first y below the first knot
            let knots = { X = [| 1.; 2.; 4. |]; Y = [| 0.1; 0.3; 0.9 |]; Tail = (0., 0.) }
            Expect.floatClose tolerance (predictIsotonic knots 0.) 0.1 "below the first knot"
            Expect.floatClose tolerance (predictIsotonic knots 1.) 0.1 "at the first knot"
            Expect.floatClose tolerance (predictIsotonic knots 1.5) 0.2 "between the first two knots"
            Expect.floatClose tolerance (predictIsotonic knots 2.) 0.3 "at the middle knot"
            Expect.floatClose tolerance (predictIsotonic knots 3.) 0.6 "between the last two knots"
            Expect.floatClose tolerance (predictIsotonic knots 4.) 0.9 "at the last knot"
            Expect.floatClose tolerance (predictIsotonic knots 7.) 0.9 "above the last knot"
            let single = { X = [| 1. |]; Y = [| 0.4 |]; Tail = (0., 0.) }
            Expect.floatClose tolerance (predictIsotonic single -3.) 0.4 "single knot below"
            Expect.floatClose tolerance (predictIsotonic single 3.) 0.4 "single knot above"

        testCase "a score that is not finite has no probability" <| fun () ->
            let knots = { X = [| 1.; 2.; 4. |]; Y = [| 0.1; 0.3; 0.9 |]; Tail = (1., 0.) }
            for s in [ nan; infinity; -infinity ] do
                Expect.isTrue (Double.IsNaN (predictIsotonic knots s)) (sprintf "nan for %g" s)

        testCase "below the first knot the logistic tail meets the first knot and falls towards 0" <| fun () ->
            let knots = { X = [| 1.; 2. |]; Y = [| 0.4; 0.8 |]; Tail = (1.5, -0.5) }
            let tail = platt (1.5, -0.5)
            Expect.floatClose tolerance (predictIsotonic knots 1.) 0.4 "at the first knot"
            Expect.floatClose { absolute = 1e-6; relative = 0. } (predictIsotonic knots 0.999999) 0.4 "continuous at the first knot"
            for s in [| 0.5; 0.; -1.; -3.; -8. |] do
                Expect.floatClose tolerance (predictIsotonic knots s) (0.4 * tail s / tail 1.) (sprintf "scaled logistic at %g" s)
            let below = [| 1.; 0.5; 0.; -1.; -3.; -8.; -20. |] |> Array.map (predictIsotonic knots)
            Expect.isTrue (below |> Array.pairwise |> Array.forall (fun (higher, lower) -> lower < higher)) "strictly decreasing for lower scores"
            Expect.isTrue (below.[1 ..] |> Array.forall (fun p -> p < 0.4)) "strictly below the first y"
            Expect.isLessThan below.[6] 1e-6 "towards 0"

        testCase "a tail without positive slope or a non-finite tail value keeps the first y" <| fun () ->
            let flat = { X = [| 1.; 2. |]; Y = [| 0.4; 0.8 |]; Tail = (0., 0.3) }
            let negative = { flat with Tail = (-1., 0.) }
            Expect.floatClose tolerance (predictIsotonic flat -5.) 0.4 "a = 0"
            Expect.floatClose tolerance (predictIsotonic negative -5.) 0.4 "a < 0"
            // platt at the first knot (-1) underflows to 0 for a huge slope, the ratio is not finite
            let huge = { flat with X = [| -1.; 2. |]; Tail = (1e6, 0.) }
            Expect.floatClose tolerance (predictIsotonic huge -5.) 0.4 "non-finite ratio"

        testCase "the tail of the fit is fitPlatt of the rows" <| fun () ->
            let rnd = Random(3)
            let scores = Array.init 400 (fun _ -> rnd.NextDouble() * 6. - 3.)
            let labels = scores |> Array.map (fun s -> rnd.NextDouble() < 1. / (1. + exp (-s)))
            let knots = fitIsotonic scores labels
            Expect.equal knots.Tail (fitPlatt scores labels) "tail"
            Expect.isGreaterThan (fst knots.Tail) 0. "positive slope"

        testCase "bins with the same mean score are merged into one knot" <| fun () ->
            // 200 rows with one score, bin size 100, 20 and 60 positives in index order: the bins have y 21/102 and 61/102
            let scores = Array.create 200 3.
            let labels = Array.init 200 (fun i -> if i < 100 then i < 20 else i - 100 < 60)
            let knots = fitIsotonic scores labels
            Expect.equal knots.X [| 3. |] "one knot"
            Expect.floatClose tolerance knots.Y.[0] ((21. / 102. + 61. / 102.) / 2.) "mean of the two bin values"

        testCase "calibrate applies the calibrator kind" <| fun () ->
            let knots = { X = [| 1.; 3. |]; Y = [| 0.2; 0.6 |]; Tail = (0., 0.) }
            Expect.floatClose tolerance (calibrate (IsotonicRegression knots) 2.) 0.4 "isotonic"
            Expect.floatClose tolerance (calibrate (PlattScaling (2., -1.)) 0.5) 0.5 "Platt at the midpoint"

        testCase "the Platt path is the logistic function of the fitted parameters" <| fun () ->
            // p(-1) = 1/3 and p(1) = 2/3 in the data: the maximum likelihood fit is a = ln 2, b = 0
            let scores = [| -1.; -1.; -1.; 1.; 1.; 1. |]
            let labels = [| false; false; true; true; true; false |]
            let a, b = fitPlatt scores labels
            Expect.floatClose { absolute = 1e-6; relative = 0. } a (log 2.) "a"
            Expect.floatClose { absolute = 1e-6; relative = 0. } b 0. "b"
            for s in [| -3.; -1.; 0.; 0.25; 1.; 4. |] do
                Expect.equal (calibrate (PlattScaling (a, b)) s) (platt (a, b) s) "calibrate equals platt"
                Expect.equal (calibrate (PlattScaling (a, b)) s) (1. / (1. + exp (-(a * s + b)))) "logistic formula"
    ]
