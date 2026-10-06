module ProteomIQon.UnitTestsPepValueCalculation

open Expecto
open ProteomIQon.PepValueCalculation

/// Model scores of the best PSM per spectrum: targets form a well separated population above
/// the decoys, as after a successful training. Scores are sums of three uniforms, close to a
/// normal distribution and bounded at three standard deviations.
let private population (seed: int) =
    let rnd = System.Random(seed)
    let normal mean sd = mean + sd * (rnd.NextDouble() + rnd.NextDouble() + rnd.NextDouble() - 1.5) * 2.
    let targets = Array.init 4000 (fun _ -> normal 15. 6., false)
    let decoys = Array.init 2000 (fun _ -> normal -5. 5., true)
    Array.append targets decoys

let private pepOf (data: (float * bool)[]) =
    let logger = NLog.LogManager.GetLogger "PepValueCalculationTests"
    initCalculatePEPValueIRLS logger 1. snd fst fst data

[<Tests>]
let pepValueTests =
    testList "initCalculatePEPValueIRLS" [
        testCase "gives a finite PEP for scores shared by several targets" <| fun _ ->
            // Tree models give many PSMs the same score
            let data = population 11 |> Array.map (fun (score, isDecoy) -> round (score * 2.) / 2., isDecoy)
            let pep = pepOf data
            for score, isDecoy in data do
                if not isDecoy then
                    let p = pep score
                    Expect.isTrue (System.Double.IsFinite p && p >= 0. && p <= 1.) (sprintf "PEP %g at the tied target score %g" p score)

        testCase "gives a PEP when every target has the same score" <| fun _ ->
            let data = population 11 |> Array.map (fun (score, isDecoy) -> (if isDecoy then score else 20.), isDecoy)
            let p = (pepOf data) 20.
            Expect.isTrue (System.Double.IsFinite p && p >= 0. && p <= 1.) (sprintf "PEP %g at the single target score" p)

        testCase "keeps low PEPs below a few high scoring decoys" <| fun _ ->
            // A handful of outliers, decoys among them, receive the highest scores
            let pocket = [| 48., true; 48.5, false; 49., true; 49.5, true; 50., false; 51., true; 52., false |]
            let data = Array.append (population 23) pocket
            let pep = pepOf data
            Expect.isTrue (pep 30. < 0.05) (sprintf "A target far above the decoys has a low PEP, got %g" (pep 30.))
            Expect.isTrue (pep 20. < 0.05) (sprintf "A target in the bulk above the decoys has a low PEP, got %g" (pep 20.))
            let withoutPocket = pepOf (population 23)
            Expect.isTrue (pep 50. > withoutPocket 50.) (sprintf "The decoys among the highest scores raise the PEP there, got %g with them and %g without" (pep 50.) (withoutPocket 50.))
            Expect.isTrue (pep 50. > 1e-4) (sprintf "Four decoys among seven PSMs at the top keep the PEP there above 1e-4, got %g" (pep 50.))
            let targetScores = data |> Array.filter (snd >> not) |> Array.map fst |> Array.sort
            for low, high in Array.pairwise targetScores do
                Expect.isTrue (pep high <= pep low) (sprintf "The PEP does not rise from score %g to %g" low high)

        testCase "gives PEPs when the scores fill fewer than four histogram bins" <| fun _ ->
            // A saturated tree model puts almost every PSM at one of two scores
            let data =
                Array.concat [ Array.replicate 600 (51.918, false); Array.replicate 5 (51.918, true)
                               Array.replicate 300 (-51.918, true); Array.replicate 40 (-51.918, false) ]
            let pep = pepOf data
            Expect.isTrue (pep 51.918 < 0.05) (sprintf "The high score with few decoys has a low PEP, got %g" (pep 51.918))
            Expect.isTrue (pep -51.918 > 0.5) (sprintf "The low score dominated by decoys has a high PEP, got %g" (pep -51.918))

        testCase "gives PEPs when all scores fill one histogram bin" <| fun _ ->
            let data = Array.append (Array.replicate 600 (51.918, false)) (Array.replicate 5 (51.918, true))
            let pep = pepOf data
            // the decoy/target ratio of the bin, with the pseudo counts of the spline start
            let expected = 5.05 / 600.05
            for score in [ -10.; 51.918; 80. ] do
                Expect.floatClose Accuracy.high (pep score) expected (sprintf "The single bin gives its decoy/target ratio at score %g" score)

        testCase "gives PEPs through the middle bin when the scores fill three histogram bins" <| fun _ ->
            let data =
                Array.concat [ Array.replicate 600 (50., false); Array.replicate 5 (50., true)
                               Array.replicate 200 (0., false); Array.replicate 20 (0., true)
                               Array.replicate 40 (-50., false); Array.replicate 300 (-50., true) ]
            let pep = pepOf data
            Expect.floatClose Accuracy.high (pep 0.) (20.05 / 200.05) "The middle bin gives its decoy/target ratio"
            Expect.isTrue (pep 50. < pep 25. && pep 25. < pep 0.) (sprintf "Between two bins the PEP lies between theirs, got %g < %g < %g" (pep 50.) (pep 25.) (pep 0.))
            Expect.equal (pep -50.) 1. "The bin dominated by decoys has PEP 1"

        testCase "gives PEP 1 without any PSM" <| fun _ ->
            let pep = pepOf [||]
            Expect.equal (pep 10.) 1. "Without PSMs every PEP is 1"
    ]
