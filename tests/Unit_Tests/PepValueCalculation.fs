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
    ]
