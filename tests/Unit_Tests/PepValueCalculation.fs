module ProteomIQon.UnitTestsPepValueCalculation

open Expecto
open ProteomIQon.PepValueCalculation

/// Model scores of the best PSM per spectrum: targets form a well separated population above
/// the decoys, as after a successful training.
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
    ]
