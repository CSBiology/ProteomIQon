module ProteomIQon.UnitTestsPeptideSpectrumMatchingTIMs

open Expecto
open ProteomIQon.SpectrumSearch

/// Candidate hyperscores of one spectrum: a random population with a decreasing tail and a
/// few high scoring candidates separated from it by a gap, as a true hit and its homologs are.
let private candidateScores (random: bool) =
    let rnd = System.Random(7)
    let population = Array.init 600 (fun _ -> 4. + 2.5 * abs (rnd.NextDouble() + rnd.NextDouble() + rnd.NextDouble() - 1.5))
    if random then population else Array.append population [| 22.; 23.5; 24. |]

[<Tests>]
let expectationModelTests =
    testList "expectationModel" [
        testCase "fits the survival function and decreases with the hyperscore" <| fun _ ->
            let expect = expectationModel (candidateScores false)
            let atMode = expect 6.
            let atTop = expect 24.
            Expect.isTrue (atMode > 10.) "A score in the bulk of the candidates is expected many times"
            Expect.isTrue (atTop < 0.01) "The top score is expected less than once in a hundred"
            Expect.isTrue (expect 12. > expect 16.) "The expectation value decreases with the score"
            Expect.isTrue (expect 200. >= 1e-15) "The expectation value is floored"
            let scale = 4. / log 10.
            let defaultLine = 10. ** (3.5 - 0.18 * 24. * scale)
            Expect.isFalse (abs (atTop - defaultLine) < 1e-12 * defaultLine) "603 candidates are fitted, the default line is not used"

        testCase "ignores the candidates above the gap when it fits the line" <| fun _ ->
            let withHits = expectationModel (candidateScores false)
            let withoutHits = expectationModel (candidateScores true)
            Expect.floatClose Accuracy.medium (withHits 15.) (withoutHits 15.) "Candidates above the gap do not change the fitted line"

        testCase "uses the default line of X!Tandem below 200 candidates" <| fun _ ->
            let expect = expectationModel (Array.init 50 (fun i -> 5. + float i * 0.1))
            let scale = 4. / log 10.
            Expect.floatClose Accuracy.high (expect 20.) (10. ** (3.5 - 0.18 * 20. * scale)) "Default line"
    ]
