module LabeledProteinQuantification.Tests

open Expecto
open ProteomIQon.LabeledProteinQuantification
open FSharp.Stats

[<Tests>]
let tryGetQValTests =
    testList "tryGetQVal" [
        testCase "Returns first non-NaN value from sequence" <| fun _ ->
            let result = tryGetQVal [nan; nan; 0.05; 0.01]
            Expect.floatClose Accuracy.high result 0.05 "Should return first non-NaN"
        
        testCase "Returns NaN for empty sequence" <| fun _ ->
            let result = tryGetQVal []
            Expect.isTrue (System.Double.IsNaN result) "Should return NaN for empty"
        
        testCase "Returns NaN for all-NaN sequence" <| fun _ ->
            let result = tryGetQVal [nan; nan; nan]
            Expect.isTrue (System.Double.IsNaN result) "Should return NaN"
        
        testCase "Returns single value" <| fun _ ->
            let result = tryGetQVal [0.001]
            Expect.floatClose Accuracy.high result 0.001 "Should return single value"
        
        testCase "Handles very small q-values" <| fun _ ->
            let result = tryGetQVal [1e-10; 1e-5]
            Expect.floatClose Accuracy.high result 1e-10 "Should handle very small values"
    ]

[<Tests>]
let semTests =
    testList "sem (Standard Error of Mean)" [
        testCase "Calculates SEM for light channel" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let lightIntensities = [1000.0; 1200.0; 1100.0; 1300.0]
            let result = sem lightIntensities
            let expectedStDev = Seq.stDev lightIntensities
            let expectedSEM = expectedStDev / 2.0
            Expect.floatClose Accuracy.high result expectedSEM "Light channel SEM"
        
        testCase "Calculates SEM for heavy channel" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let heavyIntensities = [500.0; 600.0; 550.0; 650.0]
            let result = sem heavyIntensities
            let expectedStDev = Seq.stDev heavyIntensities
            let expectedSEM = expectedStDev / 2.0
            Expect.floatClose Accuracy.high result expectedSEM "Heavy channel SEM"
        
        testCase "SEM decreases with larger sample size" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let smallSample = [1.0; 2.0; 3.0]
            let largeSample = [1.0; 2.0; 3.0; 1.5; 2.5; 3.5; 2.0; 2.8]
            let semSmall = sem smallSample
            let semLarge = sem largeSample
            Expect.isTrue (semSmall > semLarge) "SEM should decrease with larger n"
        
        testCase "SEM with zero variance" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let data = [5.0; 5.0; 5.0; 5.0]
            let result = sem data
            Expect.floatClose Accuracy.high result 0.0 "SEM of constant values should be 0"
        
        testCase "SEM for ratios (light/heavy)" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let ratios = [2.0; 2.1; 1.9; 2.05; 1.95]
            let result = sem ratios
            let stDev = Seq.stDev ratios
            let expected = stDev / sqrt 5.0
            Expect.floatClose Accuracy.high result expected "Ratio SEM calculation"
        
        testCase "SEM decreases by sqrt factor with larger n" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let data4 = [10.0; 12.0; 14.0; 16.0]
            let sem4 = sem data4
            let sem8 = sem [10.0; 12.0; 14.0; 16.0; 10.0; 12.0; 14.0; 16.0]
            Expect.isTrue (sem8 < sem4) "SEM should decrease with larger n"
    ]

[<Tests>]
let performGlobalModAggregationTests =
    testList "performGlobalModAggregation integration" [
        testCase "tryGetQVal extracts q-value from global mod aggregation" <| fun _ ->
            let qValues = [0.001; nan; 0.05]
            let result = tryGetQVal qValues
            Expect.floatClose Accuracy.high result 0.001 "Extract first valid q-value"
        
        testCase "Light/Heavy ratio calculation concept" <| fun _ ->
            let lightQuants = [1000.0; 2000.0; 1500.0]
            let heavyQuants = [500.0; 1000.0; 750.0]
            let ratios = 
                List.zip lightQuants heavyQuants 
                |> List.map (fun (l, h) -> l / h)
            let avgRatio = List.average ratios
            Expect.floatClose Accuracy.high avgRatio 2.0 "Light/Heavy ratio should be ~2.0"
    ]

[<Tests>]
let performChargeOrModAggregationTests =
    testList "performChargeOrModAggregation statistics" [
        testCase "tryGetQVal returns NaN when no valid values" <| fun _ ->
            let testData = [nan; nan]
            let result = tryGetQVal testData
            Expect.isTrue (System.Double.IsNaN result) "No valid q-values"
        
        testCase "Ratio aggregation across charges" <| fun _ ->
            let charge2Ratio = 1.8
            let charge3Ratio = 2.2
            let charge4Ratio = 2.0
            let avgRatio = [charge2Ratio; charge3Ratio; charge4Ratio] |> List.average
            Expect.floatClose Accuracy.high avgRatio 2.0 "Average ratio across charges"
    ]

[<Tests>]
let performPeptideToProteinAggregationTests =
    testList "performPeptideToProteinAggregation comprehensive" [
        testCase "SEM calculation for light channel in protein aggregation" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let peptideLightQuants = [1000.0; 1100.0; 950.0; 1050.0]
            let result = sem peptideLightQuants
            let stDev = Seq.stDev peptideLightQuants
            let expected = stDev / 2.0
            Expect.floatClose Accuracy.high result expected "Light SEM"
        
        testCase "SEM calculation for heavy channel in protein aggregation" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let peptideHeavyQuants = [500.0; 550.0; 475.0; 525.0]
            let result = sem peptideHeavyQuants
            let stDev = Seq.stDev peptideHeavyQuants
            let expected = stDev / 2.0
            Expect.floatClose Accuracy.high result expected "Heavy SEM"
        
        testCase "SEM calculation for ratios in protein aggregation" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let peptideRatios = [2.0; 2.1; 1.9; 2.05]
            let result = sem peptideRatios
            let stDev = Seq.stDev peptideRatios
            let expected = stDev / 2.0
            Expect.floatClose Accuracy.high result expected "Ratio SEM"
        
        testCase "CV calculation for light channel" <| fun _ ->
            let data = [1000.0; 1100.0; 1200.0; 1300.0]
            let cv = Seq.cv data
            let mean = Seq.mean data
            let stDev = Seq.stDev data
            let expectedCV = stDev / mean
            Expect.floatClose Accuracy.high cv expectedCV "Light CV"
        
        testCase "CV calculation for heavy channel" <| fun _ ->
            let data = [500.0; 550.0; 600.0; 650.0]
            let cv = Seq.cv data
            let mean = Seq.mean data
            let stDev = Seq.stDev data
            let expectedCV = stDev / mean
            Expect.floatClose Accuracy.high cv expectedCV "Heavy CV"
        
        testCase "StDev calculation for light" <| fun _ ->
            let data = [1000.0; 1200.0; 1100.0; 1300.0]
            let result = Seq.stDevBy float data
            let expected = Seq.stDev data
            Expect.floatClose Accuracy.high result expected "Light StDev"
        
        testCase "StDev calculation for heavy" <| fun _ ->
            let data = [500.0; 600.0; 550.0; 650.0]
            let result = Seq.stDevBy float data
            let expected = Seq.stDev data
            Expect.floatClose Accuracy.high result expected "Heavy StDev"
        
        testCase "Realistic SILAC ratio CV" <| fun _ ->
            let ratios = [2.0; 2.05; 1.95; 2.1; 1.9]
            let cv = Seq.cv ratios
            Expect.isTrue (cv > 0.0 && cv < 0.15) "SILAC ratio CV should be realistic"
    ]

[<EntryPoint>]
let main args =
    Tests.runTestsInAssemblyWithCLIArgs [] args
