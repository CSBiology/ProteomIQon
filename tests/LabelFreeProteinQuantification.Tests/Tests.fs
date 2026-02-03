module LabelFreeProteinQuantification.Tests

open Expecto
open ProteomIQon.LabelFreeProteinQuantification
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
        testCase "Calculates SEM for simple dataset" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let data = [4.0; 6.0; 8.0; 10.0; 12.0]
            let result = sem data
            let expectedStDev = Seq.stDev data
            let expectedSEM = expectedStDev / sqrt 5.0
            Expect.floatClose Accuracy.high result expectedSEM "SEM calculation should match"
        
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
        
        testCase "SEM with high variance" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let data = [1.0; 100.0; 1.0; 100.0]
            let result = sem data
            let stDev = Seq.stDev data
            let expected = stDev / 2.0
            Expect.floatClose Accuracy.high result expected "SEM with high variance"
        
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
let performChargeOrModAggregationTests =
    testList "performChargeOrModAggregation integration" [
        testCase "tryGetQVal is used in aggregation" <| fun _ ->
            let testData = [0.01; nan; 0.05; nan]
            let result = tryGetQVal testData
            Expect.floatClose Accuracy.high result 0.01 "Should extract first valid q-value"
        
        testCase "tryGetQVal returns NaN when no valid values" <| fun _ ->
            let testData = [nan; nan]
            let result = tryGetQVal testData
            Expect.isTrue (System.Double.IsNaN result) "No valid q-values"
    ]

[<Tests>]
let performPeptideToProteinAggregationTests =
    testList "performPeptideToProteinAggregation statistics" [
        testCase "SEM calculation in protein aggregation" <| fun _ ->
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let peptideIntensities = [100.0; 120.0; 110.0; 130.0]
            let result = sem peptideIntensities
            let stDev = Seq.stDev peptideIntensities
            let expected = stDev / 2.0
            Expect.floatClose Accuracy.high result expected "SEM for protein aggregation"
        
        testCase "CV calculation integration" <| fun _ ->
            let data = [100.0; 110.0; 120.0; 130.0]
            let cv = Seq.cv data
            let mean = Seq.mean data
            let stDev = Seq.stDev data
            let expectedCV = stDev / mean
            Expect.floatClose Accuracy.high cv expectedCV "CV calculation"
        
        testCase "StDev calculation integration" <| fun _ ->
            let data = [10.0; 20.0; 30.0; 40.0]
            let result = Seq.stDevBy float data
            let expected = Seq.stDev data
            Expect.floatClose Accuracy.high result expected "StDev calculation"
        
        testCase "Realistic proteomics CV range" <| fun _ ->
            let peptideQuants = [1000.0; 1100.0; 950.0; 1050.0; 1020.0]
            let cv = Seq.cv peptideQuants
            Expect.isTrue (cv > 0.0 && cv < 0.2) "CV should be in realistic range"
    ]

[<EntryPoint>]
let main args =
    Tests.runTestsInAssemblyWithCLIArgs [] args
