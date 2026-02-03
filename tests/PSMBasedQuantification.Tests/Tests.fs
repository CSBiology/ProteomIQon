module PSMBasedQuantification.Tests

open Expecto
open ProteomIQon

[<Tests>]
let klDivTests =
    testList "Kullback-Leibler Divergence" [
        testCase "identical distributions have zero KL divergence" <| fun _ ->
            let p = [|0.25; 0.25; 0.25; 0.25|]
            let q = [|0.25; 0.25; 0.25; 0.25|]
            let result = PSMBasedQuantification.klDiv p q
            Expect.floatClose Accuracy.high result 0.0 "KL divergence of identical distributions should be ~0"
        
        testCase "simple KL divergence calculation" <| fun _ ->
            let p = [|0.5; 0.5|]
            let q = [|0.25; 0.75|]
            let result = PSMBasedQuantification.klDiv p q
            Expect.isTrue (result > 0.0) "KL divergence should be positive"
            Expect.isTrue (result > 0.14 && result < 0.15) "KL divergence should be approximately 0.144"
        
        testCase "KL divergence with skewed distributions" <| fun _ ->
            let p = [|0.9; 0.1|]
            let q = [|0.5; 0.5|]
            let result = PSMBasedQuantification.klDiv p q
            Expect.isTrue (result > 0.0) "KL divergence should be positive"
        
        testCase "multi-element distributions" <| fun _ ->
            let p = [|0.1; 0.2; 0.3; 0.4|]
            let q = [|0.25; 0.25; 0.25; 0.25|]
            let result = PSMBasedQuantification.klDiv p q
            Expect.isTrue (result >= 0.0) "KL divergence should be non-negative"
        
        testCase "very different distributions have large KL divergence" <| fun _ ->
            let p = [|0.99; 0.01|]
            let q = [|0.01; 0.99|]
            let result = PSMBasedQuantification.klDiv p q
            Expect.isTrue (result > 1.0) "Very different distributions should have large KL divergence"
    ]

[<Tests>]
let weightedMeanTests =
    testList "Weighted Mean" [
        testCase "equal weights give arithmetic mean" <| fun _ ->
            let weights = [1.0; 1.0; 1.0; 1.0]
            let items = [1.0; 2.0; 3.0; 4.0]
            let result = PSMBasedQuantification.weightedMean weights items
            Expect.equal result 2.5 "Weighted mean with equal weights should be 2.5"
        
        testCase "weighted mean with different weights" <| fun _ ->
            let weights = [1.0; 2.0; 3.0]
            let items = [10.0; 20.0; 30.0]
            let result = PSMBasedQuantification.weightedMean weights items
            Expect.floatClose Accuracy.high result 23.333333 "Weighted mean should be ~23.33"
        
        testCase "single element" <| fun _ ->
            let weights = [5.0]
            let items = [42.0]
            let result = PSMBasedQuantification.weightedMean weights items
            Expect.equal result 42.0 "Weighted mean of single element should be that element"
        
        testCase "zero weight on one item" <| fun _ ->
            let weights = [0.0; 1.0; 1.0]
            let items = [100.0; 2.0; 4.0]
            let result = PSMBasedQuantification.weightedMean weights items
            Expect.equal result 3.0 "Item with zero weight should not contribute"
        
        testCase "heavily weighted single value" <| fun _ ->
            let weights = [0.01; 0.99]
            let items = [1000.0; 100.0]
            let result = PSMBasedQuantification.weightedMean weights items
            Expect.isTrue (abs(result - 100.0) < 10.0) "Should heavily weight the second value"
        
        testCase "handles fractional weights" <| fun _ ->
            let weights = [0.2; 0.3; 0.5]
            let items = [10.0; 20.0; 30.0]
            let result = PSMBasedQuantification.weightedMean weights items
            Expect.isTrue (result > 20.0 && result < 25.0) "Should handle fractional weights"
    ]

[<Tests>]
let createAveragePSMTests =
    testList "Create Average PSM" [
        testCase "creates PSM with correct fields" <| fun _ ->
            let meanPrecMz = 500.5
            let meanScanTime = 1200.0
            let weightedAvgScanTime = 1205.0
            let meanScore = 0.95
            let xXic = [|100.0; 200.0; 300.0|]
            let yXic = [|1000.0; 2000.0; 1500.0|]
            let yXicUncorr = [|1010.0; 2010.0; 1510.0|]
            
            let result = PSMBasedQuantification.createAveragePSM 
                            meanPrecMz meanScanTime weightedAvgScanTime 
                            meanScore xXic yXic yXicUncorr
            
            Expect.equal result.MeanPrecMz meanPrecMz "MeanPrecMz should match"
            Expect.equal result.MeanScanTime meanScanTime "MeanScanTime should match"
            Expect.equal result.WeightedAvgScanTime weightedAvgScanTime "WeightedAvgScanTime should match"
            Expect.equal result.MeanScore meanScore "MeanScore should match"
            Expect.equal result.X_Xic xXic "X_Xic should match"
            Expect.equal result.Y_Xic yXic "Y_Xic should match"
            Expect.equal result.Y_Xic_uncorrected yXicUncorr "Y_Xic_uncorrected should match"
        
        testCase "handles empty arrays" <| fun _ ->
            let result = PSMBasedQuantification.createAveragePSM 0.0 0.0 0.0 0.0 [||] [||] [||]
            Expect.equal result.X_Xic.Length 0 "Should handle empty arrays"
        
        testCase "handles large intensity values" <| fun _ ->
            let xXic = [|100.0; 200.0|]
            let yXic = [|1e6; 2e6|]
            let yXicUncorr = [|1.1e6; 2.1e6|]
            let result = PSMBasedQuantification.createAveragePSM 1500.0 100.0 100.5 0.99 xXic yXic yXicUncorr
            Expect.isTrue (result.Y_Xic.[0] > 1e5) "Should handle large intensities"
    ]

[<Tests>]
let getBaseLineCorrectionOffsetTests =
    testList "Baseline Correction Offset" [
        testCase "finds offset at exact target RT" <| fun _ ->
            let tarRT = 1200.0
            let x_Xic = [|1100.0; 1200.0; 1300.0|]
            let y_Xic = [|1000.0; 2000.0; 1500.0|]
            let y_Xic_uncorr = [|1010.0; 2010.0; 1510.0|]
            
            let offset = PSMBasedQuantification.getBaseLineCorrectionOffsetAt 
                            tarRT x_Xic y_Xic y_Xic_uncorr
            
            Expect.equal offset 10.0 "Offset should be 10.0"
        
        testCase "finds closest RT when not exact match" <| fun _ ->
            let tarRT = 1250.0
            let x_Xic = [|1100.0; 1200.0; 1300.0|]
            let y_Xic = [|1000.0; 2000.0; 1500.0|]
            let y_Xic_uncorr = [|1010.0; 2010.0; 1510.0|]
            
            let offset = PSMBasedQuantification.getBaseLineCorrectionOffsetAt 
                            tarRT x_Xic y_Xic y_Xic_uncorr
            
            Expect.isTrue (offset = 10.0 || offset = 10.0) "Should use closest RT"
        
        testCase "handles RT at boundaries" <| fun _ ->
            let tarRT = 1.0
            let x_Xic = [|100.0; 200.0; 300.0|]
            let y_Xic = [|50.0; 100.0; 150.0|]
            let y_Xic_uncorr = [|55.0; 110.0; 165.0|]
            
            let offset = PSMBasedQuantification.getBaseLineCorrectionOffsetAt 
                            tarRT x_Xic y_Xic y_Xic_uncorr
            
            Expect.equal offset 5.0 "Should find offset at first point"
        
        testCase "handles negative offset" <| fun _ ->
            let tarRT = 100.0
            let x_Xic = [|100.0; 200.0|]
            let y_Xic = [|100.0; 200.0|]
            let y_Xic_uncorr = [|95.0; 190.0|]
            
            let offset = PSMBasedQuantification.getBaseLineCorrectionOffsetAt 
                            tarRT x_Xic y_Xic y_Xic_uncorr
            
            Expect.equal offset -5.0 "Should handle negative offsets"
        
        testCase "handles mid-range RT values" <| fun _ ->
            let tarRT = 25.0
            let x_Xic = [|10.0; 20.0; 30.0; 40.0|]
            let y_Xic = [|100.0; 200.0; 300.0; 400.0|]
            let y_Xic_uncorr = [|110.0; 220.0; 330.0; 440.0|]
            
            let offset = PSMBasedQuantification.getBaseLineCorrectionOffsetAt 
                            tarRT x_Xic y_Xic y_Xic_uncorr
            
            // Should pick closest, which is either 20.0 or 30.0
            Expect.isTrue (offset = 20.0 || offset = 30.0) "Should use closest RT value"
    ]

[<Tests>]
let getInferredXicTests =
    testList "Inferred XIC" [
        testCase "creates InferredXic structure" <| fun _ ->
            let mockGetXic targetScanTime targetMz =
                let retData = [|1.0; 2.0; 3.0|]
                let itzData = [|100.0; 200.0; 300.0|]
                let uncorrectedItzData = [|110.0; 210.0; 310.0|]
                (retData, itzData, uncorrectedItzData)
            
            let result = PSMBasedQuantification.getInferredXic mockGetXic 1200.0 500.5
            
            Expect.equal result.X_Xic [|1.0; 2.0; 3.0|] "X_Xic should match"
            Expect.equal result.Y_Xic [|100.0; 200.0; 300.0|] "Y_Xic should match"
            Expect.equal result.Y_Xic_uncorrected [|110.0; 210.0; 310.0|] "Y_Xic_uncorrected should match"
        
        testCase "handles empty XIC data" <| fun _ ->
            let mockGetXic _ _ = ([||], [||], [||])
            let result = PSMBasedQuantification.getInferredXic mockGetXic 25.0 500.0
            Expect.equal result.X_Xic.Length 0 "Should handle empty arrays"
        
        testCase "passes correct parameters to getXic" <| fun _ ->
            let mutable capturedTime = 0.0
            let mutable capturedMz = 0.0
            let mockGetXic t mz = 
                capturedTime <- t
                capturedMz <- mz
                ([| 1.0 |], [| 10.0 |], [| 11.0 |])
            let _ = PSMBasedQuantification.getInferredXic mockGetXic 42.0 789.5
            Expect.equal capturedTime 42.0 "Should pass correct scan time"
            Expect.equal capturedMz 789.5 "Should pass correct m/z"
        
        testCase "handles high intensity values" <| fun _ ->
            let mockGetXic _ _ = 
                ([| 1.0; 2.0 |], [| 1e6; 2e6 |], [| 1.1e6; 2.1e6 |])
            let result = PSMBasedQuantification.getInferredXic mockGetXic 25.0 500.0
            Expect.isTrue (result.Y_Xic.[0] > 1e5) "Should handle high intensities"
    ]

[<Tests>]
let substractBaseLineTests =
    testList "Baseline Subtraction" [
        testCase "returns input unchanged for long arrays (>500)" <| fun _ ->
            let logger = Unchecked.defaultof<NLog.Logger>
            let baseLineParams : Domain.BaseLineCorrection = {
                MaxIterations = 10
                Lambda = 5
                P = 0.05
            }
            let yData = Array.init 600 (fun i -> float i)
            let result = PSMBasedQuantification.substractBaseLine logger baseLineParams yData
            Expect.equal result yData "Should return unchanged for arrays > 500"
        
        testCase "processes short arrays with baseline correction" <| fun _ ->
            let logger = Unchecked.defaultof<NLog.Logger>
            let baseLineParams : Domain.BaseLineCorrection = {
                MaxIterations = 10
                Lambda = 100
                P = 0.05
            }
            let yData = [| 10.0; 20.0; 30.0; 20.0; 10.0 |]
            let result = PSMBasedQuantification.substractBaseLine logger baseLineParams yData
            Expect.equal result.Length yData.Length "Should maintain array length"
            Expect.isTrue (result |> Array.forall (fun x -> x >= 0.0)) "All values should be >= 0"
        
        testCase "sets negative values to zero" <| fun _ ->
            let logger = Unchecked.defaultof<NLog.Logger>
            let baseLineParams : Domain.BaseLineCorrection = {
                MaxIterations = 5
                Lambda = 10
                P = 0.1
            }
            let yData = [| 1.0; 2.0; 1.0; 2.0; 1.0 |]
            let result = PSMBasedQuantification.substractBaseLine logger baseLineParams yData
            Expect.isTrue (result |> Array.forall (fun x -> x >= 0.0)) "Should set negative values to 0"
        
        testCase "handles uniform data" <| fun _ ->
            let logger = Unchecked.defaultof<NLog.Logger>
            let baseLineParams : Domain.BaseLineCorrection = {
                MaxIterations = 10
                Lambda = 100
                P = 0.05
            }
            let yData = Array.create 50 10.0
            let result = PSMBasedQuantification.substractBaseLine logger baseLineParams yData
            Expect.equal result.Length 50 "Should maintain length"
        
        testCase "handles various lambda values" <| fun _ ->
            let logger = Unchecked.defaultof<NLog.Logger>
            let baseLineParams1 : Domain.BaseLineCorrection = {
                MaxIterations = 10
                Lambda = 1
                P = 0.05
            }
            let baseLineParams2 : Domain.BaseLineCorrection = {
                MaxIterations = 10
                Lambda = 1000
                P = 0.05
            }
            let yData = [| 5.0; 10.0; 15.0; 10.0; 5.0 |]
            let result1 = PSMBasedQuantification.substractBaseLine logger baseLineParams1 yData
            let result2 = PSMBasedQuantification.substractBaseLine logger baseLineParams2 yData
            Expect.isTrue (result1.Length = result2.Length) "Should handle different lambda values"
    ]

[<EntryPoint>]
let main args =
    Tests.runTestsInAssemblyWithCLIArgs [] args
