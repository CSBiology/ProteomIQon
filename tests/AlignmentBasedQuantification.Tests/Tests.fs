module AlignmentBasedQuantification.Tests

open Expecto
open ProteomIQon.AlignmentBasedQuantification
open ProteomIQon.Dto
open System

[<Tests>]
let klDivTests =
    testList "klDiv (Kullback-Leibler Divergence)" [
        testCase "Identical distributions have zero divergence" <| fun _ ->
            let p = [|0.25; 0.25; 0.25; 0.25|]
            let q = [|0.25; 0.25; 0.25; 0.25|]
            let result = klDiv p q
            Expect.floatClose Accuracy.high result 0.0 "KL divergence should be 0 for identical distributions"
        
        testCase "Calculates divergence for different distributions" <| fun _ ->
            let p = [|0.5; 0.3; 0.2|]
            let q = [|0.4; 0.4; 0.2|]
            let result = klDiv p q
            Expect.isTrue (result > 0.0) "Divergence should be positive for different distributions"
        
        testCase "Non-symmetric: KL(p||q) != KL(q||p)" <| fun _ ->
            let p = [|0.7; 0.2; 0.1|]
            let q = [|0.3; 0.5; 0.2|]
            let klPQ = klDiv p q
            let klQP = klDiv q p
            Expect.notEqual klPQ klQP "KL divergence is not symmetric"
        
        testCase "Larger divergence for more different distributions" <| fun _ ->
            let p = [|0.5; 0.3; 0.2|]
            let q1 = [|0.48; 0.32; 0.2|]
            let q2 = [|0.2; 0.5; 0.3|]
            let kl1 = klDiv p q1
            let kl2 = klDiv p q2
            Expect.isTrue (kl2 > kl1) "More different distributions have larger divergence"
        
        testCase "Isotopic pattern comparison scenario" <| fun _ ->
            let observed = [|0.60; 0.30; 0.10|]
            let theoretical = [|0.65; 0.28; 0.07|]
            let result = klDiv observed theoretical
            Expect.isTrue (result >= 0.0 && result < 0.5) "Realistic isotopic pattern divergence"
    ]

[<Tests>]
let getBaseLineCorrectionOffsetAtTests =
    testList "getBaseLineCorrectionOffsetAt" [
        testCase "Returns offset at target retention time" <| fun _ ->
            let tarRT = 25.0
            let x_Xic = [|20.0; 22.0; 24.0; 26.0; 28.0; 30.0|]
            let y_Xic = [|100.0; 150.0; 200.0; 180.0; 120.0; 80.0|]
            let y_Xic_uncorrected = [|120.0; 170.0; 230.0; 210.0; 150.0; 110.0|]
            let offset = getBaseLineCorrectionOffsetAt tarRT x_Xic y_Xic y_Xic_uncorrected
            let expected = 210.0 - 180.0
            Expect.floatClose Accuracy.high offset expected "Offset should be uncorrected - corrected"
        
        testCase "Finds closest RT when exact match not available" <| fun _ ->
            let tarRT = 23.5
            let x_Xic = [|20.0; 22.0; 26.0; 28.0|]
            let y_Xic = [|100.0; 150.0; 180.0; 120.0|]
            let y_Xic_uncorrected = [|120.0; 170.0; 210.0; 150.0|]
            let offset = getBaseLineCorrectionOffsetAt tarRT x_Xic y_Xic y_Xic_uncorrected
            let expected = 170.0 - 150.0
            Expect.floatClose Accuracy.high offset expected "Should use closest RT point"
        
        testCase "Zero offset when no baseline correction" <| fun _ ->
            let tarRT = 25.0
            let x_Xic = [|20.0; 25.0; 30.0|]
            let y_Xic = [|100.0; 200.0; 150.0|]
            let y_Xic_uncorrected = [|100.0; 200.0; 150.0|]
            let offset = getBaseLineCorrectionOffsetAt tarRT x_Xic y_Xic y_Xic_uncorrected
            Expect.floatClose Accuracy.high offset 0.0 "No offset when uncorrected equals corrected"
        
        testCase "Positive offset when baseline removed" <| fun _ ->
            let tarRT = 25.0
            let x_Xic = [|25.0|]
            let y_Xic = [|150.0|]
            let y_Xic_uncorrected = [|200.0|]
            let offset = getBaseLineCorrectionOffsetAt tarRT x_Xic y_Xic y_Xic_uncorrected
            Expect.floatClose Accuracy.high offset 50.0 "Positive offset when baseline subtracted"
    ]

[<Tests>]
let getInferredXicTests =
    testList "getInferredXic" [
        testCase "Returns XIC structure with correct dimensions" <| fun _ ->
            let getXic scanTime mz ionMob =
                let x = [|20.0; 21.0; 22.0|]
                let y = [|100.0; 200.0; 150.0|]
                let yUncorr = [|120.0; 220.0; 170.0|]
                (x, y, yUncorr)
            let result = getInferredXic getXic 25.0 500.0 None
            Expect.equal result.X_Xic.Length 3 "X_Xic should have 3 points"
            Expect.equal result.Y_Xic.Length 3 "Y_Xic should have 3 points"
            Expect.equal result.Y_Xic_uncorrected.Length 3 "Y_Xic_uncorrected should have 3 points"
        
        testCase "Preserves XIC values" <| fun _ ->
            let getXic scanTime mz ionMob =
                ([|25.0; 26.0|], [|300.0; 400.0|], [|320.0; 430.0|])
            let result = getInferredXic getXic 25.5 501.5 None
            Expect.equal result.X_Xic [|25.0; 26.0|] "X values should be preserved"
            Expect.equal result.Y_Xic [|300.0; 400.0|] "Y values should be preserved"
            Expect.equal result.Y_Xic_uncorrected [|320.0; 430.0|] "Uncorrected Y should be preserved"
    ]

[<Tests>]
let lightQualityFilterTests =
    testList "lightQualityFilter with QuantificationResult" [
        testCase "Filters based on apex/quant ratio" <| fun _ ->
            let createQuant light_apex light_quant =
                {StringSequence=""; GlobalMod=0; Charge=2; PepSequenceID=1; ModSequenceID=1; PrecursorMZ=500.0;
                 MeasuredMass=1000.0; TheoMass=1000.0; AbsDeltaMass=0.0; MeanPercolatorScore=0.0; QValue=0.01; PEPValue=0.01;
                 ProteinNames=""; QuantMz_Light=500.0; Quant_Light=light_quant; MeasuredApex_Light=light_apex;
                 Seo_Light=0.0; Params_Light=[||]; Difference_SearchRT_FittedRT_Light=0.0; KLDiv_Observed_Theoretical_Light=0.0;
                 KLDiv_CorrectedObserved_Theoretical_Light=0.0; QuantMz_Heavy=0.0; Quant_Heavy=0.0; MeasuredApex_Heavy=0.0;
                 Seo_Heavy=0.0; Params_Heavy=[||]; Difference_SearchRT_FittedRT_Heavy=0.0; KLDiv_Observed_Theoretical_Heavy=0.0;
                 KLDiv_CorrectedObserved_Theoretical_Heavy=0.0; Correlation_Light_Heavy=0.0; QuantificationSource=QuantificationSource.Alignment;
                 IsotopicPatternMz_Light=[||]; IsotopicPatternIntensity_Observed_Light=[||]; IsotopicPatternIntensity_Corrected_Light=[||];
                 RtTrace_Light=[||]; IntensityTrace_Observed_Light=[||]; IntensityTrace_Corrected_Light=[||];
                 IsotopicPatternMz_Heavy=[||]; IsotopicPatternIntensity_Observed_Heavy=[||]; IsotopicPatternIntensity_Corrected_Heavy=[||];
                 RtTrace_Heavy=[||]; IntensityTrace_Observed_Heavy=[||]; IntensityTrace_Corrected_Heavy=[||];
                 AlignmentScore=0.0; AlignmentQValue=0.0; IonMobility=0.0}
            let quants = [|createQuant 100.0 100.0; createQuant 200.0 200.0; createQuant 300.0 300.0|]
            let result = lightQualityFilter -1.0 1.0 quants
            Expect.equal result.Length 3 "All should pass with ratio ~1"
        
        testCase "Accepts NaN values" <| fun _ ->
            let createQuant light_apex light_quant =
                {StringSequence=""; GlobalMod=0; Charge=2; PepSequenceID=1; ModSequenceID=1; PrecursorMZ=500.0;
                 MeasuredMass=1000.0; TheoMass=1000.0; AbsDeltaMass=0.0; MeanPercolatorScore=0.0; QValue=0.01; PEPValue=0.01;
                 ProteinNames=""; QuantMz_Light=500.0; Quant_Light=light_quant; MeasuredApex_Light=light_apex;
                 Seo_Light=0.0; Params_Light=[||]; Difference_SearchRT_FittedRT_Light=0.0; KLDiv_Observed_Theoretical_Light=0.0;
                 KLDiv_CorrectedObserved_Theoretical_Light=0.0; QuantMz_Heavy=0.0; Quant_Heavy=0.0; MeasuredApex_Heavy=0.0;
                 Seo_Heavy=0.0; Params_Heavy=[||]; Difference_SearchRT_FittedRT_Heavy=0.0; KLDiv_Observed_Theoretical_Heavy=0.0;
                 KLDiv_CorrectedObserved_Theoretical_Heavy=0.0; Correlation_Light_Heavy=0.0; QuantificationSource=QuantificationSource.Alignment;
                 IsotopicPatternMz_Light=[||]; IsotopicPatternIntensity_Observed_Light=[||]; IsotopicPatternIntensity_Corrected_Light=[||];
                 RtTrace_Light=[||]; IntensityTrace_Observed_Light=[||]; IntensityTrace_Corrected_Light=[||];
                 IsotopicPatternMz_Heavy=[||]; IsotopicPatternIntensity_Observed_Heavy=[||]; IsotopicPatternIntensity_Corrected_Heavy=[||];
                 RtTrace_Heavy=[||]; IntensityTrace_Observed_Heavy=[||]; IntensityTrace_Corrected_Heavy=[||];
                 AlignmentScore=0.0; AlignmentQValue=0.0; IonMobility=0.0}
            let quants = [|createQuant nan 100.0; createQuant 200.0 nan|]
            let result = lightQualityFilter -1.0 1.0 quants
            Expect.equal result.Length 2 "NaN values should be accepted"
    ]

[<Tests>]
let heavyQualityFilterTests =
    testList "heavyQualityFilter with QuantificationResult" [
        testCase "Filters based on heavy apex/quant ratio" <| fun _ ->
            let createQuant heavy_apex heavy_quant =
                {StringSequence=""; GlobalMod=0; Charge=2; PepSequenceID=1; ModSequenceID=1; PrecursorMZ=500.0;
                 MeasuredMass=1000.0; TheoMass=1000.0; AbsDeltaMass=0.0; MeanPercolatorScore=0.0; QValue=0.01; PEPValue=0.01;
                 ProteinNames=""; QuantMz_Light=0.0; Quant_Light=0.0; MeasuredApex_Light=0.0;
                 Seo_Light=0.0; Params_Light=[||]; Difference_SearchRT_FittedRT_Light=0.0; KLDiv_Observed_Theoretical_Light=0.0;
                 KLDiv_CorrectedObserved_Theoretical_Light=0.0; QuantMz_Heavy=500.0; Quant_Heavy=heavy_quant; MeasuredApex_Heavy=heavy_apex;
                 Seo_Heavy=0.0; Params_Heavy=[||]; Difference_SearchRT_FittedRT_Heavy=0.0; KLDiv_Observed_Theoretical_Heavy=0.0;
                 KLDiv_CorrectedObserved_Theoretical_Heavy=0.0; Correlation_Light_Heavy=0.0; QuantificationSource=QuantificationSource.Alignment;
                 IsotopicPatternMz_Light=[||]; IsotopicPatternIntensity_Observed_Light=[||]; IsotopicPatternIntensity_Corrected_Light=[||];
                 RtTrace_Light=[||]; IntensityTrace_Observed_Light=[||]; IntensityTrace_Corrected_Light=[||];
                 IsotopicPatternMz_Heavy=[||]; IsotopicPatternIntensity_Observed_Heavy=[||]; IsotopicPatternIntensity_Corrected_Heavy=[||];
                 RtTrace_Heavy=[||]; IntensityTrace_Observed_Heavy=[||]; IntensityTrace_Corrected_Heavy=[||];
                 AlignmentScore=0.0; AlignmentQValue=0.0; IonMobility=0.0}
            let quants = [|createQuant 100.0 100.0; createQuant 200.0 200.0|]
            let result = heavyQualityFilter -1.0 1.0 quants
            Expect.equal result.Length 2 "All should pass with ratio ~1"
    ]

[<Tests>]
let scanTimeDifferenceFilterTests =
    testList "scanTimeDifferenceFilter with QuantificationResult" [
        testCase "Fails with non-positive stdev factor" <| fun _ ->
            let quants = [||]
            Expect.throws (fun _ -> lightScanTimeDifferenceFilter 0.0 quants |> ignore) "Should fail with stdev <= 0"
    ]

[<EntryPoint>]
let main args =
    Tests.runTestsInAssemblyWithCLIArgs [] args
