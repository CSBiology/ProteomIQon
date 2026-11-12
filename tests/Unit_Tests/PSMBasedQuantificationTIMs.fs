module ProteomIQon.UnitTests

open Expecto 
open ProteomIQon
open ProteomIQon.PSMBasedQuantificationTIMs
open MzIO.Processing
open MzIO.Binary
open MzIO.Processing.MzIOLinq
open BioFSharp.Mz

let private peak intensity rt mobility = 
    Peak1D(intensity, 500.0, mobility)
    |> fun p -> RtIndexEntry.AsPeak2D(p, rt)

[<Tests>]
let gaborGridTest = 
    testList "tryCreateGaborGrid" [
        testCase "aggregates peaks in the same bin and ignores invalid peaks" <| fun _ -> 
            let input = 
                [|
                    [|
                        peak 2.0 10.0000 1.0000
                        peak 3.0 10.0000 1.0004
                        peak System.Double.NaN 10.0080 1.0000
                        peak 0.0 10.0120 1.0000
                    |]
                |]
            let grid = 
                Expect.wantSome
                    (tryCreateGaborGrid input)
                    "Grid should have been created"
            
            Expect.equal grid.RtBinCount 1 "Unexpected RT bin count"
            Expect.equal grid.MobilityBinCount 1 "Unexpected Ionmobility bin count"
            Expect.equal grid.Intensities.[0].[0] 5.0 "Peaks are not aggregated correctly"
    ]

[<Tests>]
let qpsmHeaderTest =
    testList "QpsmInput.checkHeader" [
        testCase "accepts the PSMStatistics layout and rejects the former FragPipe layout" <| fun _ ->
            let withHeader (header: string) f =
                let path = System.IO.Path.GetTempFileName()
                try
                    System.IO.File.WriteAllText(path, header + "\n")
                    f path
                finally
                    System.IO.File.Delete path
            withHeader (String.concat "\t" QpsmInput.expectedColumns) QpsmInput.checkHeader
            let fragPipeHeader =
                "PSMId\tPepSequenceID\tModSequenceID\tScanTime\tCharge\tPrecursorMZ\tTheoMass\tAbsDeltaMass\tIonMobility\tPeptideLength\tMissCleavages\tExpectscore\tHyperscore\tProbability\tStringSequence\tProteinNames\tGlobalMod"
            withHeader fragPipeHeader (fun path ->
                Expect.throws (fun () -> QpsmInput.checkHeader path) "the FragPipe layout is rejected")
    ]
[<Tests>]
let qpsmEmptyTest =
    testList "QpsmInput.read" [
        testCase "an empty file and a file with the header only have no PSMs" <| fun _ ->
            for content in [ ""; String.concat "\t" QpsmInput.expectedColumns + "\n" ] do
                let path = System.IO.Path.GetTempFileName()
                try
                    System.IO.File.WriteAllText(path, content)
                    Expect.isEmpty (QpsmInput.read path) "no PSMs"
                finally
                    System.IO.File.Delete path
    ]

[<Tests>]
let isotopeClusterIonMobilityTest =
    // peaks of one frame in file order: mz, intensity, ion mobility
    let spectrum =
        let peaks =
            [|
                // M0 spread over two mobility scans inside the window
                500.0000, 100., 1.00
                500.0000, 100., 1.02
                // a foreign ion at another mobility, exactly at the predicted m/z of M+1
                500.5015, 10000., 1.50
                // M+1 spread over two mobility scans inside the window, slightly off the predicted m/z
                500.5025, 60., 1.00
                500.5025, 30., 1.01
                // M+2 only at another mobility
                501.0030, 500., 1.50
                // outside the m/z range of the pattern, at the mobility of the ion
                501.2000, 70., 1.00
                499.3000, 70., 1.00
            |]
        SortedPeaks.sortPeakArraysByMz (peaks |> Array.map (fun (m, _, _) -> m)) (peaks |> Array.map (fun (_, i, _) -> i)) (peaks |> Array.map (fun (_, _, im) -> im))
    let extract (_: string) lockMz halfWindow imLow imHigh = Xic.extractPoint spectrum lockMz halfWindow imLow imHigh
    let compare imHalfWindow ionMobility pattern =
        PeakQuantification.comparePredictedAndMeasuredIsotopicClusterAt extract (fun _ -> "ms1") 0.01 imHalfWindow [| 10.; 11.; 12. |] [| 1.; 2.; 1. |] [| 1.; 2.; 1. |] pattern 11. 500.0 ionMobility
    let pattern = [| 500.0000, 0.5; 500.5015, 0.3; 501.0030, 0.2 |]
    testList "isotope pattern comparison inside the ion mobility window" [
        testCase "inside the mobility window M+1 is the sum of the true peaks at their m/z, not the foreign peak exactly at the predicted m/z" <| fun _ ->
            let c = compare 0.05 1.005 pattern
            Expect.equal c.PeakComparisons.Length 2 "M+2 has no signal inside the window"
            Expect.floatClose Accuracy.high c.PeakComparisons.[0].Mz 500.0 "observed m/z of M0"
            Expect.floatClose Accuracy.high c.PeakComparisons.[1].Mz 500.5025 "observed m/z of M+1 is the one of the true peaks"
            Expect.floatClose Accuracy.high c.PeakComparisons.[0].MeasuredIntensity (200. / 290.) "M0 is the sum of its two mobility scans"
            Expect.floatClose Accuracy.high c.PeakComparisons.[1].MeasuredIntensity (90. / 290.) "M+1 is the sum of its two mobility scans"
            Expect.floatClose Accuracy.high c.PeakComparisons.[0].PredictedRelFrequency (0.5 / 0.8) "predicted M0 frequency"
            Expect.floatClose Accuracy.high c.PeakComparisons.[1].PredictedRelFrequency (0.3 / 0.8) "predicted M+1 frequency"
            let kl = (200. / 290.) * log ((200. / 290.) / (0.5 / 0.8)) + (90. / 290.) * log ((90. / 290.) / (0.3 / 0.8))
            Expect.floatClose Accuracy.high c.KLDiv_UnCorrected kl "KL divergence"
        testCase "without a mobility window all peaks within the m/z window are summed: the foreign peak dominates M+1 and M+2 is found too" <| fun _ ->
            // expectation of the sum of Xic.extractPoint over the m/z window: M0 200 at 500.0, M+1 10000 + 60 + 30 around 500.5015, M+2 500 at 501.003
            let c = compare nan nan pattern
            Expect.equal c.PeakComparisons.Length 3 "all three isotope peaks have signal"
            Expect.floatClose Accuracy.high c.PeakComparisons.[1].Mz ((500.5015 * 10000. + 500.5025 * 90.) / 10090.) "observed m/z of M+1 is the intensity weighted m/z of the peaks in the window"
            Expect.floatClose Accuracy.high c.PeakComparisons.[1].MeasuredIntensity (10090. / 10790.) "M+1 intensity is the sum of the peaks in the window"
        testCase "a window that holds no peak and an unknown ion mobility give an empty pattern" <| fun _ ->
            Expect.isEmpty (compare 0.05 1.3 pattern).PeakComparisons "no peak at 1.3"
            Expect.isEmpty (compare 0.05 nan pattern).PeakComparisons "nan ion mobility"
        testCase "predicted peaks outside tarMz - 0.6 < m/z <= tarMz + 1.1 are left out" <| fun _ ->
            let c = compare 0.05 1.005 [| 499.3000, 0.1; 500.0000, 0.5; 500.5015, 0.3; 501.2000, 0.1 |]
            Expect.equal c.PeakComparisons.Length 2 "the peaks at 499.3 and 501.2 have signal but are outside the range"
            Expect.floatClose Accuracy.high c.PeakComparisons.[0].PredictedRelFrequency (0.5 / 0.8) "predicted M0 frequency"
            Expect.floatClose Accuracy.high c.PeakComparisons.[1].PredictedRelFrequency (0.3 / 0.8) "predicted M+1 frequency"
        testCase "predicted peaks within 2 x tolerance of each other are one group, further apart peaks stay separate" <| fun _ ->
            let merged = PeakQuantification.mergeCloseIsotopePeaks 0.01 [| 500.5015, 0.3; 500.0000, 0.5; 500.0050, 0.1 |]
            Expect.equal merged.Length 2 "two groups"
            Expect.floatClose Accuracy.high (fst merged.[0]) ((500.0 * 0.5 + 500.005 * 0.1) / 0.6) "weighted m/z"
            Expect.floatClose Accuracy.high (snd merged.[0]) 0.6 "summed frequency"
            Expect.floatClose Accuracy.high (snd merged.[1]) 0.3 "the far peak stays"
            Expect.equal (PeakQuantification.mergeCloseIsotopePeaks 0.01 [| 500.0000, 0.5; 500.0150, 0.5 |]).Length 1 "a gap of 1.5 x tolerance gives one group"
            Expect.equal (PeakQuantification.mergeCloseIsotopePeaks 0.01 [| 500.0000, 0.5; 500.0250, 0.5 |]).Length 2 "a gap of 2.5 x tolerance gives two groups"
            let c = compare 0.05 1.005 [| 500.0000, 0.5; 500.0050, 0.1 |]
            Expect.equal c.PeakComparisons.Length 1 "the merged peak is counted once"
            Expect.floatClose Accuracy.high c.PeakComparisons.[0].PredictedRelFrequency 1.0 "normalised frequency"
    ]
