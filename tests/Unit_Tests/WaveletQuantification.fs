module ProteomIQon.WaveletQuantificationTests

open System
open Expecto
open ProteomIQon
open BioFSharp.Mz.Quantification
open FSharp.Stats

/// The drop-out rule of the XIC processing before it was shared: a zero point after a point more than 100 times higher is removed.
let private oldDropOutFilter (tmp: float []) =
    tmp
    |> Array.mapi (fun i intensity ->
        if i = 0 || i = tmp.Length - 1 || intensity > 0. then Some intensity
        else
            let intensity' = tmp.[i - 1]
            if intensity' = 0. then Some intensity
            elif intensity' > (100. * (intensity + 1.)) then None
            else Some intensity)
    |> Array.choose id

let private randomXic (rnd: Random) n =
    Array.init n (fun _ ->
        match rnd.Next 4 with
        | 0 -> 0.
        | 1 -> float (rnd.Next 3)
        | 2 -> rnd.NextDouble() * 1000.
        | _ -> rnd.NextDouble() * 1e6)

let private baselineParams : Domain.BaseLineCorrection = { MaxIterations = 10; Lambda = 6; P = 0.05 }

/// The peak extraction of PSMBasedQuantificationTIMs before it was shared: the peaks with lowMz < m/z < highMz and the ion mobility
/// window are summed per exact m/z, the sum closest to the lock m/z is taken, the first in file order on a tie.
let private oldExtract (peaks: SortedPeaks.MzSortedPeaks) lockMz halfWindow imLow imHigh =
    let low, high = lockMz - halfWindow, lockMz + halfWindow
    let n = peaks.Mz.Length
    let mutable k = SortedSearch.upperBound peaks.Mz low
    let mutable found = false
    let mutable bestMz = nan
    let mutable bestIntensity = 0.
    let mutable bestDistance = nan
    let mutable bestFirst = Int32.MaxValue
    while k < n && peaks.Mz.[k] < high do
        let mz = peaks.Mz.[k]
        let mutable sum = 0.
        let mutable first = -1
        while k < n && peaks.Mz.[k] = mz do
            let im = peaks.IonMobility.[k]
            if im < imHigh && im > imLow then
                sum <- sum + peaks.Intensity.[k]
                if first < 0 then first <- peaks.Index.[k]
            k <- k + 1
        if first >= 0 then
            let distance = Math.Abs(mz - lockMz)
            let c = if found then distance.CompareTo bestDistance else -1
            if c < 0 || (c = 0 && first < bestFirst) then
                found <- true
                bestMz <- mz
                bestIntensity <- sum
                bestDistance <- distance
                bestFirst <- first
    if found then bestIntensity, bestMz else 0., nan

let private waveletParameters : Domain.WaveletParameters =
    {
        Borderpadding = None
        BorderPadMethod = Signal.Padding.BorderPaddingMethod.Zero
        InternalPaddingMethod = Signal.Padding.InternalPaddingMethod.LinearInterpolation
        HugeGapPaddingMethod = Signal.Padding.HugeGapPaddingMethod.Zero
        HugeGapPaddingDistance = 100.
        MinPeakDistance = None
        MinPeakLength = Some 0.1
        MaxPeakLength = 1.5
        NoiseQuantile = 0.01
        MinSNR = 0.01
    }

/// A fitted Gaussian shaped peak of a synthetic XIC.
let private fittedPeak center =
    let xs = [| 0. .. 0.05 .. 10. |]
    let ys = xs |> Array.map (fun x -> 1000. * exp (-((x - center) ** 2.) / (2. * 0.2 * 0.2)))
    let peaks = FSharpStats'.Wavelet.identify waveletParameters xs ys
    HULQ.quantifyPeak (peaks |> Array.maxBy (fun p -> p.Apex.YVal)), xs

[<Tests>]
let all =
    testList "Wavelet quantification helpers" [
        testCase "keepMask removes the points the drop-out filter of the XIC processing removed" <| fun _ ->
            let rnd = Random 7
            for n in [ 0; 1; 2; 3; 10; 50; 300 ] do
                for _ in 1 .. 50 do
                    let y = randomXic rnd n
                    let kept = Xic.keepMask y |> Array.zip y |> Array.filter snd |> Array.map fst
                    Expect.equal kept (oldDropOutFilter y) "same points kept"

        testCase "subtractBaseline corrects XICs of 3 to 500 points and leaves the others" <| fun _ ->
            let rnd = Random 11
            for n in [ 0; 1; 2 ] do
                let y = randomXic rnd n
                Expect.equal (Xic.subtractBaseline (Some baselineParams) y) y (sprintf "%i points: unchanged" n)
            for n in [ 3; 4; 60; 500 ] do
                let y = randomXic rnd n
                let expected =
                    BaselineAls.baselineAls baselineParams.MaxIterations baselineParams.Lambda baselineParams.P y
                    |> Array.map2 (fun v b -> if v - b < 0. then 0. else v - b) y
                let corrected = Xic.subtractBaseline (Some baselineParams) y
                Expect.equal corrected expected (sprintf "%i points: baseline subtracted" n)
                Expect.isTrue (corrected |> Array.forall (fun v -> v >= 0.)) "not negative"
            let long = randomXic rnd 501
            Expect.equal (Xic.subtractBaseline (Some baselineParams) long) long "501 points: unchanged"
            let y = randomXic rnd 20
            let copy = Xic.subtractBaseline None y
            Expect.equal copy y "without baseline correction: unchanged"
            Expect.isFalse (obj.ReferenceEquals(copy, y)) "a copy"

        testCase "baselineAls equals the FSharp.Stats baseline bit for bit" <| fun _ ->
            let rnd = Random 23
            let bits (a: float []) = a |> Array.map BitConverter.DoubleToInt64Bits
            let peakLike n =
                let center = rnd.NextDouble() * float n
                let width = 1. + rnd.NextDouble() * float n / 5.
                Array.init n (fun i ->
                    let x = float i
                    1e5 * exp (-((x - center) ** 2.) / (2. * width * width)) + 50. * rnd.NextDouble())
            for n in [ 3; 4; 5; 7; 10; 20; 60; 150; 500 ] do
                for _ in 1 .. 20 do
                    for y in [ randomXic rnd n; peakLike n ] do
                        let expected =
                            FSharp.Stats.Signal.Baseline.baselineAls' 10 6 0.05 y |> Array.ofSeq
                        let actual = BaselineAls.baselineAls 10 6 0.05 y
                        Expect.equal (bits actual) (bits expected) (sprintf "%i points: same baseline" n)

        testCase "processXic returns the points without drop-outs, corrected and raw" <| fun _ ->
            let rt = [| 1.; 2.; 3.; 4.; 5.; 6. |]
            let y = [| 0.; 500.; 0.; 0.; 20.; 0. |]
            let rt', corrected, raw = Xic.processXic None rt y
            Expect.equal rt' [| 1.; 2.; 4.; 5.; 6. |] "the zero after the 500 is removed"
            Expect.equal raw [| 0.; 500.; 0.; 20.; 0. |] "raw intensities"
            Expect.equal corrected raw "no baseline correction"

        testCase "isValidFit needs 3 finite parameters and a finite area" <| fun _ ->
            let fit ps area = HULQ.createQuantifiedPeak None [||] ps nan area 1.
            Expect.isFalse (PeakQuantification.isValidFit (fit [||] 12.5)) "trapezoid fallback: no parameters"
            Expect.isFalse (PeakQuantification.isValidFit (fit [| 1.; 2. |] 3.)) "two parameters"
            Expect.isFalse (PeakQuantification.isValidFit (fit [| 1.; nan; 3. |] 3.)) "NaN parameter"
            Expect.isFalse (PeakQuantification.isValidFit (fit [| 1.; infinity; 3. |] 3.)) "infinite parameter"
            Expect.isFalse (PeakQuantification.isValidFit (fit [| 1.; 2.; 3. |] nan)) "NaN area"
            Expect.isFalse (PeakQuantification.isValidFit (fit [| 1.; 2.; 3. |] infinity)) "infinite area"
            Expect.isTrue (PeakQuantification.isValidFit (fit [| 1.; 2.; 3. |] 3.)) "Gaussian"
            Expect.isTrue (PeakQuantification.isValidFit (fit [| 1.; 2.; 3.; 4. |] 3.)) "EMG"
            let quantified, _ = fittedPeak 5.
            Expect.isTrue (PeakQuantification.isValidFit quantified) "fitted peak"

        testCase "closestPoint equals the peak extraction of PSMBasedQuantificationTIMs" <| fun _ ->
            let rnd = Random 3
            for _ in 1 .. 300 do
                let n = rnd.Next(0, 60)
                // m/z on a grid of 0.001 with repeats (one m/z spread over mobility scans), no ties in the distance to the lock m/z
                let mzs = Array.init n (fun _ -> 500. + float (rnd.Next(0, 400)) * 0.001 + 1e-7)
                let intensities = Array.init n (fun _ -> float (rnd.Next(1, 1000)))
                let ims = Array.init n (fun _ -> if rnd.Next 20 = 0 then nan else 0.8 + rnd.NextDouble() * 0.4)
                let peaks = SortedPeaks.sortPeakArraysByMz mzs intensities ims
                let lock = 500.2 + rnd.NextDouble() * 0.0003
                let halfWindow = 0.05 + rnd.NextDouble() * 0.05
                let imLow, imHigh = 0.9, 1.1
                let intensity, mz, _ = Xic.closestPoint peaks lock halfWindow imLow imHigh
                let expectedIntensity, expectedMz = oldExtract peaks lock halfWindow imLow imHigh
                Expect.equal intensity expectedIntensity "intensity"
                Expect.isTrue (mz = expectedMz || (Double.IsNaN mz && Double.IsNaN expectedMz)) "m/z"

        testCase "closestPoint takes the edges of the window, the lower m/z of two equidistant peaks and spectra without ion mobility" <| fun _ ->
            let spectrum mzs ims =
                let mzs = Array.ofList mzs
                SortedPeaks.sortPeakArraysByMz mzs (Array.init mzs.Length (fun i -> float (i + 1) * 10.)) (Array.ofList ims)
            // a peak exactly at lockMz +- halfWindow is inside the window
            let edge = spectrum [ 500.25; 500.75 ] [ 1.; 1. ]
            let intensity, mz, _ = Xic.closestPoint edge 500.5 0.25 0.5 1.5
            Expect.equal (intensity, mz) (10., 500.25) "lower edge, and the closer of the two edge peaks by file order of equal distance"
            let upper = spectrum [ 500.75 ] [ 1. ]
            let intensity, mz, _ = Xic.closestPoint upper 500.5 0.25 0.5 1.5
            Expect.equal (intensity, mz) (10., 500.75) "upper edge"
            // two peaks with the same distance to the lock m/z: the lower m/z
            let tie = spectrum [ 500.75; 500.25 ] [ 1.; 1. ]
            let intensity, mz, _ = Xic.closestPoint tie 500.5 0.5 0.5 1.5
            Expect.equal mz 500.25 "lower m/z"
            Expect.equal intensity 20. "its intensity"
            // no ion mobility: the closest peak, no window on the ion mobility
            let noIm = spectrum [ 500.1; 500.45; 500.9 ] []
            let intensity, mz, im = Xic.closestPoint noIm 500.5 0.5 0.5 1.5
            Expect.equal (intensity, mz) (20., 500.45) "closest peak"
            Expect.isTrue (Double.IsNaN im) "no ion mobility"

        testCase "extractPoint sums the peaks of the m/z window inside the ion mobility window" <| fun _ ->
            // m/z, intensity, ion mobility
            let rows =
                [ 500.40, 5., 1.00   // outside the m/z window
                  500.48, 10., 1.00
                  500.48, 20., 1.10
                  500.49, 30., 1.00
                  500.50, 40., 0.90
                  500.50, 60., 1.00  // closest m/z, two mobility scans
                  500.50, 1000., 1.50 // outside the mobility window
                  500.52, 50., 1.20
                  500.55, 7., 1.00   // exactly at the edge: inside
                  500.60, 9., 1.00 ] // outside the m/z window
            let peaks =
                SortedPeaks.sortPeakArraysByMz
                    (rows |> List.map (fun (m, _, _) -> m) |> Array.ofList)
                    (rows |> List.map (fun (_, i, _) -> i) |> Array.ofList)
                    (rows |> List.map (fun (_, _, im) -> im) |> Array.ofList)
            let intensity, mz, im = Xic.extractPoint peaks 500.5 0.05 0.95 1.15
            // inside the mobility window (0.95, 1.15): 10, 20, 30, 60, 50 (1.20 is outside), 7 at 500.55
            let expected = [ 500.48, 10., 1.00; 500.48, 20., 1.10; 500.49, 30., 1.00; 500.50, 60., 1.00; 500.55, 7., 1.00 ]
            let total = expected |> List.sumBy (fun (_, i, _) -> i)
            Expect.equal intensity total "summed intensity"
            Expect.floatClose Accuracy.high mz (expected |> List.sumBy (fun (m, i, _) -> m * i) |> fun x -> x / total) "intensity weighted m/z"
            Expect.floatClose Accuracy.high im (expected |> List.sumBy (fun (_, i, v) -> v * i) |> fun x -> x / total) "intensity weighted ion mobility"
            // the closest group alone holds only a part of the signal
            let closest, closestMz, _ = Xic.closestPoint peaks 500.5 0.05 0.95 1.15
            Expect.equal (closest, closestMz) (60., 500.5) "closestPoint keeps the closest m/z"
            // no ion mobility window: all peaks of the m/z window
            let all, _, _ = Xic.extractPoint peaks 500.5 0.05 nan nan
            Expect.equal all (10. + 20. + 30. + 40. + 60. + 1000. + 50. + 7.) "no mobility window"
            // nothing inside
            let none, noneMz, noneIm = Xic.extractPoint peaks 500.5 0.05 2.0 3.0
            Expect.equal none 0. "no peak inside the mobility window"
            Expect.isTrue (Double.IsNaN noneMz && Double.IsNaN noneIm) "m/z and mobility are nan"
            let empty, emptyMz, _ = Xic.extractPoint peaks 700. 0.05 0.95 1.15
            Expect.equal empty 0. "no peak in the m/z window"
            Expect.isTrue (Double.IsNaN emptyMz) "m/z is nan"

        testCase "extractPoint without ion mobility keeps the closest peak" <| fun _ ->
            let spectrum mzs =
                let mzs = Array.ofList mzs
                SortedPeaks.sortPeakArraysByMz mzs (Array.init mzs.Length (fun i -> float (i + 1) * 10.)) [||]
            let peaks = spectrum [ 500.1; 500.45; 500.52; 500.9 ]
            let intensity, mz, im = Xic.extractPoint peaks 500.5 0.1 nan nan
            Expect.equal (intensity, mz) (30., 500.52) "closest peak, no sum"
            Expect.isTrue (Double.IsNaN im) "no ion mobility"
            let sameI, sameMz, _ = Xic.closestPoint peaks 500.5 0.1 0.5 1.5
            Expect.equal (intensity, mz) (sameI, sameMz) "same as closestPoint with an ion mobility window"
            let none, noneMz, _ = Xic.extractPoint peaks 600. 0.1 nan nan
            Expect.equal none 0. "no peak"
            Expect.isTrue (Double.IsNaN noneMz) "m/z is nan"

        testCase "the wavelet peak detection of a short trace, whose kernel reaches past the padded trace, finds the peak" <| fun _ ->
            for n, spacing in [ 10, 0.1; 14, 0.07; 20, 0.05; 25, 0.04; 30, 0.03 ] do
                let xs = Array.init n (fun i -> float i * spacing)
                let center = xs.[n / 2]
                let ys = xs |> Array.map (fun x -> 1000. * exp (-((x - center) ** 2.) / (2. * (2. * spacing) ** 2.)))
                let peaks = FSharpStats'.Wavelet.identify waveletParameters xs ys
                Expect.equal peaks.Length 1 (sprintf "%i points: one peak" n)
                Expect.isLessThan (abs (peaks.[0].Apex.XVal - center)) (2. * spacing) (sprintf "%i points: at the center" n)
            // a trace of 8 points with zeros: no exception
            let xs = Array.init 8 (fun i -> float i * 0.04395373458)
            let ys = [| 0.; 0.; 0.; 999.8592301; 0.; 616.0852004; 933.9307234; 0. |]
            FSharpStats'.Wavelet.identify waveletParameters xs ys |> ignore

        testCase "initIdentifyPeaks fails for Gabor3D only when the function is called" <| fun _ ->
            let gabor : Domain.Gabor3DParams = { sizeX = 3; sizeY = 3; sigmaX = 1.; sigmaY = 1.; frequency = 1.; theta = 0. }
            let identify = PeakQuantification.initIdentifyPeaks (Domain.XicProcessing.Gabor3D gabor)
            Expect.throws (fun () -> identify [| 1.; 2. |] [| 1.; 2. |] |> ignore) "fails when called"

        testCase "calcCorrelation is the correlation of the fits and nan when it cannot be computed" <| fun _ ->
            let light, xs = fittedPeak 5.
            let heavy, _ = fittedPeak 5.2
            let noModel = HULQ.createQuantifiedPeak None [||] [||] nan 1. 1.
            Expect.isTrue (Double.IsNaN (PeakQuantification.calcCorrelation xs noModel heavy)) "no model of the light peak"
            Expect.isTrue (Double.IsNaN (PeakQuantification.calcCorrelation xs light noModel)) "no model of the heavy peak"
            Expect.isTrue (Double.IsNaN (PeakQuantification.calcCorrelation [||] light heavy)) "no scan times"
            Expect.isTrue (Double.IsNaN (PeakQuantification.calcCorrelation [| 5. |] light heavy)) "one scan time"
            Expect.isTrue (Double.IsNaN (PeakQuantification.calcCorrelation [| 5.; 5.; 5. |] light heavy)) "no spacing"
            Expect.floatClose Accuracy.high (PeakQuantification.calcCorrelation xs light light) 1. "a peak with itself"
            // the formula of the unguarded version for a normal XIC
            let getValue (model: HULQ.PeakModel) ps =
                match model with
                | HULQ.PeakModel.Gaussian m -> m.GetFunctionValue (vector ps)
                | HULQ.PeakModel.EMG m -> m.GetFunctionValue (vector ps)
            let bw = [| for i in 1 .. xs.Length - 1 -> abs (xs.[i] - xs.[i - 1]) |] |> Array.median
            let grid = [| xs.[0] .. bw / 2. .. xs.[xs.Length - 1] |]
            let fQ = getValue light.Model.Value light.EstimatedParams
            let fI = getValue heavy.Model.Value heavy.EstimatedParams
            let expected = Correlation.Seq.pearson (grid |> Array.map fQ) (grid |> Array.map fI)
            Expect.equal (PeakQuantification.calcCorrelation xs light heavy) expected "same value"

        testCase "tryPeakContaining takes the peak around the scan time and none when no peak contains it" <| fun _ ->
            let peak left apex right =
                let feature x = Signal.PeakDetection.createPeakFeature -1 x 1.
                Signal.PeakDetection.createIdentifiedPeak (feature apex) None (feature left) None (feature right) false false [| left; apex; right |] [| 1.; 1.; 1. |]
            let first, second = peak 1. 2. 3., peak 5. 6. 7.
            let peaks = [| first; second |]
            Expect.equal (PeakQuantification.tryPeakContaining peaks 2.5) (Some first) "inside the first peak"
            Expect.equal (PeakQuantification.tryPeakContaining peaks 6.5) (Some second) "inside the second peak"
            Expect.isNone (PeakQuantification.tryPeakContaining peaks 4.) "between the peaks, the nearest apex is not taken"
            Expect.isNone (PeakQuantification.tryPeakContaining peaks 1.) "on the left border"
            Expect.isNone (PeakQuantification.tryPeakContaining peaks 3.) "on the right border"
            Expect.isNone (PeakQuantification.tryPeakContaining peaks 0.) "left of all peaks"
            Expect.isNone (PeakQuantification.tryPeakContaining peaks 8.) "right of all peaks"
            Expect.isNone (PeakQuantification.tryPeakContaining [||] 2.) "no peaks"

        testCase "selectPeak requires a containing peak or falls back to the closest apex" <| fun _ ->
            let peak left apex right =
                let feature x = Signal.PeakDetection.createPeakFeature -1 x 1.
                Signal.PeakDetection.createIdentifiedPeak (feature apex) None (feature left) None (feature right) false false [| left; apex; right |] [| 1.; 1.; 1. |]
            let first, second = peak 1. 2. 3., peak 5. 6. 7.
            let peaks = [| first; second |]
            Expect.equal (PeakQuantification.selectPeak true peaks 2.5) (Some first) "required, inside the first peak"
            Expect.equal (PeakQuantification.selectPeak false peaks 2.5) (Some first) "not required, inside the first peak"
            Expect.equal (PeakQuantification.selectPeak false peaks 6.5) (Some second) "not required, inside the second peak"
            Expect.isNone (PeakQuantification.selectPeak true peaks 4.5) "required, between the peaks"
            Expect.equal (PeakQuantification.selectPeak false peaks 4.5) (Some second) "not required, between the peaks the closest apex is taken"
            Expect.isNone (PeakQuantification.selectPeak true peaks 0.) "required, left of all peaks"
            Expect.equal (PeakQuantification.selectPeak false peaks 0.) (Some first) "not required, left of all peaks"
            Expect.isNone (PeakQuantification.selectPeak true peaks 8.) "required, right of all peaks"
            Expect.equal (PeakQuantification.selectPeak false peaks 8.) (Some second) "not required, right of all peaks"
            Expect.isNone (PeakQuantification.selectPeak true [||] 2.) "required, no peaks"
            Expect.isNone (PeakQuantification.selectPeak false [||] 2.) "not required, no peaks"
    ]
