/// Tests of the quality filter shared by the two PSM quantification tools.
module ProteomIQon.UnitTestsPSMQuantificationQualityFilter

open System
open Expecto
open ProteomIQon

/// A result with the given id (PepSequenceID) and light and heavy area and apex; all other fields have their default value.
let private row (id: int) (quantLight: float, apexLight: float) (quantHeavy: float, apexHeavy: float) : Dto.QuantificationResult =
    let values = dict [ "PepSequenceID", box id; "Quant_Light", box quantLight; "MeasuredApex_Light", box apexLight; "Quant_Heavy", box quantHeavy; "MeasuredApex_Heavy", box apexHeavy ]
    let fields = Reflection.FSharpType.GetRecordFields typeof<Dto.QuantificationResult>
    let defaultOf (t: Type) = if t.IsValueType then Activator.CreateInstance t else null
    Reflection.FSharpValue.MakeRecord (typeof<Dto.QuantificationResult>, fields |> Array.map (fun f -> match values.TryGetValue f.Name with | true, v -> v | _ -> defaultOf f.PropertyType)) :?> Dto.QuantificationResult

/// The ids of the rows the light filter keeps with the borders -1 and 1.
let private keptLight (rows: Dto.QuantificationResult []) =
    rows |> PSMQuantification.lightQualityFilter -1. 1. |> Array.map (fun r -> r.PepSequenceID)

let private keptHeavy (rows: Dto.QuantificationResult []) =
    rows |> PSMQuantification.heavyQualityFilter -1. 1. |> Array.map (fun r -> r.PepSequenceID)

[<Tests>]
let qualityFilterTests =
    let good = (10., 10.)
    testList "PSMQuantification quality filter" [
        testCase "rows inside the borders are kept, rows outside are dropped" <| fun () ->
            let rows =
                [| row 1 good good; row 2 good good; row 3 good good
                   row 4 (15., 10.) (15., 10.)     // log2 1.5 inside
                   row 5 (10., 40.) (10., 40.)     // log2 0.25 below
                   row 6 (40., 10.) (40., 10.) |]  // log2 4 above
            Expect.equal (keptLight rows) [| 1; 2; 3; 4 |] "light"
            Expect.equal (keptHeavy rows) [| 1; 2; 3; 4 |] "heavy"
        testCase "a row with an empty side is dropped by the filter of that side only" <| fun () ->
            let rows = [| row 1 good good; row 2 good good; row 3 good (nan, nan); row 4 (nan, nan) good; row 5 good (10., nan); row 6 (nan, 10.) good |]
            Expect.equal (keptLight rows) [| 1; 2; 3; 5 |] "light"
            Expect.equal (keptHeavy rows) [| 1; 2; 4; 6 |] "heavy"
        testCase "a row with a non-positive area or apex is dropped" <| fun () ->
            let bad = [ (0., 10.); (-5., 10.); (10., 0.); (10., -1.); (infinity, 10.); (10., infinity) ]
            let rows = Array.append [| row 1 good good; row 2 good good; row 3 good good |] (bad |> List.mapi (fun i b -> row (10 + i) b b) |> Array.ofList)
            Expect.equal (keptLight rows) [| 1; 2; 3 |] "light"
            Expect.equal (keptHeavy rows) [| 1; 2; 3 |] "heavy"
        testCase "empty input gives no rows" <| fun () ->
            Expect.isEmpty (keptLight [||]) "light"
            Expect.isEmpty (keptHeavy [||]) "heavy"
        testCase "no finite, positive median gives no rows" <| fun () ->
            let rows = [| row 1 (nan, nan) (nan, nan); row 2 (0., 0.) (-1., -1.); row 3 (nan, 10.) (10., nan) |]
            Expect.isEmpty (keptLight rows) "light"
            Expect.isEmpty (keptHeavy rows) "heavy"
        testCase "the medians ignore empty and non-positive values" <| fun () ->
            // with the usable areas 10, 10, 10, 16 the median is 10 and the row with 16 is inside (log2 1.6 < 1);
            // a median over all areas would be 5 and put it outside (log2 3.2)
            let rows = [| row 1 good good; row 2 good good; row 3 good good; row 4 (16., 10.) (16., 10.)
                          row 5 (0., 10.) (0., 10.); row 6 (0., 10.) (0., 10.); row 7 (0., 10.) (0., 10.); row 8 (-1., 10.) (-1., 10.)
                          row 9 (nan, 10.) (nan, 10.) |]
            Expect.equal (keptLight rows) [| 1; 2; 3; 4 |] "light"
            Expect.equal (keptHeavy rows) [| 1; 2; 3; 4 |] "heavy"
        ]

[<Tests>]
let labelingQualityFilterTests =
    let good = (10., 10.)
    // row 4 is outside on the light side only, row 5 on the heavy side only
    let rows = [| row 1 good good; row 2 good good; row 3 good good; row 4 (40., 10.) good; row 5 good (40., 10.) |]
    let kept labeling = PSMQuantification.labelingQualityFilter labeling 1. rows |> Array.map (fun r -> r.PepSequenceID)
    testList "PSMQuantification quality filter per labeling" [
        testCase "N15Labeling filters the heavy and the light side" <| fun () ->
            Expect.equal (kept Domain.Labeling.N15Labeling) [| 1; 2; 3 |] "kept"
        testCase "N15LabelingOnly filters the heavy side" <| fun () ->
            Expect.equal (kept Domain.Labeling.N15LabelingOnly) [| 1; 2; 3; 4 |] "kept"
        testCase "Unlabeled filters the light side" <| fun () ->
            Expect.equal (kept Domain.Labeling.Unlabeled) [| 1; 2; 3; 5 |] "kept"
        testCase "Labelshift keeps all rows" <| fun () ->
            Expect.equal (kept Domain.Labeling.Labelshift) [| 1; 2; 3; 4; 5 |] "kept"
        testCase "N15Labeling filters the heavy side first, so the light medians come from the rows it keeps" <| fun () ->
            // the heavy filter drops row 3, the light areas of rows 1, 2 and 4 then have the median 10 and put row 4 (25) outside
            // (ratio 2.5); the light filter first would have the median 17.5 over all four rows and keep row 4 (ratio 1.43)
            let rows = [| row 1 good good; row 2 good good; row 3 (40., 10.) (40., 10.); row 4 (25., 10.) good |]
            Expect.equal (PSMQuantification.labelingQualityFilter Domain.Labeling.N15Labeling 1. rows |> Array.map (fun r -> r.PepSequenceID)) [| 1; 2 |] "kept"
            Expect.equal (rows |> PSMQuantification.lightQualityFilter -1. 1. |> PSMQuantification.heavyQualityFilter -1. 1. |> Array.map (fun r -> r.PepSequenceID)) [| 1; 2; 4 |] "the other order keeps row 4"
        ]
