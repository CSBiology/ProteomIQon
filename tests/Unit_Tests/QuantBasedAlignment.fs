/// Tests of the ion mobility QuantBasedAlignment aligns.
module ProteomIQon.UnitTestsQuantBasedAlignment

open System
open System.IO
open Expecto
open FSharpAux.IO.SchemaReader
open ProteomIQon
open ProteomIQon.Dto

[<Tests>]
let tests =
    testList "QuantBasedAlignment ion mobility" [

        testCase "the ion mobility model, the .align rows and the .alignmetric rows use the ion mobility centre of the rows" <| fun _ ->
            // in copies of the test .quant files the ion mobility of every row is 0.9. The centre is 0.8 plus a small value that depends
            // on the ion in the source runs (trunc2, trunc3) and 1.0 plus the same value in the target run (trunc1), a shift of 0.2,
            // except for the rows of charge 3, which get no centre
            let directory = Path.Combine (Path.GetTempPath (), sprintf "qba_ionmobility_%s" (Guid.NewGuid().ToString "N"))
            let inDirectory = Path.Combine (directory, "in")
            let outDirectory = Path.Combine (directory, "out")
            Directory.CreateDirectory inDirectory |> ignore
            Directory.CreateDirectory outDirectory |> ignore
            let files =
                [| "trunc1", 1.0; "trunc2", 0.8; "trunc3", 0.8 |]
                |> Array.map (fun (name, centre) ->
                    let lines = File.ReadAllLines (Path.Combine (__SOURCE_DIRECTORY__, "..", "data", "QuantBasedAlignment", "in", name + ".quant"))
                    let header = lines.[0].Split '\t'
                    let column (name: string) = Array.IndexOf (header, name)
                    let changed =
                        lines
                        |> Array.mapi (fun i line ->
                            if i = 0 then line
                            else
                                let fields = line.Split '\t'
                                fields.[column "IonMobility"] <- "0.9"
                                let ion = (fields.[column "PepSequenceID"] |> int) % 10
                                fields.[column "IonMobilityCentre"] <- (if fields.[column "Charge"] = "3" then "NaN" else sprintf "%.4f" (centre + 0.001 * float ion))
                                String.Join ("\t", fields))
                    Expect.isTrue (changed |> Array.exists (fun line -> line.Split('\t').[column "Charge"] = "3")) (sprintf "%s has rows of charge 3, which have a finite ion mobility and no centre" name)
                    let path = Path.Combine (inDirectory, name + ".quant")
                    File.WriteAllLines (path, changed)
                    path)
            try
                QuantBasedAlignment.alignRuns false outDirectory 1 1 [| files.[0] |] [| files.[1]; files.[2] |]
                let rows =
                    Csv.CsvReader<AlignmentResult>().ReadFile(Path.Combine (outDirectory, "trunc1.align"), '\t', false, 1) |> Array.ofSeq
                Expect.isTrue (rows.Length > 0) "transferred ions"
                for row in rows do
                    if row.Charge = 3 then
                        Expect.isTrue (Double.IsNaN row.IonMobility) "a row without a centre has no ion mobility"
                        Expect.isTrue (Double.IsNaN row.PredictedIonMobility) "and no predicted ion mobility"
                    else
                        // the source centre, not the ion mobility 0.9 of the rows and not a centre of the target run
                        Expect.isTrue (row.IonMobility >= 0.8 && row.IonMobility < 0.81) (sprintf "the centre of the source row: %g" row.IonMobility)
                let predicted = rows |> Array.filter (fun r -> Double.IsFinite r.PredictedIonMobility)
                Expect.isTrue (predicted.Length > 0) "ion mobility predictions from the centres"
                for row in predicted do
                    // the centre of the source plus the shift to the target run, so near the centres of the target run
                    Expect.isTrue (abs (row.PredictedIonMobility - 1.0) < 0.02) (sprintf "the prediction lies at the centres of the target run: %g" row.PredictedIonMobility)
                let metrics =
                    Csv.CsvReader<AlignmentMetricsDTO>().ReadFile(Path.Combine (outDirectory, "trunc1.alignmetric"), '\t', false, 1) |> Array.ofSeq
                Expect.isTrue (metrics.Length > 0) "metric rows"
                for row in metrics do
                    if row.Charge = 3 then Expect.isTrue (Double.IsNaN row.IonMobility) "no ion mobility without a centre in the metric rows"
                    else Expect.isTrue (row.IonMobility >= 0.8 && row.IonMobility < 0.81) (sprintf "the centre of the source row in the metric rows: %g" row.IonMobility)
            finally
                Directory.Delete (directory, true)
    ]
