/// Tests of RatioLFQ on peptide ion level: reading the ion table, the pairwise profile of the runs,
/// the reference level and the centring of runs.
module ProteomIQon.UnitTestsRatioLFQ

open System
open System.IO
open Expecto
open ProteomIQon
open ProteomIQon.Domain

let private logTolerance = { absolute = 1e-9; relative = 0. }
let private linearTolerance = { absolute = 0.; relative = 1e-9 }
let private roundedTolerance = { absolute = 0.; relative = 1e-6 }

let private header = "FileName\tProteinGroup\tStringSequence\tQuant_Light\tQuant_Heavy\tCharge"

/// The header and the rows, every ion has charge 2.
let private ionTable (rows: string list) = header :: (rows |> List.map (fun row -> if row = "" then row else row + "\t2"))

/// Reads the lines as an ion table from a temporary file. CsvReader keeps the file open when reading stops with an error.
let private readLines (lines: string list) =
    let path = Path.GetTempFileName()
    try
        File.WriteAllLines(path, lines)
        RatioLFQ.readIonTable path
    finally
        try File.Delete path with :? IOException -> ()

let private parse (rows: string list) = readLines (ionTable rows)

let private parameters reference minShared center : RatioLFQParams =
    {
        Reference      = reference
        MinSharedIons  = minShared
        CentreRuns     = center
        PairWeighting  = PairWeighting.SharedIonCount
        ReferenceLevel = ReferenceLevel.MeanOfIonMedians
    }

/// Checks equal length, equal nan positions and close values.
let private expectClose (tolerance: Accuracy) (expected: float []) (actual: float []) message =
    Expect.equal actual.Length expected.Length (sprintf "%s: length" message)
    Array.iteri2 (fun i (e: float) (a: float) ->
        if Double.IsNaN e then Expect.isTrue (Double.IsNaN a) (sprintf "%s: run %i is nan" message i)
        else Expect.floatClose tolerance a e (sprintf "%s: run %i" message i)) expected actual

/// The worked example: one protein, runs R1 to R3, ions a, b, c.
let private workedExample =
    parse [
        "R1\tP\ta\t2000\t1000"; "R2\tP\ta\t4000\t1000"; "R3\tP\ta\t8000\t1000"
        "R1\tP\tb\t400\t400";   "R2\tP\tb\t1200\t600"
        "R1\tP\tc\t1600\t200";  "R3\tP\tc\t1200\t300"
    ]

let private nanMatrix ions runs = Array.init ions (fun _ -> Array.create runs nan)

/// A table with one protein group per row of ratios (light over heavy, heavy fixed at 1000), one ion each.
let private singleIonTable (ratios: float [] []) =
    [ for p in 0 .. ratios.Length - 1 do
        for r in 0 .. ratios.[p].Length - 1 do
            yield sprintf "R%i\tP%i\tseq\t%s\t1000" r p ((1000. * ratios.[p].[r]).ToString Globalization.CultureInfo.InvariantCulture) ]
    |> parse

let private repoPath (parts: string list) = Path.Combine(Array.ofList (Pipeline.repoRoot :: parts))

/// Runs ratioLFQ on the rows in a temporary directory and returns the lines of RatioLFQ.txt.
let private runRatioLFQ (rows: string list) =
    let work = Path.Combine(Path.GetTempPath(), sprintf "ratioLFQUnit_%s" (Guid.NewGuid().ToString "N"))
    Directory.CreateDirectory work |> ignore
    try
        let table = Path.Combine(work, "ionTable.txt")
        File.WriteAllLines(table, ionTable rows)
        RatioLFQ.ratioLFQ (parameters ReferenceChannel.Heavy 1 false) (Path.Combine(work, "out")) table
        File.ReadAllLines(Path.Combine(work, "out", "RatioLFQ.txt"))
    finally
        Directory.Delete(work, true)

[<Tests>]
let ratioLFQTests =
    testList "RatioLFQ" [
        testCase "worked example" <| fun _ ->
            let check minShared reference ratio intensity =
                let r = (RatioLFQ.compute (parameters reference minShared false) workedExample).[0]
                expectClose roundedTolerance ratio r.Ratio "ratio"
                expectClose roundedTolerance intensity r.Intensity "intensity"
            check 1 ReferenceChannel.Heavy [| 2.181015; 3.363586; 4. |] [| 1272.2590; 1962.0916; 2333.3333 |]
            check 1 ReferenceChannel.Light [| 2.181015; 3.363586; 4. |] [| 947.5708; 614.4237; 516.6667 |]
            check 2 ReferenceChannel.Heavy [| 2.; 4.; 2.828427 |] [| 1166.6667; 2333.3333; 1649.9158 |]
            let r2 = (RatioLFQ.compute (parameters ReferenceChannel.Heavy 2 false) workedExample).[0]
            Expect.equal r2.RunGroups 1 "the pairs with R1 still connect the runs"
            Expect.equal workedExample.Runs [| "R1"; "R2"; "R3" |] "runs"

        testCase "unweighted pairs give the equal weight fit of the worked example" <| fun _ ->
            let unweighted = { parameters ReferenceChannel.Heavy 1 false with PairWeighting = PairWeighting.Unweighted }
            let r = (RatioLFQ.compute unweighted workedExample).[0]
            // z = (-0.5, 0, 0.5), level 1.5, ratios 2^(1, 1.5, 2), reference level 583.3333
            expectClose roundedTolerance [| 2.; 2.828427; 4. |] r.Ratio "ratio"
            expectClose roundedTolerance [| 1166.6667; 1649.9158; 2333.3333 |] r.Intensity "intensity"
            let weighted = (RatioLFQ.compute (parameters ReferenceChannel.Heavy 1 false) workedExample).[0]
            Expect.isGreaterThan (abs (weighted.Ratio.[0] - r.Ratio.[0])) 0.1 "the shared ion count changes the ratio"

        testCase "complete matrix with identical ion profiles gives the plain ratios" <| fun _ ->
            let y = [| [| 1.; 2.5; -1.; 0.3 |]; [| 1.; 2.5; -1.; 0.3 |]; [| 1.; 2.5; -1.; 0.3 |] |]
            expectClose logTolerance (RatioLFQ.plainRatio y) (RatioLFQ.pairwiseProfile 1 y) "x = m"

        testCase "single ion gives the plain ratios" <| fun _ ->
            let y = [| [| 0.5; nan; 2.; -3. |] |]
            let x = RatioLFQ.pairwiseProfile 1 y
            expectClose logTolerance [| 0.5; nan; 2.; -3. |] x "x = m"

        testCase "additive data with missing ions recovers the run differences" <| fun _ ->
            let mu = [| 0.; 1.5; -0.7; 2.2 |]
            let delta = [| 0.; 3.; -2. |]
            let present = [| [| true; true; true; true |]; [| true; true; false; false |]; [| false; true; true; true |] |]
            let y = Array.init 3 (fun p -> Array.init 4 (fun i -> if present.[p].[i] then mu.[i] + delta.[p] else nan))
            let x = RatioLFQ.pairwiseProfile 1 y
            for i in 0 .. 3 do
                for j in 0 .. 3 do
                    Expect.floatClose logTolerance (x.[i] - x.[j]) (mu.[i] - mu.[j]) (sprintf "difference %i %i" i j)

        testCase "two run groups without a shared ion are anchored at their plain ratios" <| fun _ ->
            let table =
                parse [ "R1\tP\ta\t2000\t1000"; "R2\tP\ta\t4000\t1000"
                        "R3\tP\tb\t32000\t1000"; "R4\tP\tb\t128000\t1000" ]
            let r = (RatioLFQ.compute (parameters ReferenceChannel.Heavy 1 false) table).[0]
            Expect.equal r.RunGroups 2 "run groups"
            expectClose linearTolerance [| 2.; 4.; 32.; 128. |] r.Ratio "ratio"
            Expect.isFalse (r.Intensity |> Array.exists Double.IsNaN) "no nan"

        testCase "protein in one run gives the plain ratio and its intensity" <| fun _ ->
            let table =
                parse [ "R1\tP\ta\t2000\t1000"; "R1\tP\tb\t4000\t500"
                        "R2\tQ\ta\t1000\t1000" ]
            let results = RatioLFQ.compute (parameters ReferenceChannel.Heavy 1 false) table
            let p = results.[0]
            // log2 ratios 1 and 3, median 2
            Expect.floatClose linearTolerance p.Ratio.[0] 4. "ratio"
            Expect.isTrue (Double.IsNaN p.Ratio.[1] && Double.IsNaN p.Intensity.[1]) "absent run is nan"
            Expect.equal p.RunGroups 1 "run groups"
            // reference level: mean of the two ion medians 1000 and 500
            Expect.floatClose linearTolerance p.Intensity.[0] (4. * 750.) "intensity"

        testCase "the ratio does not depend on the reference" <| fun _ ->
            let table =
                parse [ "R1\tP\ta\t2000\t1000"; "R2\tP\ta\t4000\t1000"; "R3\tP\ta\t8000\t1000"
                        "R1\tP\tb\t400\t400";   "R2\tP\tb\t1200\t600"
                        "R1\tQ\tc\t1600\t200";  "R3\tQ\tc\t1200\t300" ]
            let heavy = RatioLFQ.compute (parameters ReferenceChannel.Heavy 1 false) table
            let light = RatioLFQ.compute (parameters ReferenceChannel.Light 1 false) table
            Array.iter2 (fun (h: RatioLFQ.ProteinResult) (l: RatioLFQ.ProteinResult) ->
                expectClose logTolerance h.Ratio l.Ratio h.ProteinGroup) heavy light

        testCase "centring gives every run a median log2 ratio of 0" <| fun _ ->
            let table = singleIonTable [| [| 1.; 2. |]; [| 4.; 3. |]; [| 0.5; 8. |]; [| 2.; 1.5 |]; [| 3.; 5. |] |]
            let results = RatioLFQ.compute (parameters ReferenceChannel.Heavy 1 true) table
            for run in 0 .. 1 do
                let med = results |> Array.map (fun r -> Math.Log2 r.Ratio.[run]) |> Numerics.median
                Expect.floatClose logTolerance med 0. (sprintf "median of run %i" run)

        testCase "referenceLevel applies both rules to the ion medians" <| fun _ ->
            let reference = [| [| 10.; 30.; nan |]; [| 100.; nan; 300. |]; [| 5.; 5.; 5. |]; [| 1000.; 2000.; nan |]; [| nan; nan; nan |] |]
            // ion medians 20, 200, 5, 1500: mean 1725 / 4 = 431.25, the three largest are 1500, 200 and 20 with mean 1720 / 3
            let level rule = RatioLFQ.referenceLevel rule reference
            Expect.floatClose linearTolerance (level ReferenceLevel.MeanOfIonMedians) 431.25 "mean"
            Expect.floatClose linearTolerance (level ReferenceLevel.TopThreeIonMedians) (1720. / 3.) "top three"
            Expect.floatClose linearTolerance (RatioLFQ.referenceLevel ReferenceLevel.TopThreeIonMedians [| [| 10.; 30. |]; [| 100.; nan |] |]) 60. "two ions"
            Expect.floatClose linearTolerance (RatioLFQ.referenceLevel ReferenceLevel.MeanOfIonMedians [| [| 7.; nan |] |]) 7. "one ion"

        testCase "the reference level rule changes the intensities and not the ratios" <| fun _ ->
            // the worked example plus ion d with heavy 5000 in R1 and R3: heavy ion medians 1000, 500, 250 and 5000,
            // so the mean is 1687.5 and the three largest give 6500 / 3
            let table =
                parse [
                    "R1\tP\ta\t2000\t1000"; "R2\tP\ta\t4000\t1000"; "R3\tP\ta\t8000\t1000"
                    "R1\tP\tb\t400\t400";   "R2\tP\tb\t1200\t600"
                    "R1\tP\tc\t1600\t200";  "R3\tP\tc\t1200\t300"
                    "R1\tP\td\t5000\t5000"; "R3\tP\td\t10000\t5000"
                ]
            let result rule = (RatioLFQ.compute { parameters ReferenceChannel.Heavy 1 false with ReferenceLevel = rule } table).[0]
            let mean, top = result ReferenceLevel.MeanOfIonMedians, result ReferenceLevel.TopThreeIonMedians
            expectClose logTolerance mean.Ratio top.Ratio "ratio"
            let differs = Array.exists2 (fun (m: float) (t: float) -> abs (m - t) > 1e-6 * abs m) mean.Intensity top.Intensity
            Expect.isTrue differs "the intensities differ"
            Array.iteri2 (fun i (r: float) (v: float) -> Expect.floatClose linearTolerance v (r * 1687.5) (sprintf "mean intensity run %i" i)) mean.Ratio mean.Intensity
            Array.iteri2 (fun i (r: float) (v: float) -> Expect.floatClose linearTolerance v (r * 6500. / 3.) (sprintf "top three intensity run %i" i)) top.Ratio top.Intensity

        testCase "readIonTable reads and rejects as specified" <| fun _ ->
            Expect.throws (fun () -> readLines [ "FileName\tProteinGroup\tStringSequence\tQuant_Light"; "R1\tP\ta\t1" ] |> ignore) "missing column"
            let message = (try readLines [ "FileName\tProteinGroup\tStringSequence\tQuant_Light" ] |> ignore; "" with e -> e.Message)
            Expect.stringContains message "Quant_Heavy" "the message names the column"
            let duplicate = (try parse [ "R1\tGROUPX\tDUPLICATEDK\t1000\t1000"; "R1\tGROUPX\tDUPLICATEDK\t2000\t1000" ] |> ignore; "" with e -> e.Message)
            Expect.stringContains duplicate "R1" "duplicate names the run"
            Expect.stringContains duplicate "DUPLICATEDK" "duplicate names the ion"
            Expect.stringContains duplicate "GROUPX" "duplicate names the protein group"
            // the same ion in two runs is fine
            let plain = parse [ "R1\tP\ta\t1000\t1000"; "R2\tP\ta\t2000\t1000" ]
            Expect.equal plain.Measurements.[0].Length 2 "one ion in two runs"
            // the charge separates ions
            let charged = readLines [ header; "R1\tP\ta\t1000\t1000\t2"; "R1\tP\ta\t1000\t1000\t3" ]
            Expect.equal (charged.Measurements.[0] |> Array.map (fun m -> m.Ion)) [| 0; 1 |] "two charges, two ions"
            // invalid rows: zero, negative, NaN, empty intensity, empty protein group
            let withInvalid =
                parse [ "R1\tP\ta\t0\t1000"; "R1\tP\tb\t-5\t1000"; "R1\tP\tc\tNaN\t1000"; "R1\tP\td\t\t1000"
                        "R1\t\te\t1000\t1000"; "R2\tP\tf\t1000\t1000"; ""; "R3\tOnlyInvalid\tg\t0\t0" ]
            Expect.equal withInvalid.RowsRead 7 "blank line is not a row"
            Expect.equal withInvalid.RowsIgnored 6 "ignored rows"
            Expect.equal withInvalid.ProteinGroups [| "P" |] "a protein group with only invalid rows is absent"
            Expect.equal withInvalid.Runs [| "R1"; "R2"; "R3" |] "runs by first appearance over all rows"
            Expect.equal withInvalid.Measurements.[0].[0].Run 1 "the valid row is in R2"
            let empty = parse []
            Expect.equal (empty.Runs, empty.ProteinGroups, empty.Measurements.Length, empty.RowsRead) ([||], [||], 0, 0) "header only"
            // quoted run names and protein groups lose their quotes
            let quoted =
                parse [ "\"sample 03\"\t\"Prot07;Prot08\"\ta\t1000\t1000"; "R2\tProt09\tb\t1000\t1000" ]
            Expect.equal quoted.Runs [| "sample 03"; "R2" |] "quoted run name with a space"
            Expect.equal quoted.ProteinGroups [| "Prot07;Prot08"; "Prot09" |] "quoted protein groups"
            // one sequence under two protein groups in one run belongs to both
            let shared = parse [ "R1\tG1\tSHAREDK\t1000\t1000"; "R1\tG2\tSHAREDK\t2000\t1000" ]
            Expect.equal shared.ProteinGroups [| "G1"; "G2" |] "both protein groups"
            Expect.equal (shared.Measurements |> Array.map Array.length) [| 1; 1 |] "one measurement in each protein group"

        testCase "a table without a valid measurement gives the header only" <| fun _ ->
            let lines = runRatioLFQ [ "R1\tP\ta\t0\t1000"; "R2\tQ\tb\t1000\t" ]
            Expect.equal lines [| "ProteinGroup\tR1.Ratio_LightByHeavy\tR1.LFQ_Intensity\tR2.Ratio_LightByHeavy\tR2.LFQ_Intensity" |] "header only"

        testCase "a run without a valid cell keeps its columns with empty cells" <| fun _ ->
            let lines = runRatioLFQ [ "R1\tP\ta\t2000\t1000"; "R2\tP\ta\t0\t1000"; "R3\tP\ta\t4000\t1000" ]
            Expect.equal lines.[0] "ProteinGroup\tR1.Ratio_LightByHeavy\tR1.LFQ_Intensity\tR2.Ratio_LightByHeavy\tR2.LFQ_Intensity\tR3.Ratio_LightByHeavy\tR3.LFQ_Intensity" "header"
            let cells = lines.[1].Split '\t'
            Expect.equal (cells.[0], cells.[3], cells.[4]) ("P", "", "") "R2 cells are empty"
            let value i = Double.Parse(cells.[i], Globalization.CultureInfo.InvariantCulture)
            expectClose linearTolerance [| 2.; 2000.; 4.; 4000. |] [| value 1; value 2; value 5; value 6 |] "R1 and R3"

        testCase "a table with log2 light and linear heavy intensities is rejected as transformed" <| fun _ ->
            let work = Path.Combine(Path.GetTempPath(), sprintf "ratioLFQUnit_%s" (Guid.NewGuid().ToString "N"))
            Directory.CreateDirectory work |> ignore
            try
                let table = Path.Combine(work, "ionTable.txt")
                File.WriteAllLines(table, ionTable [ "R1\tP\ta\t12.5\t5000"; "R2\tP\ta\t13.1\t6000" ])
                let message = (try RatioLFQ.ratioLFQ (parameters ReferenceChannel.Heavy 1 false) (Path.Combine(work, "out")) table; "" with e -> e.Message)
                Expect.stringContains message "transformed" "the message says transformed"
                Expect.isFalse (File.Exists(Path.Combine(work, "out", "RatioLFQ.txt"))) "no output file"
            finally
                Directory.Delete(work, true)

        testCase "compute rejects MinSharedIons 0" <| fun _ ->
            Expect.throws (fun () -> RatioLFQ.compute (parameters ReferenceChannel.Heavy 0 false) workedExample |> ignore) "MinSharedIons = 0"

        testCase "a one-run table written by LabeledProteinQuantification gives the median log2 ratio of each protein group" <| fun _ ->
            let path = repoPath [ "tests"; "data"; "LabeledProteinQuantification"; "out"; "corrFilterChargeAgg"; "minimalGlobModReference.txt" ]
            let lines = File.ReadAllLines path
            let table = RatioLFQ.readIonTable path
            Expect.equal table.Runs.Length 1 "one run"
            let columns = lines.[0].Split '\t'
            let iGroup, iLight, iHeavy = Array.IndexOf(columns, "ProteinGroup"), Array.IndexOf(columns, "Quant_Light"), Array.IndexOf(columns, "Quant_Heavy")
            let unquote (c: string) = if c.Length >= 2 && c.StartsWith "\"" && c.EndsWith "\"" then c.[1 .. c.Length - 2].Replace("\"\"", "\"") else c
            let tryNumber (c: string) =
                match Double.TryParse(c, Globalization.NumberStyles.Float, Globalization.CultureInfo.InvariantCulture) with
                | true, v when v > 0. && not (Double.IsInfinity v) -> Some v
                | _ -> None
            let expected =
                lines.[1..]
                |> Array.map (fun l -> l.Split '\t' |> Array.map unquote)
                |> Array.filter (fun c -> (tryNumber c.[iLight]).IsSome && (tryNumber c.[iHeavy]).IsSome)
                |> Array.groupBy (fun c -> c.[iGroup])
                |> Array.map (fun (g, rows) -> g, rows |> Array.map (fun c -> Math.Log2((tryNumber c.[iLight]).Value / (tryNumber c.[iHeavy]).Value)) |> Numerics.median |> fun m -> Math.Pow(2., m))
            let results = RatioLFQ.compute (parameters ReferenceChannel.Heavy 1 false) table
            Expect.equal (results |> Array.map (fun r -> r.ProteinGroup)) (expected |> Array.map fst) "protein groups"
            Array.iter2 (fun (r: RatioLFQ.ProteinResult) (_, e) -> Expect.floatClose linearTolerance r.Ratio.[0] e r.ProteinGroup) results expected

        testCase "the charge aggregated table of LabeledProteinQuantification is rejected" <| fun _ ->
            let path = repoPath [ "tests"; "data"; "LabeledProteinQuantification"; "out"; "corrFilterChargeAgg"; "minimalChargeReference.txt" ]
            let message = (try RatioLFQ.readIonTable path |> ignore; "" with e -> e.Message)
            Expect.stringContains message "column(s) Charge" "the message names the column"
    ]
