/// Tests of the readers of the spectrum descriptions and the sorted-array helpers the quantification and search tools share.
module ProteomIQon.UnitTestsSpectrumIndex

open System
open System.Data.SQLite
open Expecto
open ProteomIQon

let private runId = "sample=0"

let private description id msLevel (scanTime: string option) =
    let scans =
        match scanTime with
        | Some t -> sprintf ",\"Scans\":{\"properties\":{\"s\":{\"properties\":{\"MS:1000016\":{\"Values\":[\"%s\"]}}}}}" t
        | None -> ""
    sprintf "{\"ID\":\"%s\",\"properties\":{\"MS:1000511\":{\"Values\":[\"%i\"]}}%s}" id msLevel scans

/// A copy of the MS2 only fixture with the given descriptions added as spectra.
let private withDescriptions (descriptions: (string * string) list) (f: string -> unit) =
    let path = IO.Path.Combine(IO.Path.GetTempPath(), sprintf "spectrumIndex_%s.mzlite" (Guid.NewGuid().ToString "N"))
    IO.File.Copy(IO.Path.Combine(AppContext.BaseDirectory, "minimalTIMs.mzlite"), path)
    try
        do
            use connection = new SQLiteConnection(sprintf "Data Source=%s" path)
            connection.Open()
            for (spectrumId, json) in descriptions do
                use cmd = new SQLiteCommand("INSERT INTO Spectrum (RunID, SpectrumID, Description, PeakArray, PeakData) VALUES (@run, @id, @description, '{}', x'00')", connection)
                cmd.Parameters.AddWithValue("@run", runId) |> ignore
                cmd.Parameters.AddWithValue("@id", spectrumId) |> ignore
                cmd.Parameters.AddWithValue("@description", json) |> ignore
                cmd.ExecuteNonQuery() |> ignore
        f path
    finally
        SQLiteConnection.ClearAllPools()
        try IO.File.Delete path with _ -> ()

let private withReader (descriptions: (string * string) list) (f: MzIO.MzSQL.MzSQL -> unit) =
    withDescriptions descriptions (fun path ->
        use reader = Core.MzIO.Reader.getReader path
        Core.MzIO.Reader.openConnection reader
        f (reader :?> MzIO.MzSQL.MzSQL))

let private goodMs1 =
    [
        "ms1_c", description "ms1_c" 1 (Some "3.0")
        "ms1_a", description "ms1_a" 1 (Some "1.0")
        "ms1_b", description "ms1_b" 1 (Some "2.0")
        "ms3_extra", description "ms3_extra" 3 (Some "2.5")
    ]

/// Descriptions of spectra the tools cannot use: not JSON, no scan time, NaN and infinite scan time, blank ID, and one that is not a MS1.
let private badRows =
    [
        "bad_json", "{\"ID\":\"bad_json\",\"properties\":{"
        "bad_missing", description "bad_missing" 1 None
        "bad_nan", description "bad_nan" 1 (Some "NaN")
        "bad_inf", description "bad_inf" 1 (Some "Infinity")
        "bad_blank", description " " 1 (Some "4.0")
        "bad_ms2_json", "not a json"
    ]

[<Tests>]
let all =
    testList "SpectrumIndex" [
        testCase "the MS1 index is sorted by scan time and hands the other spectra to the callback" <| fun _ ->
            withReader goodMs1 (fun reader ->
                let others = ResizeArray<string>()
                let warnings = ResizeArray<string>()
                let index = SpectrumIndex.readMs1Index reader runId warnings.Add (fun h -> others.Add h.ID)
                Expect.equal index [| 1.0, "ms1_a"; 2.0, "ms1_b"; 3.0, "ms1_c" |] "sorted by scan time"
                Expect.isEmpty warnings "no warning for clean descriptions"
                Expect.equal others.Count 63 "62 MS2 spectra of the fixture and the added MS3"
                Expect.contains others "ms3_extra" "other spectra are given to the callback")
        testCase "descriptions that cannot be used are skipped with one warning instead of failing" <| fun _ ->
            withReader (goodMs1 @ badRows) (fun reader ->
                let warnings = ResizeArray<string>()
                let index = SpectrumIndex.readMs1Index reader runId warnings.Add ignore
                Expect.equal index [| 1.0, "ms1_a"; 2.0, "ms1_b"; 3.0, "ms1_c" |] "only the usable MS1 spectra"
                Expect.equal warnings.Count 1 "one warning"
                Expect.stringContains warnings.[0] "6 spectrum descriptions" "the warning counts the skipped spectra")
        testCase "the retention time index of PSMBasedQuantification skips bad descriptions" <| fun _ ->
            withReader (goodMs1 @ badRows) (fun reader ->
                let warnings = ResizeArray<string>()
                let rtIndex, ms1s = PSMBasedQuantification.Query.readMs1Index reader runId warnings.Add
                Expect.equal ms1s [| 1.0, "ms1_a"; 2.0, "ms1_b"; 3.0, "ms1_c" |] "only the usable MS1 spectra"
                Expect.equal (Seq.length rtIndex) 3 "the retention time index has them"
                Expect.equal warnings.Count 1 "one warning")
        testCase "the spectrum index of PSMBasedQuantificationTIMs skips bad descriptions and keeps the ion mobility of the PSM spectra" <| fun _ ->
            withReader (goodMs1 @ badRows) (fun reader ->
                let warnings = ResizeArray<string>()
                let psmSpectrumIds = Collections.Generic.HashSet<string>([ "ms3_extra"; "bad_json" ])
                let rtIndex, ms1s, mobilities = PSMBasedQuantificationTIMs.QpsmInput.readSpectrumIndex reader runId warnings.Add psmSpectrumIds
                Expect.equal ms1s [| 1.0, "ms1_a"; 2.0, "ms1_b"; 3.0, "ms1_c" |] "only the usable MS1 spectra"
                Expect.equal (Seq.length rtIndex) 3 "the retention time index has them"
                Expect.equal (mobilities.Keys |> Seq.toArray) [| "ms3_extra" |] "the readable PSM spectrum, without ion mobility"
                Expect.isTrue (Double.IsNaN mobilities.["ms3_extra"]) "NaN without ion mobility"
                Expect.equal warnings.Count 1 "one warning")
        testCase "the MS2 headers skip descriptions that cannot be read and keep the order of the fixture" <| fun _ ->
            withReader [ "bad_ms2_json", "not a json"; "ms1_a", description "ms1_a" 1 (Some "1.0") ] (fun reader ->
                let log = ResizeArray<string>()
                let headers = SpectrumIndex.readMs2Headers reader runId log.Add
                Expect.equal headers.Length 62 "the MS2 spectra of the fixture"
                Expect.equal log.Count 1 "the bad description is logged"
                Expect.isTrue (headers |> Array.pairwise |> Array.forall (fun (a, b) -> (a.ScanTime, a.Position) <= (b.ScanTime, b.Position))) "scan time order, file order on ties")
        testCase "closestIndex equals the first minimum of Array.minBy, also with equal scan times and equal distances" <| fun _ ->
            let rng = Random 11
            for _ in 1 .. 300 do
                let times = Array.init (1 + rng.Next 12) (fun _ -> float (rng.Next 10)) |> Array.sort
                for target in [ -3.; 0.; 0.5; 1.; 2.5; 4.5; 9.5; 12. ] @ [ for _ in 1 .. 5 -> float (rng.Next 20) / 2. - 2. ] do
                    let expected = times |> Array.mapi (fun i t -> i, t) |> Array.minBy (fun (_, t) -> abs (t - target)) |> fst
                    Expect.equal (SortedSearch.closestIndex times target) expected (sprintf "%A %f" times target)
        testCase "initClosestMs1 gives the item of the closest scan time" <| fun _ ->
            let closest = SpectrumIndex.initClosestMs1 [| 1.0, "a"; 2.0, "b"; 2.0, "c"; 4.0, "d" |]
            Expect.equal (closest 0.) "a" "before the first"
            Expect.equal (closest 1.5) "a" "equal distance: the earlier one"
            Expect.equal (closest 2.) "b" "the first of equal scan times"
            Expect.equal (closest 3.) "b" "equal distance to 2 and 4: the earlier one"
            Expect.equal (closest 100.) "d" "after the last"
        testCase "the m/z sort drops peaks without m/z and keeps the file order of equal m/z" <| fun _ ->
            let sorted = TimsPeaks.sortPeakArraysByMz [| 5.; nan; 3.; 5.; 3.; 4. |] [| 10.; 20.; 30.; 40.; 50.; 60. |] [| 0.1; 0.2; 0.3; 0.4; 0.5; 0.6 |]
            Expect.equal sorted.Mz [| 3.; 3.; 4.; 5.; 5. |] "sorted by m/z, NaN dropped"
            Expect.equal sorted.Index [| 2; 4; 5; 0; 3 |] "equal m/z in file order"
            Expect.equal sorted.Intensity [| 30.; 50.; 60.; 10.; 40. |] "intensities follow"
            Expect.equal sorted.IonMobility [| 0.3; 0.5; 0.6; 0.1; 0.4 |] "ion mobilities follow"
            Expect.isEmpty (TimsPeaks.sortPeakArraysByMz [| 2.; 1. |] [| 1.; 2. |] [||]).IonMobility "no ion mobility stays empty"
        testCase "mapiOrderedBy runs in the order of the key and returns the results in the order of the input" <| fun _ ->
            let calls = ResizeArray<int * float>()
            let keys = [| 3.; nan; 1.; 3.; 2. |]
            let results = PSMBasedQuantificationTIMs.mapiOrderedBy id (fun position key -> calls.Add((position, key)); sprintf "%i:%g" position key) keys
            Expect.equal (calls |> Seq.map fst |> Seq.toList) [ 0; 1; 2; 3; 4 ] "the positions in the order of the key"
            Expect.equal (calls |> Seq.map snd |> Seq.toList |> List.map (fun k -> if Double.IsNaN k then -1. else k)) [ -1.; 1.; 2.; 3.; 3. ] "NaN first, equal keys in the order of the input"
            Expect.equal results [| "3:3"; "0:NaN"; "1:1"; "4:3"; "2:2" |] "the results stay in the order of the input"
            Expect.isEmpty (PSMBasedQuantificationTIMs.mapiOrderedBy id (fun _ _ -> 0) [||]) "nothing to map"
        testCase "the peak cache keeps at most its capacity, drops the least recently used and counts the reads" <| fun _ ->
            let reads = ResizeArray<string>()
            let read (id: string) =
                reads.Add id
                TimsPeaks.sortPeakArraysByMz [| float id.Length |] [| 1. |] [||]
            let cache = PSMBasedQuantificationTIMs.Query.PeakCache(read, 10)
            for i in 1 .. 10 do cache.Get (string (i * 11)) |> ignore
            Expect.equal cache.Reads 10 "ten different spectra"
            cache.Get "11" |> ignore
            Expect.equal cache.Reads 10 "a hit is not a read"
            // the tenth that was used least recently is dropped when the cache is full: "22"
            cache.Get "x" |> ignore
            Expect.equal cache.Reads 11 "a new spectrum"
            cache.Get "11" |> ignore
            Expect.equal cache.Reads 11 "recently used spectra stay"
            cache.Get "22" |> ignore
            Expect.equal cache.Reads 12 "the least recently used one was dropped and is read again"
            Expect.equal (reads |> Seq.filter ((=) "22") |> Seq.length) 2 "read twice"
    ]
