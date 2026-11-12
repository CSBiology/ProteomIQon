/// Tests of the encoding of the isotope patterns, their table in the peptide data base, the look up of the quantification tools and the store of PSMStatistics.
module ProteomIQon.UnitTestsIsotopePatterns

open System
open System.Data.SQLite
open System.IO
open Expecto
open BioFSharp.Mz
open ProteomIQon

let private dataDirectory = Path.Combine(Pipeline.repoRoot, "tests", "data")

/// The data base with isotopic modifications of the PSMBasedQuantification test and its ions, and the data base without of the PSMStatistics test.
let private labelledDb = Path.Combine(dataDirectory, "PSMBasedQuantification", "in", "Minimal.db")
let private labelledQpsm = Path.Combine(dataDirectory, "PSMBasedQuantification", "in", "minimal.qpsm")
let private unlabelledDb = Path.Combine(dataDirectory, "PSMStatistics", "in", "MinimalEstimate.db")
let private unlabelledQpsm = Path.Combine(dataDirectory, "PSMStatistics", "out", "estimateOut", "minimalReference.qpsm")

let private readIons = Pipeline.readIons

/// A logger that keeps its messages as "LEVEL message".
let private memoryLogger () =
    let target = new NLog.Targets.MemoryTarget("memory")
    target.Layout <- NLog.Layouts.Layout.FromString "${level:uppercase=true} ${message}"
    let config = NLog.Config.LoggingConfiguration()
    config.AddRule(NLog.LogLevel.Trace, NLog.LogLevel.Fatal, target)
    let factory = new NLog.LogFactory()
    factory.Configuration <- config
    factory.GetLogger "IsotopePatternsTest", target

/// Calls f with the path of a temporary copy of the data base.
let private withCopy (source: string) (f: string -> unit) =
    let path = Path.Combine(Path.GetTempPath(), sprintf "isotopePatterns_%s.db" (Guid.NewGuid().ToString "N"))
    File.Copy(source, path)
    try
        f path
    finally
        SQLiteConnection.ClearAllPools()
        GC.Collect()
        GC.WaitForPendingFinalizers()
        try
            File.SetAttributes(path, FileAttributes.Normal)
            File.Delete path
        with _ -> ()

/// The bytes of a file that may still be open: BioFSharp.Mz does not dispose the commands of its look ups, so SQLite keeps the file open until they are finalized.
let private readBytes (path: string) =
    use stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)
    use memory = new MemoryStream()
    stream.CopyTo memory
    memory.ToArray()

let private scalar (dbPath: string) (query: string) =
    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" dbPath)
    cn.Open()
    use cmd = new SQLiteCommand(query, cn)
    Convert.ToInt64 (cmd.ExecuteScalar())

let private execute (dbPath: string) (query: string) =
    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" dbPath)
    cn.Open()
    use cmd = new SQLiteCommand(query, cn)
    cmd.ExecuteNonQuery() |> ignore

let private sameBits (a: (float * float)[]) (b: (float * float)[]) =
    a.Length = b.Length
    && Array.forall2 (fun (x1: float, y1: float) (x2: float, y2: float) -> BitConverter.DoubleToInt64Bits x1 = BitConverter.DoubleToInt64Bits x2 && BitConverter.DoubleToInt64Bits y1 = BitConverter.DoubleToInt64Bits y2) a b

/// The requests of the quantification tools for the ions: both global modifications of a data base with isotopic modifications, otherwise the one of the ion.
let private requests (labelled: bool) (ions: (string * int * int)[]) =
    ions
    |> Array.collect (fun (sequence, globalMod, charge) -> (if labelled then [| 0; 1 |] else [| globalMod |]) |> Array.map (fun gMod -> sequence, gMod, charge))
    |> Array.distinct

/// Looks up every request through the IsotopePatternLookUp on the data base file and returns the counters and whether all patterns are bitwise equal to the computed ones.
let private lookUpAll (logger: NLog.Logger) (dbPath: string) (requests: (string * int * int)[]) =
    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" dbPath)
    cn.Open()
    let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBySequenceAndGMod cn (SearchDB.getSDBParamsByCn cn)
    use patterns = new IsotopePatterns.IsotopePatternLookUp(logger, cn)
    let equal =
        requests
        |> Array.forall (fun (sequence, gMod, charge) ->
            let peptide = lookUp sequence gMod
            sameBits (patterns.Get (peptide, charge)) (IsotopePatterns.compute charge peptide.BioSequence))
    patterns.Read, patterns.Computed, equal

[<Tests>]
let tests =
    testList "IsotopePatterns" [
        testCase "encode and decode keep the bits" <| fun _ ->
            let nanWithPayload = BitConverter.Int64BitsToDouble 0x7FF8000000000123L
            let peaks = [| -0.0, 0.0; 5e-324, -5e-324; nanWithPayload, nan; infinity, -infinity; 1.0000000000000002, 0.1; 1234.5678901234567, 0.001 |]
            let decoded = IsotopePatterns.decode (IsotopePatterns.encode peaks)
            Expect.isSome decoded "the bytes decode"
            Expect.isTrue (sameBits peaks decoded.Value) "all bits are kept"
            Expect.equal (IsotopePatterns.encode [| 1.0, 2.0 |]) [| 0uy; 0uy; 0uy; 0uy; 0uy; 0uy; 0xF0uy; 0x3Fuy; 0uy; 0uy; 0uy; 0uy; 0uy; 0uy; 0uy; 0x40uy |] "float64 little-endian pairs"
            Expect.equal (IsotopePatterns.decode [||]) (Some [||]) "no peaks"
            Expect.isNone (IsotopePatterns.decode (Array.zeroCreate 24)) "a length that is no multiple of 16 is not a pattern"
            Expect.isNone (IsotopePatterns.decode null) "null is not a pattern"

        testCase "the stored patterns equal the computed ones, a second store writes nothing" <| fun _ ->
            withCopy labelledDb (fun db ->
                let logger, target = memoryLogger ()
                let ions = readIons labelledQpsm
                Expect.isGreaterThan ions.Length 0 "the .qpsm holds ions"
                IsotopePatterns.store logger db ions
                let all = requests true ions
                let stored, computed, equal = lookUpAll logger db all
                Expect.equal stored all.Length "every pattern of both global modifications is read from the data base"
                Expect.equal computed 0 "no pattern is computed"
                Expect.isTrue equal "the stored patterns are bitwise equal to the computed ones"
                let keys =
                    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" db)
                    cn.Open()
                    let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBySequenceAndGMod cn (SearchDB.getSDBParamsByCn cn)
                    all |> Array.map (fun (sequence, gMod, charge) -> (lookUp sequence gMod).ModSequenceID, charge) |> Array.distinct
                Expect.equal (scalar db "SELECT COUNT(*) FROM IsotopePattern") (int64 keys.Length) "one row per ModSequence row and charge"
                Expect.equal (scalar db "SELECT COUNT(*) FROM IsotopePatternSettings") 1L "one settings row"
                Expect.isEmpty (target.Logs |> Seq.filter (fun x -> x.StartsWith "WARN") |> Array.ofSeq) "no warning"
                // a second store finds everything and does not touch the file
                let bytes = readBytes db
                let written = File.GetLastWriteTimeUtc db
                Threading.Thread.Sleep 50
                IsotopePatterns.store logger db ions
                Expect.equal (readBytes db) bytes "the data base is unchanged"
                Expect.equal (File.GetLastWriteTimeUtc db) written "the modification time is unchanged"
                Expect.stringContains (String.concat "\n" target.Logs) "no pattern missing" "the log tells that nothing was missing")

        testCase "a data base without the table is not changed by the look up" <| fun _ ->
            withCopy labelledDb (fun db ->
                let logger, target = memoryLogger ()
                let ions = readIons labelledQpsm
                let all = requests true ions
                let bytes = readBytes db
                let stored, computed, equal = lookUpAll logger db all
                Expect.equal stored 0 "nothing is read"
                Expect.equal computed all.Length "everything is computed"
                Expect.isTrue equal "the computed patterns are the same"
                Expect.equal (readBytes db) bytes "the data base is unchanged"
                Expect.stringContains (String.concat "\n" target.Logs) "holds no isotope patterns" "the log tells that the patterns are computed")

        testCase "a peptide that is not in the data base is computed and not stored" <| fun _ ->
            withCopy labelledDb (fun db ->
                let logger, target = memoryLogger ()
                IsotopePatterns.store logger db (Array.append (readIons labelledQpsm) [| "PEPTIDEZZZK", 0, 2 |])
                use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" db)
                cn.Open()
                let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBySequenceAndGMod cn (SearchDB.getSDBParamsByCn cn)
                let unknown = lookUp "PEPTIDEZZZK" 0
                Expect.isLessThan unknown.ModSequenceID 0 "the peptide is not in the data base"
                use patterns = new IsotopePatterns.IsotopePatternLookUp(logger, cn)
                patterns.Get (unknown, 2) |> ignore
                Expect.equal (patterns.Read, patterns.Computed) (0, 1) "the pattern is computed"
                Expect.stringContains (String.concat " " target.Logs) "2 not in the data base" "the log counts the ions that are not in the data base"
                Expect.equal (scalar db "SELECT COUNT(*) FROM IsotopePattern WHERE ModSequenceID < 0") 0L "no row of a missing peptide")

        testCase "patterns of another settings row are ignored" <| fun _ ->
            withCopy labelledDb (fun db ->
                let logger, _ = memoryLogger ()
                let ions = readIons labelledQpsm
                IsotopePatterns.store logger db ions
                execute db "UPDATE IsotopePatternSettings SET Settings = 'another version'"
                let all = requests true ions
                let stored, computed, equal = lookUpAll logger db all
                Expect.equal (stored, computed) (0, all.Length) "the patterns of the other settings are not used"
                Expect.isTrue equal "the computed patterns are the same"
                // the store adds the patterns under its own settings row
                IsotopePatterns.store logger db ions
                Expect.equal (scalar db "SELECT COUNT(*) FROM IsotopePatternSettings") 2L "a settings row per version"
                let stored, computed, _ = lookUpAll logger db all
                Expect.equal (stored, computed) (all.Length, 0) "the patterns of the own settings are used")

        testCase "a data base that cannot be written gives a warning" <| fun _ ->
            withCopy labelledDb (fun db ->
                let logger, target = memoryLogger ()
                File.SetAttributes(db, FileAttributes.ReadOnly)
                let bytes = readBytes db
                IsotopePatterns.store logger db (readIons labelledQpsm)
                Expect.equal (readBytes db) bytes "the data base is unchanged"
                Expect.isTrue (target.Logs |> Seq.exists (fun x -> x.StartsWith "WARN" && x.Contains "could not be stored")) "a warning")

        testCase "parallel stores with overlapping ions store each pattern once" <| fun _ ->
            withCopy labelledDb (fun db ->
                let logger, target = memoryLogger ()
                let ions = readIons labelledQpsm
                // the ions of the stores overlap, as the ions of the files of one PSMStatistics run do
                let parts = [| ions; ions; Array.rev ions; ions |]
                parts |> Array.Parallel.iter (fun part -> IsotopePatterns.store logger db part)
                let all = requests true ions
                let stored, computed, equal = lookUpAll logger db all
                Expect.equal (stored, computed) (all.Length, 0) "all patterns are stored"
                Expect.isTrue equal "the stored patterns are the computed ones"
                Expect.equal (scalar db "SELECT COUNT(*) FROM IsotopePattern") (scalar db "SELECT COUNT(*) FROM (SELECT DISTINCT SettingsID, ModSequenceID, Charge FROM IsotopePattern)") "each key once"
                Expect.equal (scalar db "SELECT COUNT(*) FROM IsotopePatternSettings") 1L "one settings row"
                Expect.isEmpty (target.Logs |> Seq.filter (fun x -> x.StartsWith "WARN") |> Array.ofSeq) "no warning")

        testCase "a data base without isotopic modifications stores the global modification of the ion only" <| fun _ ->
            withCopy unlabelledDb (fun db ->
                let logger, _ = memoryLogger ()
                let ions = readIons unlabelledQpsm
                IsotopePatterns.store logger db (Array.append ions [| "PEPTIDEK", 1, 2 |])
                let all = requests false ions
                let stored, computed, equal = lookUpAll logger db all
                Expect.equal (stored, computed) (all.Length, 0) "the patterns of the ions are stored"
                Expect.isTrue equal "the stored patterns are the computed ones"
                Expect.equal (scalar db "SELECT COUNT(*) FROM IsotopePattern p JOIN ModSequence m ON m.ID = p.ModSequenceID WHERE m.GlobalMod <> 0") 0L "no pattern of GlobalMod 1")

        testCase "AlignmentBasedQuantification isotope peaks equal with stored and computed patterns" <| fun _ ->
            withCopy labelledDb (fun db ->
                let logger, _ = memoryLogger ()
                let ions = readIons labelledQpsm
                IsotopePatterns.store logger db ions
                use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" db)
                cn.Open()
                let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBySequenceAndGMod cn (SearchDB.getSDBParamsByCn cn)
                use patterns = new IsotopePatterns.IsotopePatternLookUp(logger, cn)
                let all = requests true ions
                for (sequence, gMod, charge) in all do
                    let peptide = lookUp sequence gMod
                    let computed = AlignmentBasedQuantification.isotopePeaks (fun peptide charge -> IsotopePatterns.compute charge peptide.BioSequence) charge peptide
                    let stored = AlignmentBasedQuantification.isotopePeaks (fun peptide charge -> patterns.Get (peptide, charge)) charge peptide
                    Expect.equal stored computed (sprintf "%s %i %i" sequence gMod charge)
                Expect.equal (patterns.Read, patterns.Computed) (all.Length, 0) "all patterns are read from the data base")

        testCase "the look up holds no lock between two calls" <| fun _ ->
            withCopy labelledDb (fun db ->
                let logger, _ = memoryLogger ()
                let ions = readIons labelledQpsm
                IsotopePatterns.store logger db ions
                use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" db)
                cn.Open()
                let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBySequenceAndGMod cn (SearchDB.getSDBParamsByCn cn)
                use patterns = new IsotopePatterns.IsotopePatternLookUp(logger, cn)
                let sequence, gMod, charge = ions.[0]
                patterns.Get (lookUp sequence gMod, charge) |> ignore
                // a writer with a timeout of one second gets the lock
                use writer = new SQLiteConnection(sprintf "Data Source=%s;Version=3;Default Timeout=1" db)
                writer.Open()
                use tr = writer.BeginTransaction(Data.IsolationLevel.Serializable)
                use cmd = new SQLiteCommand("DELETE FROM IsotopePatternSettings WHERE ID < 0", writer, tr)
                cmd.ExecuteNonQuery() |> ignore
                tr.Commit())
    ]
