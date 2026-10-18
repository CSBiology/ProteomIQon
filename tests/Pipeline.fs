module Pipeline

open Expecto
open System
open System.IO
open ProteomIQon.Core.InputPaths
open BioFSharp.Mz
open System.Data.SQLite
open FSharpAux.IO

/// Directory of the compiled test assembly (tests/bin/<config>/<tfm>)
let baseDir = AppContext.BaseDirectory

let repoRoot = Path.GetFullPath(Path.Combine(baseDir, "..", "..", "..", ".."))

/// Locates a tool dll in its project build output under src/<projectDir>/bin,
/// preferring the most recently built configuration.
let toolDll (projectDir: string) =
    let assemblyName = sprintf "ProteomIQon.%s.dll" (projectDir.Replace("-", "_"))
    let candidates =
        [ for config in ["Release"; "Debug"] ->
            Path.Combine(repoRoot, "src", projectDir, "bin", config, "net10.0", assemblyName) ]
        |> List.filter File.Exists
    match candidates with
    | [] -> failwithf "%s not found, build it first (dotnet build src/%s)" assemblyName projectDir
    | existing -> existing |> List.maxBy File.GetLastWriteTimeUtc

let runDotNet cmd workingDir =
    let psi =
        Diagnostics.ProcessStartInfo(
            FileName = "dotnet",
            Arguments = cmd,
            WorkingDirectory = workingDir,
            UseShellExecute = false,
            RedirectStandardOutput = true,
            RedirectStandardError = true
        )
    use proc = Diagnostics.Process.Start psi
    // read both streams concurrently to avoid pipe-buffer deadlocks
    let stdOutTask = proc.StandardOutput.ReadToEndAsync()
    let stdErrTask = proc.StandardError.ReadToEndAsync()
    proc.WaitForExit()
    let stdOut = stdOutTask.Result
    let stdErr = stdErrTask.Result
    if proc.ExitCode <> 0 then failwithf "'dotnet %s' failed in %s:\n%s\n%s" cmd workingDir stdOut stdErr

/// Compares two tab separated result files cell by cell. Numeric cells are
/// compared with a relative tolerance, because iterative fits (PEP values,
/// q-values) pick up last-digit differences from the platform's math library
/// and do not reproduce bit for bit between Windows and Linux. Every other
/// cell has to match exactly.
let tsvFilesEqual (relTolerance: float) (referencePath: string) (testPath: string) =
    let tryFloat (s: string) =
        match Double.TryParse(s, Globalization.NumberStyles.Float, Globalization.CultureInfo.InvariantCulture) with
        | true, v -> Some v
        | _ -> None
    let cellsEqual (r: string) (t: string) =
        r = t ||
        match tryFloat r, tryFloat t with
        | Some r, Some t ->
            (Double.IsNaN r && Double.IsNaN t) ||
            abs (r - t) <= relTolerance * max (abs r) (abs t)
        | _ -> false
    let reference = File.ReadAllLines referencePath
    let test = File.ReadAllLines testPath
    reference.Length = test.Length &&
    Array.forall2 (fun (r: string) (t: string) ->
        let rc = r.Split '\t'
        let tc = t.Split '\t'
        rc.Length = tc.Length && Array.forall2 cellsEqual rc tc
    ) reference test

/// The problems of a .qpsm of an estimated threshold, empty when there are none. The trained
/// model, and with it every score, differs between operating systems and with the row order, so
/// the result is checked by properties any correct one has. The bounds cover the spread of runs
/// with shuffled rows, rounded features or another seed: the accepted count stays within 12% of
/// the reference, at most 15% of all accepted PSMs are accepted by only one of the two, every
/// accepted PSM passes the q-value and PEP cutoffs with finite values, the PEPs sum to no more
/// false PSMs than the q-value cutoff allows, and the PEP does not rise as the model score rises.
/// Small results may differ by five PSMs.
let qpsmProblems (qValueThreshold: float) (pepThreshold: float) (referencePath: string) (testPath: string) =
    let read (path: string) =
        let lines = File.ReadAllLines path |> Array.filter (fun line -> line <> "")
        let header = lines.[0].Split '\t'
        let missing = [ "PSMId"; "StringSequence"; "ModelScore"; "QValue"; "PEPValue" ] |> List.filter (fun name -> not (Array.contains name header))
        if not missing.IsEmpty then Error (sprintf "%s lacks the columns %s" path (String.concat ", " missing)) else
        let column name = Array.findIndex ((=) name) header
        let psmId, sequence, score, qValue, pep = column "PSMId", column "StringSequence", column "ModelScore", column "QValue", column "PEPValue"
        let number (s: string) =
            match Double.TryParse(s, Globalization.NumberStyles.Float, Globalization.CultureInfo.InvariantCulture) with
            | true, v -> v
            | _ -> nan
        lines.[1..]
        |> Array.map (fun line ->
            let cells = line.Split '\t'
            (cells.[psmId], cells.[sequence]), number cells.[score], number cells.[qValue], number cells.[pep])
        |> Ok
    match read referencePath, read testPath with
    | Error problem, _ | _, Error problem -> [ problem ]
    | Ok reference, Ok test ->
        let keys rows = rows |> Array.map (fun (key, _, _, _) -> key) |> Set.ofArray
        let union = Set.union (keys reference) (keys test)
        let oneSided = union.Count - (Set.intersect (keys reference) (keys test)).Count
        let summedPep = test |> Array.sumBy (fun (_, _, _, pep) -> pep)
        let byScore = test |> Array.sortBy (fun (_, score, _, _) -> score)
        // PEPs are compared with the relative tolerance of the exact comparisons
        let differs (a: float) (b: float) = abs (a - b) > 1e-9 * max (abs a) (abs b)
        let problems =
            [
                if abs (float test.Length - float reference.Length) > max 5. (0.12 * float reference.Length) then
                    sprintf "%i accepted PSMs, the reference has %i" test.Length reference.Length
                if float oneSided > max 5. (0.15 * float union.Count) then
                    sprintf "%i of %i accepted PSMs are accepted by only the result or the reference" oneSided union.Count
                if summedPep > qValueThreshold * float test.Length then
                    sprintf "The PEPs sum to %g expected false PSMs among %i accepted ones" summedPep test.Length
                for (id, _), score, qValue, pep in test do
                    if not (Array.forall Double.IsFinite [| score; qValue; pep |]) then sprintf "%s has a value that is not finite" id
                    elif qValue > qValueThreshold || pep > pepThreshold then sprintf "%s has q-value %g and PEP %g" id qValue pep
                for (_, lowScore, _, lowPep), ((highId, _), highScore, _, highPep) in Array.pairwise byScore do
                    if highScore > lowScore && highPep > lowPep && differs highPep lowPep then
                        sprintf "The PEP rises from %g at score %g to %g at score %g (%s)" lowPep lowScore highPep highScore highId
                    elif highScore = lowScore && differs highPep lowPep then
                        sprintf "The score %g has the PEPs %g and %g (%s)" highScore lowPep highPep highId
            ]
        if problems.Length > 10 then List.truncate 10 problems @ [ sprintf "and %i more" (problems.Length - 10) ] else problems

/// A copy of a peptide data base in the temp folder, for the tools that write to it.
let tempCopy (source: string) =
    let path = Path.Combine(Path.GetTempPath(), sprintf "pipeline_%s_%s" (Guid.NewGuid().ToString "N") (Path.GetFileName source))
    File.Copy(source, path)
    path

/// Deletes a temporary copy; SQLite keeps the file open until the commands of the look ups of BioFSharp.Mz, which are not disposed, are finalized.
let deleteTemp (path: string) =
    GC.Collect()
    GC.WaitForPendingFinalizers()
    try File.Delete path with _ -> ()

/// The bytes of a file that may still be open in another process.
let fileBytes (path: string) =
    use stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)
    use memory = new MemoryStream()
    stream.CopyTo memory
    memory.ToArray()

/// The distinct (sequence, global modification, charge) ions of a .qpsm file.
let readIons (qpsmPath: string) =
    FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.PSMStatisticsResult>(SchemaMode=FSharpAux.IO.SchemaReader.Csv.Fill).ReadFile(qpsmPath,'\t',false,1)
    |> Seq.map (fun x -> x.StringSequence, x.GlobalMod, x.Charge)
    |> Seq.distinct
    |> Array.ofSeq

/// How many isotope patterns of the (sequence, global modification, charge) ions of a .qpsm file the IsotopePatternLookUp reads from the
/// peptide data base and how many it computes, and the number of requests that are in the data base: both global modifications of every ion if the
/// data base has isotopic modifications, else the global modification of the ion.
let isotopePatternCounts (dbPath: string) (qpsmPath: string) =
    let ions = readIons qpsmPath
    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" dbPath)
    cn.Open()
    let sdbParams = SearchDB.getSDBParamsByCn cn
    let lookUp = SearchDB.getThreadSafePeptideLookUpFromFileBySequenceAndGMod cn sdbParams
    let requests =
        ions
        |> Array.collect (fun (sequence, globalMod, charge) -> (if sdbParams.IsotopicMod.IsEmpty then [| globalMod |] else [| 0; 1 |]) |> Array.map (fun gMod -> sequence, gMod, charge))
        |> Array.distinct
    // a peptide that is not in the data base has the ID -1, its pattern is never stored
    let found =
        requests
        |> Array.map (fun (sequence, gMod, charge) -> lookUp sequence gMod, charge)
        |> Array.filter (fun (peptide, _) -> peptide.ModSequenceID >= 0)
        |> Array.distinctBy (fun (peptide, charge) -> peptide.ModSequenceID, charge)
    use patterns = new ProteomIQon.IsotopePatterns.IsotopePatternLookUp(NLog.LogManager.CreateNullLogger(), cn)
    for (peptide, charge) in found do
        patterns.Get (peptide, charge) |> ignore
    patterns.Read, patterns.Computed, found.Length

/// Whether the data base has the table of the isotope patterns.
let hasIsotopePatternTable (dbPath: string) =
    use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" dbPath)
    cn.Open()
    use cmd = new SQLiteCommand("SELECT COUNT(*) FROM sqlite_master WHERE name = 'IsotopePattern'", cn)
    Convert.ToInt32 (cmd.ExecuteScalar()) > 0

let selectAllModSequence cn =
    let querystring = "SELECT * FROM ModSequence"
    use cmd = new SQLiteCommand(querystring, cn)
    use reader = cmd.ExecuteReader()
    [
        while reader.Read() do 
            yield (reader.GetValue(0), reader.GetValue(1),reader.GetValue(2), reader.GetValue(3), reader.GetValue(4), reader.GetValue(5))
    ]

let selectAllCleavageIndex cn =
    let querystring = "SELECT * FROM CleavageIndex"
    use cmd = new SQLiteCommand(querystring, cn)
    use reader = cmd.ExecuteReader()
    [
        while reader.Read() do 
            yield (reader.GetValue(0), reader.GetValue(1),reader.GetValue(2), reader.GetValue(3), reader.GetValue(4), reader.GetValue(5))
    ]

let selectAllPepSequence cn =
    let querystring = "SELECT * FROM PepSequence"
    use cmd = new SQLiteCommand(querystring, cn)
    use reader = cmd.ExecuteReader()
    [
        while reader.Read() do 
            yield (reader.GetValue(0), reader.GetValue(1))
    ]

let selectAllProtein cn =
    let querystring = "SELECT * FROM Protein"
    use cmd = new SQLiteCommand(querystring, cn)
    use reader = cmd.ExecuteReader()
    [
        while reader.Read() do 
            yield (reader.GetValue(0), reader.GetValue(1))
    ]

/// The parameters with the second derivative peak detection instead of the wavelet. FSharpStats'.Wavelet.identify ignores its parameters and pads the XIC borders
/// with unseeded random numbers, so the quantification with the wavelet differs from run to run. The second derivative has no randomness, but initSpline shuffles the ions
/// of the m/z correction once there are 500 PSMs or more; the 3 PSMs of the test data stay below that.
let secondDerivativeParams (json: string) =
    Text.RegularExpressions.Regex.Replace(
        json,
        "\"XicProcessing\":\s*\{\s*\"Case\":\s*\"Wavelet\".*?\}\s*\]\s*\}",
        "\"XicProcessing\": {\"Case\": \"SecondDerivative\", \"Fields\": [{\"MinSNR\": 0.01, \"PolynomOrder\": 3, \"WindowSize\": {\"Case\": \"Fixed\", \"Fields\": [7]}}]}",
        Text.RegularExpressions.RegexOptions.Singleline)

/// Runs a quantification tool on a peptide data base without the isotope patterns and on a copy that PSMStatistics' store filled with the patterns of the
/// ions of the .qpsm: the logs of the tools show that the first computed every pattern and the second read them, and the .quant files have to be byte identical.
/// The wavelet pads randomly, so with deterministic = false the .quant files only have to hold the same ions.
let isotopePatternEquivalenceCase name (tool: string) (dataDirectory: string) (paramsFile: string) (editParams: string -> string) (deterministic: bool) =
    testCase name <| fun _ ->
        let relToDirectory = getRelativePath baseDir
        let inDirectory = relToDirectory (sprintf "../../../data/%s/in/" dataDirectory)
        let work = Path.Combine(Path.GetTempPath(), sprintf "isotopePatterns_%s" (Guid.NewGuid().ToString "N"))
        Directory.CreateDirectory work |> ignore
        // copies, as the tools lock the mzlite
        let mzlite = Path.Combine(work, "minimal.mzlite")
        let qpsm = Path.Combine(work, "minimal.qpsm")
        let quantParams = Path.Combine(work, paramsFile)
        File.Copy(Path.Combine(inDirectory, "minimal.mzlite"), mzlite)
        File.Copy(Path.Combine(inDirectory, "minimal.qpsm"), qpsm)
        File.WriteAllText(quantParams, editParams (File.ReadAllText (Path.Combine(inDirectory, paramsFile))))
        let dbName = Directory.GetFiles(inDirectory, "*.db") |> Array.exactlyOne
        let dbWithout = Path.Combine(work, "without.db")
        let dbWith = Path.Combine(work, "with.db")
        File.Copy(dbName, dbWithout)
        File.Copy(dbName, dbWith)
        ProteomIQon.IsotopePatterns.store (NLog.LogManager.CreateNullLogger()) dbWith (readIons qpsm)
        let quantExe = toolDll tool
        let run (db: string) (outDirectory: string) =
            runDotNet (sprintf "%s -i %s -ii %s -o %s -p %s -d %s" quantExe mzlite qpsm outDirectory quantParams db) baseDir
            let log = File.ReadAllText (Path.Combine(outDirectory, "minimal_log.txt"))
            let m = Text.RegularExpressions.Regex.Match(log, @"Isotope patterns: (\d+) read from the peptide data base, (\d+) computed")
            Expect.isTrue m.Success "the log has the counts of the isotope patterns"
            File.ReadAllBytes (Path.Combine(outDirectory, "minimal.quant")), int m.Groups.[1].Value, int m.Groups.[2].Value
        let quantWithout, storedWithout, computedWithout = run dbWithout (Path.Combine(work, "without"))
        let quantWith, storedWith, computedWith = run dbWith (Path.Combine(work, "with"))
        let isFilled = hasIsotopePatternTable dbWith
        GC.Collect()
        GC.WaitForPendingFinalizers()
        try Directory.Delete(work, true) with _ -> ()
        Expect.isTrue isFilled "the store filled the data base"
        Expect.equal storedWithout 0 "without the table nothing is read"
        Expect.isGreaterThan computedWithout 0 "without the table the patterns are computed"
        Expect.isGreaterThan storedWith 0 "with the table the patterns are read"
        Expect.equal computedWith 0 "with the table no pattern is computed"
        Expect.equal storedWith computedWithout "the same patterns are requested"
        if deterministic then
            Expect.isTrue (quantWith = quantWithout) "the .quant files are byte identical"
        else
            // the wavelet pads randomly: the same ions in the same order
            let ionsOf (quant: byte[]) =
                Text.Encoding.UTF8.GetString(quant).Split('\n')
                |> Array.map (fun line -> String.Join("\t", line.Split('\t') |> Array.truncate 3))
            Expect.equal (ionsOf quantWith) (ionsOf quantWithout) "the .quant files have the same ions"

/// The spectra of an mzML file, decoded without MzIO: ID, MS level, scan start time and its unit, m/z of the first selected ion
/// (nan without one), and the m/z, intensity and ion mobility arrays (32 or 64 bit floats, zlib compressed or not).
let readMzMLSpectra (path: string) =
    let ns = Xml.Linq.XNamespace.Get "http://psi.hupo.org/ms/mzml"
    let attribute (name: string) (element: Xml.Linq.XElement) = element.Attribute(Xml.Linq.XName.Get name).Value
    let tryCv (accession: string) (element: Xml.Linq.XElement) =
        element.Elements(ns + "cvParam") |> Seq.tryFind (fun p -> attribute "accession" p = accession)
    let decode (array: Xml.Linq.XElement) =
        let bytes = Convert.FromBase64String (array.Element(ns + "binary").Value)
        let raw =
            if (tryCv "MS:1000574" array).IsSome then
                use input = new MemoryStream(bytes)
                use zlib = new Compression.ZLibStream(input, Compression.CompressionMode.Decompress)
                use output = new MemoryStream()
                zlib.CopyTo output
                output.ToArray()
            else bytes
        if (tryCv "MS:1000521" array).IsSome then Array.init (raw.Length / 4) (fun i -> float (BitConverter.ToSingle(raw, 4 * i)))
        else Array.init (raw.Length / 8) (fun i -> BitConverter.ToDouble(raw, 8 * i))
    Xml.Linq.XDocument.Load(path).Descendants(ns + "spectrum")
    |> Seq.map (fun spectrum ->
        let arrays = spectrum.Descendants(ns + "binaryDataArray") |> Array.ofSeq
        let arrayOf accession = arrays |> Array.find (fun a -> (tryCv accession a).IsSome) |> decode
        let scanTime = spectrum.Descendants(ns + "scan") |> Seq.pick (tryCv "MS:1000016")
        let selectedIonMz =
            spectrum.Descendants(ns + "selectedIon") |> Seq.tryHead |> Option.bind (tryCv "MS:1000744")
            |> Option.map (attribute "value" >> float) |> Option.defaultValue nan
        {|
            ID            = attribute "id" spectrum
            MsLevel       = tryCv "MS:1000511" spectrum |> Option.get |> attribute "value" |> int
            ScanTime      = float (attribute "value" scanTime)
            ScanTimeUnit  = attribute "unitAccession" scanTime
            SelectedIonMz = selectedIonMz
            Mz            = arrayOf "MS:1000514"
            Intensity     = arrayOf "MS:1000515"
            IonMobility   = arrayOf "MS:1003006"
        |})
    |> Array.ofSeq

/// The spectra of an mzlite file read with MzIO, by ID: MS level, scan time and its unit, precursor m/z and peaks.
let readMzliteSpectra (path: string) =
    use reader = new MzIO.MzSQL.MzSQL(path)
    reader.Connection.Open()
    reader.ReadMassSpectra "sample=0"
    |> Seq.map (fun spectrum ->
        let scanTimeUnit =
            spectrum.Scans.GetProperties false
            |> Seq.pick (fun scan ->
                match scan.Value with
                | :? MzIO.Model.Scan as s -> s.TryGetValue MzIO.MetaData.PSIMSExtension.PSIMS_Scan.ScanStartTime
                | _ -> None)
            |> fun p -> MzIO.Model.CvParam.tryGetCvUnitAccession (p :?> MzIO.Model.CvParam.IParamBase<IConvertible>)
        spectrum.ID,
        (MzIO.Processing.MassSpectrum.getMsLevel spectrum, MzIO.Processing.MassSpectrum.getScanTime spectrum, scanTimeUnit,
         MzIO.Processing.MassSpectrum.getPrecursorMZ spectrum, reader.ReadSpectrumPeaks(spectrum.ID).Peaks |> Array.ofSeq))
    |> Map.ofSeq

/// Deletes a directory whose SQLite files may still be open: MzIO does not dispose its commands, so a file stays open until they are
/// finalized.
let deleteDirectory (path: string) =
    let rec attempt n =
        GC.Collect()
        GC.WaitForPendingFinalizers()
        GC.Collect()
        try Directory.Delete(path, true)
        with :? IOException | :? UnauthorizedAccessException when n > 0 ->
            Threading.Thread.Sleep 500
            attempt (n - 1)
    attempt 10

[<Tests>]
let pipelineTests =
    testList "Pipeline" [
        testCase "PeptideDB" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let fastaPath = relToDirectory "../../../data/PeptideDB/in/example.fasta"
            let peptideDBParams = relToDirectory "../../../data/PeptideDB/in/peptideDBParams.json"
            let outDirectory = relToDirectory "../../../data/PeptideDB/out/"
            let pepDBExe = toolDll "PeptideDB"
            // run tool
            runDotNet (sprintf "%s -i %s -o %s -p %s" pepDBExe fastaPath outDirectory peptideDBParams) baseDir
            let referenceDB =
                let dbPath = relToDirectory "../../../data/PeptideDB/out/MinimalReference.db"
                use cn = SearchDB.getDBConnection dbPath
                let res = selectAllCleavageIndex cn, selectAllModSequence cn, selectAllPepSequence cn, selectAllProtein cn
                res
            let testDB =
                let dbPath = relToDirectory "../../../data/PeptideDB/out/Minimal.db"
                use cn = SearchDB.getDBConnection dbPath
                let res = selectAllCleavageIndex cn, selectAllModSequence cn, selectAllPepSequence cn, selectAllProtein cn
                res
            let compare = referenceDB = testDB
            // cleanup
            File.Delete (relToDirectory "../../../data/PeptideDB/out/Minimal.db")
            File.Delete (relToDirectory "../../../data/PeptideDB/out/PeptideDB_log.txt")
            File.Delete (relToDirectory "../../../data/PeptideDB/out/PeptideDB_Minimal_log.txt")
            Expect.isTrue compare "Peptide databases are different"
    
        testCase "PeptideSpectrumMatching" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let db = relToDirectory "../../../data/PeptideSpectrumMatching/in/Minimal.db"
            let mzlite = relToDirectory "../../../data/PeptideSpectrumMatching/in/minimal.mzlite"
            let psmParams = relToDirectory "../../../data/PeptideSpectrumMatching/in/defaultParams.json"
            let outDirectory = relToDirectory "../../../data/PeptideSpectrumMatching/out/"
            let psmExe = toolDll "PeptideSpectrumMatching"
            // run tool
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmExe mzlite outDirectory psmParams db) baseDir
            let compare =
                tsvFilesEqual 1e-6
                    (relToDirectory "../../../data/PeptideSpectrumMatching/out/minimalReference.psm")
                    (relToDirectory "../../../data/PeptideSpectrumMatching/out/minimal.psm")
            // cleanup
            File.Delete (relToDirectory "../../../data/PeptideSpectrumMatching/out/minimal.psm")
            File.Delete (relToDirectory "../../../data/PeptideSpectrumMatching/out/minimal_log.txt")
            File.Delete (relToDirectory "../../../data/PeptideSpectrumMatching/out/PeptideSpectrumMatching_log.txt")
            Expect.isTrue compare "PSMs are different"

        testCase "PeptideSpectrumMatchingTIMs" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let db = relToDirectory "../../../data/PeptideSpectrumMatching/in/Minimal.db"
            let mzlite = relToDirectory "../../../data/PeptideSpectrumMatching/in/minimal.mzlite"
            let psmParams = relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/in/defaultParams.json"
            let outDirectory = relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/out/"
            let psmExe = toolDll "PeptideSpectrumMatchingTIMs"
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmExe mzlite outDirectory psmParams db) baseDir
            let compare =
                tsvFilesEqual 1e-6
                    (relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/out/minimalReference.psm")
                    (relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/out/minimal.psm")
            File.Delete (relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/out/minimal.psm")
            File.Delete (relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/out/minimal_log.txt")
            File.Delete (relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/out/PeptideSpectrumMatchingTIMs_log.txt")
            Expect.isTrue compare "PSMs are different"

        testCase "PeptideSpectrumMatchingTIMsIonMobility" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let db = relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/in/MinimalTIMs.db"
            let mzlite = relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/in/minimalTIMs.mzlite"
            let psmParams = relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/in/defaultParams.json"
            // its own folder, because the test without ion mobility runs at the same time and writes
            // and deletes the tool log in out/
            let outDirectory = relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/out/ionMobility/"
            let psmExe = toolDll "PeptideSpectrumMatchingTIMs"
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmExe mzlite outDirectory psmParams db) baseDir
            let compare =
                tsvFilesEqual 1e-6
                    (relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/out/minimalTIMsReference.psm")
                    (relToDirectory "../../../data/PeptideSpectrumMatchingTIMs/out/ionMobility/minimalTIMs.psm")
            Directory.Delete(outDirectory, true)
            Expect.isTrue compare "PSMs are different"

        testCase "PSMStatistics" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            // the tool stores the isotope patterns in the peptide data base, so it runs on copies; the tracked data bases must not change
            let originalEstimate = relToDirectory "../../../data/PSMStatistics/in/MinimalEstimate.db"
            let originalFixed = relToDirectory "../../../data/PSMStatistics/in/MinimalFixed.db"
            let originalBytes = fileBytes originalEstimate, fileBytes originalFixed
            let dbEstimate = tempCopy originalEstimate
            let dbFixed = tempCopy originalFixed
            let psmEstimate = relToDirectory "../../../data/PSMStatistics/in/minimalEstimate.psm"
            let psmFixed = relToDirectory "../../../data/PSMStatistics/in/minimalFixed.psm"
            let psmStatsParamsEstimate = relToDirectory "../../../data/PSMStatistics/in/pSMStatisticsParamsEstimate.json"
            let psmStatsParamsFixed = relToDirectory "../../../data/PSMStatistics/in/pSMStatisticsParamsFixed.json"
            let outDirectoryEstimate = relToDirectory "../../../data/PSMStatistics/out/estimateOut"
            let outDirectoryFixed = relToDirectory "../../../data/PSMStatistics/out/fixedOut"
            let psmStatsExe = toolDll "PSMStatistics"
            // run tool
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmStatsExe psmEstimate outDirectoryEstimate psmStatsParamsEstimate dbEstimate) baseDir
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmStatsExe psmFixed outDirectoryFixed psmStatsParamsFixed dbFixed) baseDir
            let estimateCounts = isotopePatternCounts dbEstimate (relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalEstimate.qpsm")
            let fixedCounts = isotopePatternCounts dbFixed (relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalFixed.qpsm")
            let unchanged = (fileBytes originalEstimate, fileBytes originalFixed) = originalBytes
            // the cutoffs of pSMStatisticsParamsEstimate.json
            let estimateProblems =
                qpsmProblems 0.01 0.05
                    (relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalReference.qpsm")
                    (relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalEstimate.qpsm")
            let fixedEqual =
                tsvFilesEqual 1e-9
                    (relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalReference.qpsm")
                    (relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalFixed.qpsm")
            // cleanup
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalFixed.qpsm")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalFixed_log.txt")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/fixedOut/PSMStatistics_log.txt")
            Directory.Delete (relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalFixed_plots")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalEstimate.qpsm")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalEstimate_log.txt")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/estimateOut/PSMStatistics_log.txt")
            Directory.Delete (relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalEstimate_plots")
            deleteTemp dbEstimate
            deleteTemp dbFixed
            Expect.isTrue (List.isEmpty estimateProblems) (String.concat "; " estimateProblems)
            Expect.isTrue fixedEqual "QPSMs are different"
            Expect.isTrue unchanged "the data bases of the test data are not changed"
            for name, (stored, computed, requests) in [ "estimate", estimateCounts; "fixed", fixedCounts ] do
                Expect.isGreaterThan requests 0 (sprintf "%s: the .qpsm holds ions" name)
                Expect.equal (stored, computed) (requests, 0) (sprintf "%s: the data base holds the patterns of the ions of the .qpsm" name)

        testCase "PSMStatisticsTIMs" <| fun _ ->
            // a .psm of PeptideSpectrumMatchingTIMs, which carries the hyperscore and the expectation value
            let relToDirectory = getRelativePath baseDir
            let db = tempCopy (relToDirectory "../../../data/PSMStatistics/in/MinimalTIMs.db")
            let psm = relToDirectory "../../../data/PSMStatistics/in/minimalTIMs.psm"
            let psmStatsParams = relToDirectory "../../../data/PSMStatistics/in/pSMStatisticsParamsEstimate.json"
            let outDirectory = relToDirectory "../../../data/PSMStatistics/out/timsOut"
            let psmStatsExe = toolDll "PSMStatistics"
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmStatsExe psm outDirectory psmStatsParams db) baseDir
            let log = File.ReadAllText (relToDirectory "../../../data/PSMStatistics/out/timsOut/minimalTIMs_log.txt")
            let stored, computed, requests = isotopePatternCounts db (relToDirectory "../../../data/PSMStatistics/out/timsOut/minimalTIMs.qpsm")
            // the cutoffs of pSMStatisticsParamsEstimate.json
            let problems =
                qpsmProblems 0.01 0.05
                    (relToDirectory "../../../data/PSMStatistics/out/timsOut/minimalReference.qpsm")
                    (relToDirectory "../../../data/PSMStatistics/out/timsOut/minimalTIMs.qpsm")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/timsOut/minimalTIMs.qpsm")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/timsOut/minimalTIMs_log.txt")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/timsOut/PSMStatistics_log.txt")
            Directory.Delete (relToDirectory "../../../data/PSMStatistics/out/timsOut/minimalTIMs_plots")
            deleteTemp db
            Expect.isGreaterThan requests 0 "the .qpsm holds ions"
            Expect.equal (stored, computed) (requests, 0) "the data base holds the patterns of both global modifications of the ions"
            Expect.stringContains log "Hyperscore and expectation value present: true" "the hyperscore and the expectation value are used"
            Expect.isTrue (List.isEmpty problems) (String.concat "; " problems)

        testCase "PSMStatisticsCrossValidation" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let db = tempCopy (relToDirectory "../../../data/PSMStatistics/in/MinimalTIMs.db")
            let psm = relToDirectory "../../../data/PSMStatistics/in/minimalTIMs.psm"
            let psmStatsParams = relToDirectory "../../../data/PSMStatistics/in/pSMStatisticsParamsCrossValidation.json"
            let outDirectory = relToDirectory "../../../data/PSMStatistics/out/timsCvOut"
            let psmStatsExe = toolDll "PSMStatistics"
            if Directory.Exists outDirectory then Directory.Delete(outDirectory, true)
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmStatsExe psm outDirectory psmStatsParams db) baseDir
            let log = File.ReadAllText (Path.Combine(outDirectory, "minimalTIMs_log.txt"))
            let qpsmLines = File.ReadAllLines (Path.Combine(outDirectory, "minimalTIMs.qpsm")) |> Array.filter (fun line -> line <> "")
            Directory.Delete(outDirectory, true)
            deleteTemp db
            Expect.stringContains log "Cross-validated scores" "the scans are scored by cross-validation"
            Expect.isGreaterThan qpsmLines.Length 1 "the qpsm holds PSMs"

        testCase "PSMBasedQuantification" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let db = relToDirectory "../../../data/PSMBasedQuantification/in/Minimal.db"
            let qpsm = relToDirectory "../../../data/PSMBasedQuantification/in/minimal.qpsm"
            let mzlite = relToDirectory "../../../data/PSMBasedQuantification/in/minimal.mzlite"
            let quantParams = relToDirectory "../../../data/PSMBasedQuantification/in/QuantificationParams.json"
            let outDirectory = relToDirectory "../../../data/PSMBasedQuantification/out"
            let quantExe = toolDll "PSMBasedQuantification"
            // cleanup
            try
                File.Delete (relToDirectory "../../../data/PSMBasedQuantification/out/minimal.quant")
                File.Delete (relToDirectory "../../../data/PSMBasedQuantification/out/minimal_log.txt")
                File.Delete (relToDirectory "../../../data/PSMBasedQuantification/out/PSMBasedQuantification_log.txt")
                Directory.Delete (relToDirectory "../../../data/PSMBasedQuantification/out/minimal_plots")
            with
            | _ -> ()
            // run tool
            runDotNet (sprintf "%s -i %s -ii %s -o %s -p %s -d %s -dc" quantExe mzlite qpsm outDirectory quantParams db) baseDir
            let referenceQuant =
                let quantPath = relToDirectory "../../../data/PSMBasedQuantification/out/minimalReference.quant"
                FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.QuantificationResult>().ReadFile(quantPath,'\t',false,1)
                |> Array.ofSeq
            let testQuant =
                let quantPath = relToDirectory "../../../data/PSMBasedQuantification/out/minimal.quant"
                FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.QuantificationResult>().ReadFile(quantPath,'\t',false,1)
                |> Array.ofSeq
            let nanCasesForUnlabledData lowerLimit upperLimit actual reference = 
                if Double.IsNaN actual && Double.IsNaN reference  then 
                    true
                elif Double.IsNaN actual || Double.IsNaN reference then
                    false
                else 
                    actual >= reference * lowerLimit &&
                    actual <= reference * upperLimit
            let compare,fields = 
                let unequalFields =
                    Array.map2 (fun (reference: ProteomIQon.Dto.QuantificationResult) (test: ProteomIQon.Dto.QuantificationResult) ->
                        [|
                            "StringSequence",
                                test.StringSequence = reference.StringSequence
                            "GlobalMod",
                                test.GlobalMod = reference.GlobalMod
                            "Charge",
                                test.Charge = reference.Charge
                            "PepSequenceID",
                                test.PepSequenceID = reference.PepSequenceID
                            "ModSequenceID",
                                test.ModSequenceID = reference.ModSequenceID
                            "PrecursorMZ",
                                nanCasesForUnlabledData 0.99 1.01 test.PrecursorMZ reference.PrecursorMZ
                            "MeasuredMass",
                                nanCasesForUnlabledData 0.99 1.01 test.MeasuredMass reference.MeasuredMass
                            "TheoMass"
                                ,test.TheoMass = reference.TheoMass
                            "AbsDeltaMass",
                                nanCasesForUnlabledData 0.7 1.3 test.AbsDeltaMass reference.AbsDeltaMass
                            "ProteinNames",
                                test.ProteinNames = reference.ProteinNames
                            "QuantMz_Light",
                                nanCasesForUnlabledData 0.99 1.01 test.QuantMz_Light reference.QuantMz_Light
                            "Quant_Light",
                                nanCasesForUnlabledData 0.99 1.01 test.Quant_Light reference.Quant_Light
                            "MeasuredApex_Light",
                                nanCasesForUnlabledData 0.99 1.01 test.MeasuredApex_Light reference.MeasuredApex_Light
                            //"Seo_Light",
                            //    test.Seo_Light >= reference.Seo_Light * 0.99 &&
                            //    test.Seo_Light <= reference.Seo_Light * 1.01
                            "Difference_SearchRT_FittedRT_Light",
                                nanCasesForUnlabledData 0.9 1.1 (abs test.Difference_SearchRT_FittedRT_Light) (abs reference.Difference_SearchRT_FittedRT_Light)
                            "KLDiv_CorrectedObserved_Theoretical_Light",
                                nanCasesForUnlabledData 0.9 1.1 test.KLDiv_CorrectedObserved_Theoretical_Light reference.KLDiv_CorrectedObserved_Theoretical_Light
                            "KLDiv_Observed_Theoretical_Light",
                                nanCasesForUnlabledData 0.99 1.01 test.KLDiv_Observed_Theoretical_Light reference.KLDiv_Observed_Theoretical_Light
                            "QuantMz_Heavy",
                                nanCasesForUnlabledData 0.99 1.01 test.QuantMz_Heavy reference.QuantMz_Heavy
                            "Quant_Heavy",
                                nanCasesForUnlabledData 0.99 1.01 test.Quant_Heavy reference.Quant_Heavy
                            "MeasuredApex_Heavy",
                                nanCasesForUnlabledData 0.99 1.01 test.MeasuredApex_Heavy reference.MeasuredApex_Heavy
                            //"Seo_Heavy",
                            //    test.Seo_Heavy >= reference.Seo_Heavy * 0.99 &&
                            //    test.Seo_Heavy <= reference.Seo_Heavy * 1.01
                            "Difference_SearchRT_FittedRT_Heavy",
                                nanCasesForUnlabledData 0.9 1.1 (abs test.Difference_SearchRT_FittedRT_Heavy) (abs reference.Difference_SearchRT_FittedRT_Heavy)
                            "KLDiv_CorrectedObserved_Theoretical_Heavy",
                                nanCasesForUnlabledData 0.9 1.1 test.KLDiv_CorrectedObserved_Theoretical_Heavy reference.KLDiv_CorrectedObserved_Theoretical_Heavy
                            "KLDiv_Observed_Theoretical_Heavy",
                                nanCasesForUnlabledData 0.9 1.1 test.KLDiv_Observed_Theoretical_Heavy reference.KLDiv_Observed_Theoretical_Heavy
                            "Correlation_Light_Heavy",
                                nanCasesForUnlabledData 0.99 1.01 test.Correlation_Light_Heavy reference.Correlation_Light_Heavy
                            "QuantificationSource",
                                test.QuantificationSource = reference.QuantificationSource
                        |]
                    ) referenceQuant testQuant
                    |> Array.collect (
                        Array.filter (fun (field,equal) -> equal = false)
                    )
                    |> Array.distinct
                    |> Array.map fst
                if unequalFields |> Array.isEmpty then true , ""
                else false, unequalFields |> String.concat ";"
            Expect.isTrue compare (sprintf "Quants are different in the following fields: %s" fields)
        testCase "PSMStatisticsStoreIsotopePatternsFalse" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let original = relToDirectory "../../../data/PSMStatistics/in/MinimalFixed.db"
            let db = tempCopy original
            let psm = relToDirectory "../../../data/PSMStatistics/in/minimalFixed.psm"
            let outDirectory = Path.Combine(Path.GetTempPath(), sprintf "psmStatistics_%s" (Guid.NewGuid().ToString "N"))
            Directory.CreateDirectory outDirectory |> ignore
            // the fixed cutoffs of the test with the parameter off
            let psmStatsParams = Path.Combine(outDirectory, "params.json")
            File.WriteAllText(psmStatsParams, (File.ReadAllText (relToDirectory "../../../data/PSMStatistics/in/pSMStatisticsParamsFixed.json")).Replace("\"StoreIsotopePatterns\":true", "\"StoreIsotopePatterns\":false"))
            let before = fileBytes db
            let psmStatsExe = toolDll "PSMStatistics"
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmStatsExe psm outDirectory psmStatsParams db) baseDir
            let qpsmLines = File.ReadAllLines (Path.Combine(outDirectory, "minimalFixed.qpsm")) |> Array.filter (fun line -> line <> "")
            let unchanged = fileBytes db = before
            let hasTable = hasIsotopePatternTable db
            deleteTemp db
            try Directory.Delete(outDirectory, true) with _ -> ()
            Expect.isGreaterThan qpsmLines.Length 1 "the .qpsm holds PSMs"
            Expect.isTrue unchanged "the data base is not written to"
            Expect.isFalse hasTable "the data base has no table of isotope patterns"

        testCase "PSMStatisticsIsotopePatternsOfSeveralFiles" <| fun _ ->
            // two files with the same ions: the union is stored once
            let relToDirectory = getRelativePath baseDir
            let db = tempCopy (relToDirectory "../../../data/PSMStatistics/in/MinimalFixed.db")
            let inDirectory = Path.Combine(Path.GetTempPath(), sprintf "psmStatistics_%s" (Guid.NewGuid().ToString "N"))
            let outDirectory = Path.Combine(inDirectory, "out")
            Directory.CreateDirectory inDirectory |> ignore
            File.Copy(relToDirectory "../../../data/PSMStatistics/in/minimalFixed.psm", Path.Combine(inDirectory, "a.psm"))
            File.Copy(relToDirectory "../../../data/PSMStatistics/in/minimalFixed.psm", Path.Combine(inDirectory, "b.psm"))
            let psmStatsExe = toolDll "PSMStatistics"
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s -c 2" psmStatsExe inDirectory outDirectory (relToDirectory "../../../data/PSMStatistics/in/pSMStatisticsParamsFixed.json") db) baseDir
            let stored, computed, requests = isotopePatternCounts db (Path.Combine(outDirectory, "a.qpsm"))
            let rows =
                use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3" db)
                cn.Open()
                use cmd = new SQLiteCommand("SELECT COUNT(*) FROM IsotopePattern", cn)
                Convert.ToInt32 (cmd.ExecuteScalar())
            deleteTemp db
            try Directory.Delete(inDirectory, true) with _ -> ()
            Expect.isGreaterThan requests 0 "the .qpsm holds ions"
            Expect.equal (stored, computed) (requests, 0) "the data base holds the patterns of the ions"
            Expect.equal rows requests "each pattern is stored once"

        // sequenced, as the tools lock the data bases they read and the runs compare the CPU bound tools
        isotopePatternEquivalenceCase "PSMBasedQuantificationIsotopePatterns" "PSMBasedQuantification" "PSMBasedQuantification" "QuantificationParams.json" secondDerivativeParams true |> testSequenced
        isotopePatternEquivalenceCase "PSMBasedQuantificationIsotopePatternsLabelFree" "PSMBasedQuantification" "PSMBasedQuantification" "QuantificationParams.json" (fun json -> (secondDerivativeParams json).Replace("N15Labeling", "Unlabeled")) true |> testSequenced
        isotopePatternEquivalenceCase "PSMBasedQuantificationWaveletIsotopePatterns" "PSMBasedQuantification" "PSMBasedQuantification" "QuantificationParams.json" id false |> testSequenced
        isotopePatternEquivalenceCase "PSMBasedQuantificationTIMsIsotopePatterns" "PSMBasedQuantificationTIMs" "PSMBasedQuantificationTIMs" "QuantificationTIMsParams.json" id true |> testSequenced
        isotopePatternEquivalenceCase "PSMBasedQuantificationTIMsLabelledIsotopePatterns" "PSMBasedQuantificationTIMs" "PSMBasedQuantificationTIMs" "QuantificationTIMsParams.json" (fun json -> json.Replace("Unlabeled", "N15Labeling")) true |> testSequenced
        isotopePatternEquivalenceCase "PSMBasedQuantificationTIMsWaveletIsotopePatterns" "PSMBasedQuantificationTIMs" "PSMBasedQuantificationTIMs" "QuantificationTIMsWaveletParams.json" id false |> testSequenced

        testCase "MsFraggerToPSM" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let tsv = relToDirectory "../../../data/MsFraggerToPSM/in/minimal.tsv"
            let db = relToDirectory "../../../data/PSMBasedQuantificationTIMs/in/minimal.db"
            let outDirectory = relToDirectory "../../../data/MsFraggerToPSM/out/test/"
            let msFraggerExe = toolDll "MsFraggerToPSM"
            if Directory.Exists outDirectory then Directory.Delete(outDirectory, true)
            Directory.CreateDirectory outDirectory |> ignore
            // a copy, because the PSMBasedQuantificationTIMs test locks the mzlite at the same time
            let mzlite = Path.Combine(outDirectory, "minimal.mzlite")
            File.Copy(relToDirectory "../../../data/PSMBasedQuantificationTIMs/in/minimal.mzlite", mzlite)
            runDotNet (sprintf "%s -i %s -ii %s -o %s -d %s" msFraggerExe tsv mzlite outDirectory db) baseDir
            let testPath = relToDirectory "../../../data/MsFraggerToPSM/out/test/minimal.qpsm"
            // the output has the column layout of PSMStatistics
            let header = (File.ReadLines testPath |> Seq.head).Split '\t'
            let expectedHeader =
                Reflection.FSharpType.GetRecordFields(typeof<ProteomIQon.Dto.PSMStatisticsResult>)
                |> Array.map (fun field -> field.Name)
            let compare =
                tsvFilesEqual 1e-9
                    (relToDirectory "../../../data/MsFraggerToPSM/out/minimalReference.qpsm")
                    testPath
            Directory.Delete(outDirectory, true)
            Expect.equal header expectedHeader "the header is the one of PSMStatisticsResult"
            Expect.isTrue compare "QPSMs are different"

        let quantTIMsCase name paramsFile outSubDir referenceFile =
            testCase name <| fun _ ->
                let relToDirectory = getRelativePath baseDir
                let db = relToDirectory "../../../data/PSMBasedQuantificationTIMs/in/minimal.db"
                let qpsm = relToDirectory "../../../data/PSMBasedQuantificationTIMs/in/minimal.qpsm"
                let mzlite = relToDirectory "../../../data/PSMBasedQuantificationTIMs/in/minimal.mzlite"
                let quantParamsTIMs = relToDirectory "../../../data/PSMBasedQuantificationTIMs/in/" + paramsFile
                let outDirectory = relToDirectory ("../../../data/PSMBasedQuantificationTIMs/out" + outSubDir)
                let quantTIMsExe = toolDll "PSMBasedQuantificationTIMs"
                // cleanup
                try
                    File.Delete (Path.Combine(outDirectory, "minimal.quant"))
                    File.Delete (Path.Combine(outDirectory, "minimal_log.txt"))
                    File.Delete (Path.Combine(outDirectory, "PSMBasedQuantification_log.txt"))
                    Directory.Delete (Path.Combine(outDirectory, "minimal_plots"))
                with
                | _ -> ()
                // run tool
                runDotNet (sprintf "%s -i %s -ii %s -o %s -p %s -d %s -dc" quantTIMsExe mzlite qpsm outDirectory quantParamsTIMs db) baseDir
                let referenceQuantTIMs =
                    let quantPath = relToDirectory "../../../data/PSMBasedQuantificationTIMs/out/" + referenceFile
                    FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.QuantificationResult>().ReadFile(quantPath,'\t',false,1)
                    |> Array.ofSeq
                let testQuantTIMs =
                    let quantPath = Path.Combine(outDirectory, "minimal.quant")
                    FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.QuantificationResult>().ReadFile(quantPath,'\t',false,1)
                    |> Array.ofSeq
                let nanCasesForUnlabledData lowerLimit upperLimit actual reference = 
                    if Double.IsNaN actual && Double.IsNaN reference  then 
                        true
                    elif Double.IsNaN actual || Double.IsNaN reference then
                        false
                    else 
                        actual >= reference * lowerLimit &&
                        actual <= reference * upperLimit
                let compare,fields = 
                    let unequalFields =
                        Array.map2 (fun (reference: ProteomIQon.Dto.QuantificationResult) (test: ProteomIQon.Dto.QuantificationResult) ->
                            [|
                                "StringSequence",
                                    test.StringSequence = reference.StringSequence
                                "GlobalMod",
                                    test.GlobalMod = reference.GlobalMod
                                "Charge",
                                    test.Charge = reference.Charge
                                "PepSequenceID",
                                    test.PepSequenceID = reference.PepSequenceID
                                "ModSequenceID",
                                    test.ModSequenceID = reference.ModSequenceID
                                "PrecursorMZ",
                                    nanCasesForUnlabledData 0.99 1.01 test.PrecursorMZ reference.PrecursorMZ
                                "MeasuredMass",
                                    nanCasesForUnlabledData 0.99 1.01 test.MeasuredMass reference.MeasuredMass
                                "TheoMass"
                                    ,test.TheoMass = reference.TheoMass
                                "AbsDeltaMass",
                                    nanCasesForUnlabledData 0.7 1.3 test.AbsDeltaMass reference.AbsDeltaMass
                                "MeanPercolatorScore",
                                    nanCasesForUnlabledData 0.99 1.01 test.MeanPercolatorScore reference.MeanPercolatorScore
                                "QValue",
                                    nanCasesForUnlabledData 0.99 1.01 test.QValue reference.QValue
                                "PEPValue",
                                    nanCasesForUnlabledData 0.99 1.01 test.PEPValue reference.PEPValue
                                "ProteinNames",
                                    test.ProteinNames = reference.ProteinNames
                                "QuantMz_Light",
                                    nanCasesForUnlabledData 0.99 1.01 test.QuantMz_Light reference.QuantMz_Light
                                "Quant_Light",
                                    nanCasesForUnlabledData 0.99 1.01 test.Quant_Light reference.Quant_Light
                                "MeasuredApex_Light",
                                    nanCasesForUnlabledData 0.99 1.01 test.MeasuredApex_Light reference.MeasuredApex_Light
                                //"Seo_Light",
                                //    test.Seo_Light >= reference.Seo_Light * 0.99 &&
                                //    test.Seo_Light <= reference.Seo_Light * 1.01
                                "Difference_SearchRT_FittedRT_Light",
                                    nanCasesForUnlabledData 0.9 1.1 (abs test.Difference_SearchRT_FittedRT_Light) (abs reference.Difference_SearchRT_FittedRT_Light)
                                "KLDiv_CorrectedObserved_Theoretical_Light",
                                    nanCasesForUnlabledData 0.9 1.1 test.KLDiv_CorrectedObserved_Theoretical_Light reference.KLDiv_CorrectedObserved_Theoretical_Light
                                "KLDiv_Observed_Theoretical_Light",
                                    nanCasesForUnlabledData 0.99 1.01 test.KLDiv_Observed_Theoretical_Light reference.KLDiv_Observed_Theoretical_Light
                                "QuantMz_Heavy",
                                    nanCasesForUnlabledData 0.99 1.01 test.QuantMz_Heavy reference.QuantMz_Heavy
                                "Quant_Heavy",
                                    nanCasesForUnlabledData 0.99 1.01 test.Quant_Heavy reference.Quant_Heavy
                                "MeasuredApex_Heavy",
                                    nanCasesForUnlabledData 0.99 1.01 test.MeasuredApex_Heavy reference.MeasuredApex_Heavy
                                //"Seo_Heavy",
                                //    test.Seo_Heavy >= reference.Seo_Heavy * 0.99 &&
                                //    test.Seo_Heavy <= reference.Seo_Heavy * 1.01
                                "Difference_SearchRT_FittedRT_Heavy",
                                    nanCasesForUnlabledData 0.9 1.1 (abs test.Difference_SearchRT_FittedRT_Heavy) (abs reference.Difference_SearchRT_FittedRT_Heavy)
                                "KLDiv_CorrectedObserved_Theoretical_Heavy",
                                    nanCasesForUnlabledData 0.9 1.1 test.KLDiv_CorrectedObserved_Theoretical_Heavy reference.KLDiv_CorrectedObserved_Theoretical_Heavy
                                "KLDiv_Observed_Theoretical_Heavy",
                                    nanCasesForUnlabledData 0.9 1.1 test.KLDiv_Observed_Theoretical_Heavy reference.KLDiv_Observed_Theoretical_Heavy
                                "Correlation_Light_Heavy",
                                    nanCasesForUnlabledData 0.99 1.01 test.Correlation_Light_Heavy reference.Correlation_Light_Heavy
                                "QuantificationSource",
                                    test.QuantificationSource = reference.QuantificationSource
                                "IonMobility",
                                    nanCasesForUnlabledData 0.99 1.01 test.IonMobility reference.IonMobility
                            |]
                        ) referenceQuantTIMs testQuantTIMs
                        |> Array.collect (
                            Array.filter (fun (field,equal) -> equal = false)
                        )
                        |> Array.distinct
                        |> Array.map fst
                    if unequalFields |> Array.isEmpty then true , ""
                    else false, unequalFields |> String.concat ";"
                Expect.isTrue compare (sprintf "Quants are different in the following fields: %s" fields)

        // sequenced, as both read minimal.mzlite, which the tool locks
        quantTIMsCase "PSMBasedQuantificationTIMs" "QuantificationTIMsParams.json" "" "minimalReference.quant" |> testSequenced
        // wavelet path with baseline correction, in its own folder
        quantTIMsCase "PSMBasedQuantificationTIMsWavelet" "QuantificationTIMsWaveletParams.json" "/wavelet" "minimalWaveletReference.quant" |> testSequenced
            

        testCase "QuantBasedAlignment" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let target = relToDirectory "../../../data/QuantBasedAlignment/in"
            let sources = relToDirectory "../../../data/QuantBasedAlignment/in"
            let outDirectory = relToDirectory "../../../data/QuantBasedAlignment/out"
            let relQuantPath =
                if System.Runtime.InteropServices.RuntimeInformation.IsOSPlatform(System.Runtime.InteropServices.OSPlatform.Linux) then
                    toolDll "QuantBasedAlignment_linux-x64"
                elif System.Runtime.InteropServices.RuntimeInformation.IsOSPlatform(System.Runtime.InteropServices.OSPlatform.Windows) then
                    toolDll "QuantBasedAlignment_win-x64"
                else    
                    failwith "not supported OS to test QuantBasedAlignment"                
            let quantExe = relQuantPath
            // cleanup
            try
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc1.align")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc1.alignmetric")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc1_log.txt")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/QuantBasedAlignment_log.txt")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc1_Metrics.html")   
                
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc2.align")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc2.alignmetric")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc2_log.txt")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc2_Metrics.html")   

                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc3.align")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc3.alignmetric")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc3_log.txt")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc3_Metrics.html")   
            with
            | _ -> ()
            // run tool
            runDotNet (sprintf "%s -i %s -ii %s -o %s" quantExe target sources outDirectory) baseDir
            let referenceAlign =
                let quantPath = relToDirectory "../../../data/QuantBasedAlignment/out/trunc1_Reference.align"
                FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.AlignmentResult>().ReadFile(quantPath,'\t',false,1)
                |> Array.ofSeq
            let testAlign =
                let quantPath = relToDirectory "../../../data/QuantBasedAlignment/out/trunc1.align"
                FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.AlignmentResult>().ReadFile(quantPath,'\t',false,1)
                |> Array.ofSeq
            let compare = 
                    referenceAlign
                    |> Array.filter (fun (reference: ProteomIQon.Dto.AlignmentResult) ->
                            testAlign
                            |> Array.exists (fun (test: ProteomIQon.Dto.AlignmentResult) ->
                                test.StringSequence = reference.StringSequence
                                test.GlobalMod = reference.GlobalMod
                                test.Charge = reference.Charge
                                test.PepSequenceID = reference.PepSequenceID
                                test.ModSequenceID = reference.ModSequenceID
                            )
                    )  
                    |> Array.length
                    |> (=) referenceAlign.Length           
            Expect.isTrue compare (sprintf "Quants are different in the following fields")


        testCase "QuantBasedAlignment_v2" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let quantFiles = relToDirectory "../../../data/QuantBasedAlignment/in"
            // its own folder, because the QuantBasedAlignment test writes files of the same names
            let outDirectory = relToDirectory "../../../data/QuantBasedAlignment_v2/out/test/"
            let alignExe = toolDll "QuantBasedAlignment_v2"
            runDotNet (sprintf "%s -i %s -ii %s -o %s" alignExe quantFiles quantFiles outDirectory) baseDir
            let compare =
                tsvFilesEqual 1e-6
                    (relToDirectory "../../../data/QuantBasedAlignment_v2/out/trunc1_Reference.align")
                    (relToDirectory "../../../data/QuantBasedAlignment_v2/out/test/trunc1.align")
            Directory.Delete(outDirectory, true)
            Expect.isTrue compare "Alignments are different"

        testCase "MzMLToMzLiteIonMobility" <| fun _ ->
            // the first 40 spectra (11 MS1, 29 MS2) of a timsTOF DDA run converted by msconvert: m/z, intensity and ion mobility as zlib
            // compressed 64 bit floats, the layout the raw reader of MzIO reads; spectrum 5 is stored uncompressed and spectrum 12 as
            // 32 bit floats, so the XML reader reads them
            let relToDirectory = getRelativePath baseDir
            let mzML = relToDirectory "../../../data/MzMLToMzLiteIonMobility/in/minimal.mzML"
            let convParams = relToDirectory "../../../data/MzMLToMzLiteIonMobility/in/mzMLToMzLiteParams.json"
            let outDirectory = relToDirectory "../../../data/MzMLToMzLiteIonMobility/out/"
            if Directory.Exists outDirectory then Directory.Delete(outDirectory, true)
            runDotNet (sprintf "%s -i %s -o %s -p %s" (toolDll "MzMLToMzLiteIonMobility") mzML outDirectory convParams) baseDir
            let expected = readMzMLSpectra mzML
            let converted = readMzliteSpectra (Path.Combine(outDirectory, "minimal.mzlite"))
            deleteDirectory outDirectory
            Expect.equal expected.Length 40 "the mzML holds 40 spectra"
            Expect.equal (converted |> Map.keys |> Set.ofSeq) (expected |> Array.map (fun s -> s.ID) |> Set.ofArray) "the mzlite holds the spectra of the mzML"
            for s in expected do
                let msLevel, scanTime, scanTimeUnit, precursorMz, peaks = converted.[s.ID]
                Expect.equal msLevel s.MsLevel (sprintf "%s: MS level" s.ID)
                Expect.equal s.ScanTimeUnit "UO:0000010" (sprintf "%s: the scan time of the mzML is in seconds" s.ID)
                Expect.equal scanTimeUnit (Some "UO:0000031") (sprintf "%s: the scan time of the mzlite is in minutes" s.ID)
                Expect.floatClose Accuracy.veryHigh scanTime (s.ScanTime / 60.) (sprintf "%s: scan time" s.ID)
                if s.MsLevel = 2 then Expect.equal precursorMz s.SelectedIonMz (sprintf "%s: precursor m/z" s.ID)
                // the peaks in the order of the file, values unchanged
                Expect.equal (peaks |> Array.map (fun p -> p.Mz)) s.Mz (sprintf "%s: m/z" s.ID)
                Expect.equal (peaks |> Array.map (fun p -> p.Intensity)) s.Intensity (sprintf "%s: intensities" s.ID)
                Expect.equal (peaks |> Array.map (fun p -> p.IonMobility)) (s.IonMobility |> Array.map Some) (sprintf "%s: ion mobilities" s.ID)

        testCase "ProteinInference" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let db = relToDirectory "../../../data/ProteinInference/in/Minimal.db"
            let qpsm = relToDirectory "../../../data/ProteinInference/in/minimal.qpsm"
            let proteinInferenceParamsStorey = relToDirectory "../../../data/ProteinInference/in/ProteinInferenceParamsStorey.json"
            let proteinInferenceParamsMAYU = relToDirectory "../../../data/ProteinInference/in/ProteinInferenceParamsMAYU.json"
            let outDirectoryStorey = relToDirectory "../../../data/ProteinInference/out/storeyOut"
            let outDirectoryMAYU = relToDirectory "../../../data/ProteinInference/out/mayuOut"
            let protInfExe = toolDll "ProteinInference"
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" protInfExe qpsm outDirectoryStorey proteinInferenceParamsStorey db) baseDir
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" protInfExe qpsm outDirectoryMAYU proteinInferenceParamsMAYU db) baseDir
            let compare =
                tsvFilesEqual 1e-9
                    (relToDirectory "../../../data/ProteinInference/out/storeyOut/minimalReference.prot")
                    (relToDirectory "../../../data/ProteinInference/out/storeyOut/minimal.prot") &&
                tsvFilesEqual 1e-9
                    (relToDirectory "../../../data/ProteinInference/out/mayuOut/minimalReference.prot")
                    (relToDirectory "../../../data/ProteinInference/out/mayuOut/minimal.prot")
            // cleanup
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/minimal.prot")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/ProteinInference_createClassItemCollection_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/ProteinInference_inferProteins_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/ProteinInference_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/ProteinInference_readAndInferFile_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/QValueGraph.html")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/minimal.prot")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/ProteinInference_createClassItemCollection_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/ProteinInference_inferProteins_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/ProteinInference_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/ProteinInference_readAndInferFile_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/QValueGraph.html")
            Expect.isTrue compare "Prots are different"
        testCase "LabelFreeProteinQuantification" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let quantAndProt = relToDirectory "../../../data/LabelFreeProteinQuantification/in/minimal.quantAndProt"
            let labelFreeQuantificationParams = "../../../data/LabelFreeProteinQuantification/in/LabelFreeQuantificationParams.json"
            let labelFreeQuantificationParamsChargeAgg = "../../../data/LabelFreeProteinQuantification/in/LabelFreeQuantificationParams_ChargeAgg.json"
            let labelFreeQuantificationParamsChargeAggModAgg = "../../../data/LabelFreeProteinQuantification/in/LabelFreeQuantificationParams_ChargeAgg_ModAgg.json"
            let labelFreeQuantificationParamsTransformFilterSum = "../../../data/LabelFreeProteinQuantification/in/LabelFreeQuantificationParams_Transform_Filter_Sum.json"
            let outDirectory = "../../../data/LabelFreeProteinQuantification/out/normal"
            let outDirectoryChargeAgg = "../../../data/LabelFreeProteinQuantification/out/chargeAgg"
            let outDirectoryChargeAggModAgg = "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg"
            let outDirectoryTransformFilterSum = "../../../data/LabelFreeProteinQuantification/out/transformFilterSum"
            let labelFreeExe = toolDll "LabelFreeProteinQuantification"
            runDotNet (sprintf "%s -i %s -o %s -p %s" labelFreeExe quantAndProt outDirectory labelFreeQuantificationParams) baseDir
            runDotNet (sprintf "%s -i %s -o %s -p %s" labelFreeExe quantAndProt outDirectoryChargeAgg labelFreeQuantificationParamsChargeAgg) baseDir
            runDotNet (sprintf "%s -i %s -o %s -p %s" labelFreeExe quantAndProt outDirectoryChargeAggModAgg labelFreeQuantificationParamsChargeAggModAgg) baseDir
            runDotNet (sprintf "%s -i %s -o %s -p %s" labelFreeExe quantAndProt outDirectoryTransformFilterSum labelFreeQuantificationParamsTransformFilterSum) baseDir
            let referenceLabelFree =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/minimalReference.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/minimalProteinReference.txt"
                File.ReadAllLines labelFreeProtein
            let labelFree =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/LabelFreeQuant.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/ProteinAggregation.txt"
                File.ReadAllLines labelFreeProtein
            let referenceLabelFreeChargeAgg =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/minimalReference.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/minimalProteinReference.txt"
                File.ReadAllLines labelFreeProtein,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/minimalChargeReference.txt"
                File.ReadAllLines labelFreeCharge
            let labelFreeChargeAgg =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/LabelFreeQuant.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/ProteinAggregation.txt"
                File.ReadAllLines labelFreeProtein,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/ChargeAggregation.txt"
                File.ReadAllLines labelFreeCharge
            let referenceLabelFreeChargeAggModAgg =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/minimalReference.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/minimalProteinReference.txt"
                File.ReadAllLines labelFreeProtein,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/minimalChargeReference.txt"
                File.ReadAllLines labelFreeCharge,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/minimalModificationReference.txt"
                File.ReadAllLines labelFreeCharge
            let labelFreeChargeAggModAgg =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/LabelFreeQuant.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ProteinAggregation.txt"
                File.ReadAllLines labelFreeProtein,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ChargeAggregation.txt"
                File.ReadAllLines labelFreeCharge,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ModificationAggregation.txt"
                File.ReadAllLines labelFreeCharge
            let referenceLabelFreeTransformFilterSum =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/minimalReference.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/minimalProteinReference.txt"
                File.ReadAllLines labelFreeProtein
            let labelFreeTransformFilterSum =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/LabelFreeQuant.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/ProteinAggregation.txt"
                File.ReadAllLines labelFreeProtein
            let compare = referenceLabelFree = labelFree && referenceLabelFreeChargeAgg = labelFreeChargeAgg && referenceLabelFreeChargeAggModAgg = labelFreeChargeAggModAgg && referenceLabelFreeTransformFilterSum = labelFreeTransformFilterSum
            // cleanup
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/LabelFreeQuant.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/LabelFreeQuant.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/ChargeAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/LabelFreeQuant.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ChargeAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ModificationAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/LabelFreeQuant.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/ProteinAggregation.txt")
            Expect.isTrue compare "Output files are not identical"

        testCase "LabeledProteinQuantification" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let quantAndProt = relToDirectory "../../../data/LabeledProteinQuantification/in/minimal.quantAndProt"
            let labeledQuantificationParams = "../../../data/LabeledProteinQuantification/in/LabeledQuantificationParams.json"
            let labeledQuantificationParamsCorrChargeAgg = "../../../data/LabeledProteinQuantification/in/LabeledQuantificationParams_CorrFilter_ChargeAgg.json"
            let labeledQuantificationParamsCorrChargeAggModAgg = "../../../data/LabeledProteinQuantification/in/LabeledQuantificationParams_CorrFilter_ChargeAgg_ModAgg.json"
            let labeledQuantificationParamsTransformFilterSum = "../../../data/LabeledProteinQuantification/in/LabeledQuantificationParams_Transform_Filter_Sum.json"
            let labeledQuantificationParamsCorr = "../../../data/LabeledProteinQuantification/in/LabeledQuantificationParams_withCorrFilter.json"
            let outDirectory = "../../../data/LabeledProteinQuantification/out/normal"
            let outDirectoryCorrChargeAgg = "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg"
            let outDirectoryCorrChargeAggModAgg = "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg"
            let outDirectoryTransformFilterSum = "../../../data/LabeledProteinQuantification/out/transformFilterSum"
            let outDirectoryCorr = "../../../data/LabeledProteinQuantification/out/corrFilter"
            let labeledExe = toolDll "LabeledProteinQuantification"
            runDotNet (sprintf "%s -i %s -o %s -p %s" labeledExe quantAndProt outDirectory labeledQuantificationParams) baseDir
            runDotNet (sprintf "%s -i %s -o %s -p %s" labeledExe quantAndProt outDirectoryCorrChargeAgg labeledQuantificationParamsCorrChargeAgg) baseDir
            runDotNet (sprintf "%s -i %s -o %s -p %s" labeledExe quantAndProt outDirectoryCorrChargeAggModAgg labeledQuantificationParamsCorrChargeAggModAgg) baseDir
            runDotNet (sprintf "%s -i %s -o %s -p %s" labeledExe quantAndProt outDirectoryTransformFilterSum labeledQuantificationParamsTransformFilterSum) baseDir
            runDotNet (sprintf "%s -i %s -o %s -p %s" labeledExe quantAndProt outDirectoryCorr labeledQuantificationParamsCorr) baseDir
            let referenceLabeled =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/minimalReference.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/minimalProteinReference.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/minimalGlobModReference.txt"
                File.ReadAllLines labeledGlobMod
            let labeled =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/LabeledQuant.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/ProteinAggregation.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/GlobModAggregation.txt"
                File.ReadAllLines labeledGlobMod
            let referenceLabeledCorrChargeAgg =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/minimalReference.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/minimalProteinReference.txt"
                File.ReadAllLines labeledProtein,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/minimalChargeReference.txt"
                File.ReadAllLines labeledCharge,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/minimalGlobModReference.txt"
                File.ReadAllLines labeledGlobMod
            let labeledCorrChargeAgg =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/LabeledQuant.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/ProteinAggregation.txt"
                File.ReadAllLines labeledProtein,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/ChargeAggregation.txt"
                File.ReadAllLines labeledCharge,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/GlobModAggregation.txt"
                File.ReadAllLines labeledGlobMod
            let referenceLabeledCorrChargeAggModAgg =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/minimalReference.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/minimalProteinReference.txt"
                File.ReadAllLines labeledProtein,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/minimalChargeReference.txt"
                File.ReadAllLines labeledCharge,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/minimalModificationReference.txt"
                File.ReadAllLines labeledCharge,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/minimalGlobModReference.txt"
                File.ReadAllLines labeledGlobMod
            let labeledCorrChargeAggModAgg =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/LabeledQuant.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ProteinAggregation.txt"
                File.ReadAllLines labeledProtein,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ChargeAggregation.txt"
                File.ReadAllLines labeledCharge,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ModificationAggregation.txt"
                File.ReadAllLines labeledCharge,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/GlobModAggregation.txt"
                File.ReadAllLines labeledGlobMod
            let referenceLabeledTransformFilterSum =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/minimalReference.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/minimalProteinReference.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/minimalGlobModReference.txt"
                File.ReadAllLines labeledGlobMod
            let labeledTransformFilterSum =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/LabeledQuant.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/ProteinAggregation.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/GlobModAggregation.txt"
                File.ReadAllLines labeledGlobMod
            let referenceLabeledCorr =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/minimalReference.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/minimalProteinReference.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/minimalGlobModReference.txt"
                File.ReadAllLines labeledGlobMod
            let labeledCorr =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/LabeledQuant.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/ProteinAggregation.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/GlobModAggregation.txt"
                File.ReadAllLines labeledGlobMod
            let compare = referenceLabeled = labeled && referenceLabeledCorrChargeAgg = labeledCorrChargeAgg && referenceLabeledCorrChargeAggModAgg = labeledCorrChargeAggModAgg && referenceLabeledTransformFilterSum = labeledTransformFilterSum && referenceLabeledCorr = labeledCorr
            // cleanup
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/normal/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/normal/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/normal/LabeledQuant.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/normal/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/normal/GlobModAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/LabeledQuant.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/GlobModAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/LabeledQuant.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/GlobModAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/ChargeAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/LabeledQuant.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/GlobModAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ChargeAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ModificationAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/LabeledQuant.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/GlobModAggregation.txt")
            Expect.isTrue compare "Output files are not identical"
        testCase "AddDeducedPeptides" <| fun _ ->
            let relToDirectory = getRelativePath baseDir
            let quantDirectory = relToDirectory "../../../data/AddDeducedPeptides/in/quant"
            let protDirectory = relToDirectory "../../../data/AddDeducedPeptides/in/prot"
            let outDirectory = relToDirectory "../../../data/AddDeducedPeptides/out"
            let addDeducedPeptidesExe = toolDll "AddDeducedPeptides"
            // run tool
            runDotNet (sprintf "%s -i %s -ii %s -o %s" addDeducedPeptidesExe quantDirectory protDirectory outDirectory) baseDir
            let reference1 =
                let qpsmPath = relToDirectory "../../../data/AddDeducedPeptides/out/minimalReference1.prot"
                File.ReadAllLines qpsmPath
            let test1 =
                let qpsmPath = relToDirectory "../../../data/AddDeducedPeptides/out/minimalQuant1.prot"
                File.ReadAllLines qpsmPath
            let reference2 =
                let qpsmPath = relToDirectory "../../../data/AddDeducedPeptides/out/minimalReference2.prot"
                File.ReadAllLines qpsmPath
            let test2 =
                let qpsmPath = relToDirectory "../../../data/AddDeducedPeptides/out/minimalQuant2.prot"
                File.ReadAllLines qpsmPath

            let compare =
                reference1 = test1 && reference2 = test2
            // cleanup
            File.Delete (relToDirectory "../../../data/AddDeducedPeptides/out/minimalQuant1.prot")
            File.Delete (relToDirectory "../../../data/AddDeducedPeptides/out/minimalQuant2.prot")
            File.Delete (relToDirectory "../../../data/AddDeducedPeptides/out/AddDeducedPeptides_log.txt")
            Expect.isTrue compare "Deduced Proteins are different"
    ]
