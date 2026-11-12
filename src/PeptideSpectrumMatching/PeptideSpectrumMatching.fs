namespace ProteomIQon

open System
open System.Collections.Generic
open System.IO
open BioFSharp
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open MzIO
open MzIO.IO
open MzIO.MzSQL
open MzIO.Processing
open Domain
open Core

/// Scores the MS2 spectra of a run against a peptide data base. Spectrum headers are parsed from the
/// description text one spectrum at a time, the peptide data base is read once into a mass sorted
/// table, the parsed candidates and their ion series are reused between neighbouring spectra, and
/// the result rows are written directly.
/// Conventions: spectra with equal scan times keep their file order, and candidates enter the
/// scoring by ascending RoundedMass, then ModSequenceID, so equal scores keep that order. The
/// X!Tandem-like scores of a row belong to the peptide of the row.
module PeptideSpectrumMatching =

    open BioFSharp.Mz.ChargeState

    /// Time spent per phase of a run. Every run has its own counters and is processed on one
    /// thread, so plain accumulators are enough.
    module Profile =
        let names = [| "headers"; "ms1 peaks"; "ms2 peaks"; "charge states"; "p values"; "psm peaks"; "candidates"; "sequest"; "xscoring"; "write" |]
        let create () : int64[] = Array.zeroCreate names.Length
        let add (ticks: int64[]) (phase: int) (start: int64) =
            ticks.[phase] <- ticks.[phase] + (Diagnostics.Stopwatch.GetTimestamp() - start)
        let summary (ticks: int64[]) =
            let seconds (t: int64) = float t / float Diagnostics.Stopwatch.Frequency
            names
            |> Array.mapi (fun i name -> sprintf "%s %.1f s" name (seconds ticks.[i]))
            |> String.concat ", "

    /// Charge state candidates of the MS2 spectra, each assigned to the last MS1 at or before its scan time.
    /// ms1s are (scan time, ID) sorted by scan time and ms2s sorted by scan time (SpectrumIndex.readMs1AndMs2).
    let getPrecursorCharge (chParams:ChargeState.ChargeDetermParams) rnd (ticks: int64[]) (ms1s: (float * string)[]) (ms2s: SpectrumIndex.Ms2Header[]) (inReader: IMzIODataReader) =
        let ms1ScanTimes = ms1s |> Array.map fst
        let ms2AssignedToMS1 =
            ms2s
            |> Array.choose (fun ms2 ->
                // the last MS1 at or before the MS2
                match SortedSearch.upperBound ms1ScanTimes ms2.ScanTime - 1 with
                | -1 -> None
                | ms1 -> Some (ms1, ms2))
            |> Array.groupBy fst
            |> Array.map (fun (ms1, ms1And2) -> ms1s.[ms1], ms1And2 |> Array.map snd)
        let ms2PossibleChargestates =
            ms2AssignedToMS1
            |> Array.map (fun ((_, ms1Id), ms2s) ->
                let t = Diagnostics.Stopwatch.GetTimestamp()
                let mzdata, intensityData = MzIO.Peaks.unzipIMzliteArray (inReader.ReadSpectrumPeaks(ms1Id).Peaks)
                Profile.add ticks 1 t
                ms2s
                |> Array.filter (fun ms2 ->
                    let t = Diagnostics.Stopwatch.GetTimestamp()
                    let hasPeaks = inReader.ReadSpectrumPeaks(ms2.Id).Peaks |> Seq.isEmpty = false
                    Profile.add ticks 2 t
                    hasPeaks)
                |> Array.map (fun ms2 ->
                    let t = Diagnostics.Stopwatch.GetTimestamp()
                    let assignedCharges = ChargeState.putativePrecursorChargeStatesBy chParams mzdata intensityData ms1Id ms2.Id ms2.PrecursorMz
                    Profile.add ticks 3 t
                    match assignedCharges with
                    | [] ->
                        [
                            for i = 2 to 3 do
                                let precursorMz = ms2.PrecursorMz
                                let mass = Mass.ofMZ precursorMz (float i)
                                let score = getScore 10 1 100.
                                yield createAssignedCharge ms1Id ms2.Id precursorMz i mass 100. score [0.] 0 0. (Set[]) (Some 1.)
                        ]
                    | _ -> assignedCharges))
            |> Array.filter (fun x -> Array.isEmpty x |> not)
            |> Array.concat
            |> Array.filter (fun assignedCharges -> assignedCharges <> [])
            |> List.ofArray
        let peakPosStdDev =
            ms2PossibleChargestates
            |> List.filter (fun assignedCharges -> assignedCharges.Head.PositionMetricPValue.IsNone)
            |> List.map (fun assignedCharges -> assignedCharges.Head)
            |> ChargeState.peakPosStdDevBy
        let precursorMzOf =
            let byId = Dictionary<string, float>()
            for ms2 in ms2s do byId.[ms2.Id] <- ms2.PrecursorMz
            fun (id: string) -> byId.[id]
        let init = ChargeState.initMzDevOfRndSpec rnd {chParams with ExpectedMaximumCharge=8} peakPosStdDev
        let t = Diagnostics.Stopwatch.GetTimestamp()
        let positionMetricScoredCharges =
            ms2PossibleChargestates
            |> List.map (fun assignedCharges ->
                let ms1ID = assignedCharges.Head.PrecursorSpecID
                let ms2ID = assignedCharges.Head.ProductSpecID
                assignedCharges
                |> List.map (fun putativeCharge ->
                    match putativeCharge.PositionMetricPValue with
                    | Some _ -> putativeCharge
                    | None ->
                        let pValue = ChargeState.empiricalPValueOfSim init (putativeCharge.SubSetLength, float putativeCharge.PrecCharge) putativeCharge.MZChargeDev
                        {putativeCharge with PositionMetricPValue = Some pValue})
                |> List.filter (fun testIt -> testIt.PositionMetricPValue.Value <= 0.05)
                |> ChargeState.removeSubSetsOfBestHit
                |> (fun charges ->
                    match charges with
                    | [] ->
                        [
                            for i = 2 to 3 do
                                let precMz = precursorMzOf ms2ID
                                let mass = Mass.ofMZ precMz (float i)
                                let score = ChargeState.getScore 10 1 100.
                                yield ChargeState.createAssignedCharge ms1ID ms2ID precMz i mass 100. score [0.] 0 0. (Set[]) None
                        ]
                    | _ -> charges))
        Profile.add ticks 4 t
        positionMetricScoredCharges, peakPosStdDev

    /// The ModSequence table sorted by RoundedMass, then by ID.
    type PeptideTable =
        {
            ModSequenceID : int[]
            PepSequenceID : int[]
            RealMass      : float[]
            RoundedMass   : int64[]
            Sequence      : string[]
            GlobalMod     : int[]
        }

    let loadPeptideTable (cn: System.Data.SQLite.SQLiteConnection) =
        let count =
            use cmd = new System.Data.SQLite.SQLiteCommand("SELECT COUNT(*) FROM ModSequence", cn)
            Convert.ToInt32(cmd.ExecuteScalar())
        let modSequenceID = Array.zeroCreate<int> count
        let pepSequenceID = Array.zeroCreate<int> count
        let realMass = Array.zeroCreate<float> count
        let roundedMass = Array.zeroCreate<int64> count
        let sequence = Array.zeroCreate<string> count
        let globalMod = Array.zeroCreate<int> count
        use cmd = new System.Data.SQLite.SQLiteCommand("SELECT ID, PepSequenceID, RealMass, RoundedMass, Sequence, GlobalMod FROM ModSequence ORDER BY RoundedMass, ID", cn)
        use reader = cmd.ExecuteReader()
        let mutable row = 0
        while reader.Read() do
            modSequenceID.[row] <- reader.GetInt32 0
            pepSequenceID.[row] <- reader.GetInt32 1
            realMass.[row] <- reader.GetDouble 2
            roundedMass.[row] <- reader.GetInt64 3
            sequence.[row] <- reader.GetString 4
            globalMod.[row] <- reader.GetInt32 5
            row <- row + 1
        if row <> count then failwithf "The ModSequence table returned %i rows, %i were expected." row count
        {
            ModSequenceID = modSequenceID
            PepSequenceID = pepSequenceID
            RealMass = realMass
            RoundedMass = roundedMass
            Sequence = sequence
            GlobalMod = globalMod
        }

    /// The rows [first, last) whose rounded mass lies between the rounded bounds, as the query
    /// "RoundedMass BETWEEN lower AND upper" selects them.
    let rowRange (table: PeptideTable) (lowerMass: float) (upperMass: float) =
        let lower = Convert.ToInt64(lowerMass * 1000000.)
        let upper = Convert.ToInt64(upperMass * 1000000.)
        let bound (atMost: bool) (value: int64) =
            let mutable lo = 0
            let mutable hi = table.RoundedMass.Length
            while lo < hi do
                let m = lo + ((hi - lo) >>> 1)
                let r = table.RoundedMass.[m]
                if (if atMost then r <= value else r < value) then lo <- m + 1 else hi <- m
            lo
        bound false lower, bound true upper

    /// The candidates of the mass windows of one spectrum together with their ion series. A
    /// candidate is parsed and fragmented once while it stays in reach of the ascending windows of
    /// a charge state: the spectra of a charge state are processed in precursor mass order, so a
    /// row below the lowest window of the current spectrum is not needed again until the next
    /// charge state, and is dropped.
    type CandidateCache(table: PeptideTable, parse: int -> string -> AminoAcids.AminoAcid list, calcIonSeries: AminoAcids.AminoAcid list -> Fragmentation.FragmentMasses) =
        let cache = Dictionary<int, LookUpResult<AminoAcids.AminoAcid> * Fragmentation.FragmentMasses>()
        let stale = ResizeArray<int>()
        member _.Clear() = cache.Clear()
        /// Drops the rows below the given row.
        member _.DropBelow (row: int) =
            stale.Clear()
            for key in cache.Keys do
                if key < row then stale.Add key
            for key in stale do cache.Remove key |> ignore
        /// The rows of a window in table order.
        member _.Window (first: int, last: int) =
            [
                for row = first to last - 1 do
                    match cache.TryGetValue row with
                    | true, entry -> yield entry
                    | _ ->
                        let bioSequence = parse table.GlobalMod.[row] table.Sequence.[row]
                        let lookUp = createLookUpResult table.ModSequenceID.[row] table.PepSequenceID.[row] table.RealMass.[row] table.RoundedMass.[row] table.Sequence.[row] bioSequence table.GlobalMod.[row]
                        let entry = lookUp, calcIonSeries lookUp.BioSequence
                        cache.[row] <- entry
                        yield entry
            ]

    /// Writes one result row, tab separated, in the field order of the record. Numbers are
    /// written culture independent, floats with the shortest text that reads back the same value.
    let private writeRow (w: TextWriter) (r: Dto.PeptideSpectrumMatchingResult) =
        let inv = Globalization.CultureInfo.InvariantCulture
        let sep () = w.Write '\t'
        let int (x: int) = w.Write(x.ToString(inv)); sep ()
        let float (x: float) = w.Write(x.ToString(inv)); sep ()
        w.Write r.PSMId; sep ()
        int r.GlobalMod
        int r.PepSequenceID
        int r.ModSequenceID
        int r.Label
        int r.ScanNr
        float r.ScanTime
        int r.Charge
        float r.PrecursorMZ
        float r.TheoMass
        float r.AbsDeltaMass
        int r.PeptideLength
        int r.MissCleavages
        float r.SequestScore
        float r.SequestNormDeltaBestToRest
        float r.SequestNormDeltaNext
        float r.AndroScore
        float r.AndroNormDeltaBestToRest
        float r.AndroNormDeltaNext
        float r.XtandemScore
        float r.XtandemNormDeltaBestToRest
        float r.XtandemNormDeltaNext
        w.Write r.StringSequence; sep ()
        float r.IonMobility
        float r.Hyperscore
        float r.Expectscore
        int r.MatchedIons
        w.Write(r.TotalIons.ToString(inv))
        w.WriteLine()

    /// Scores every spectrum at its assigned charges and writes the rows. Returns the number of
    /// spectrum and charge attempts and the number of those that failed.
    let psm (processParams:PeptideSpectrumMatchingParams) (ticks: int64[]) (table: PeptideTable) (candidates: CandidateCache) (scanTimeOf: string -> float) (reader: IMzIODataReader) (outFilePath: string) (ms2sAndAssignedCharges: AssignedCharge list list) =
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension outFilePath)
        let mutable attempts = 0
        let mutable failed = 0
        use resultWriter = new StreamWriter(outFilePath, false)
        Reflection.FSharpType.GetRecordFields(typeof<Dto.PeptideSpectrumMatchingResult>)
        |> Seq.map (fun field -> field.Name)
        |> String.concat "\t"
        |> resultWriter.WriteLine
        let ms2IDAssignedCharge =
            ms2sAndAssignedCharges
            |> List.mapi (fun i assignedCharges ->
                assignedCharges |> List.map (fun assCh -> i, assCh.ProductSpecID, assCh))
            |> List.concat
            |> List.groupBy (fun (ascendingID, ms2Id, assCH) -> assCH.PrecCharge)
            |> List.map (fun (ch, ms2IDassCHL) ->
                ch, ms2IDassCHL |> List.sortBy (fun (ascendingID, ms2ID, assCh) -> assCh.PutMass))
            |> List.sortBy fst

        ms2IDAssignedCharge
        |> List.iter (fun (ch, ms2IdAssCh) -> logger.Trace (sprintf "%i spectra with charge %i" ms2IdAssCh.Length ch))

        ms2IDAssignedCharge
        |> List.iter (fun (ch, ms2IdAssCh) ->
            logger.Trace (sprintf "%i spectra with charge %i processed" ms2IdAssCh.Length ch)
            candidates.Clear()
            ms2IdAssCh
            |> List.iteri (fun i (ascendingID, ms2Id, assCh) ->
                if i%10000 = 0 then logger.Trace (sprintf "%i" i)
                attempts <- attempts + 1
                try
                    let scanTime = scanTimeOf ms2Id
                    let t = Diagnostics.Stopwatch.GetTimestamp()
                    let recSpec =
                        MzIO.Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(ms2Id).Peaks)
                        |> fun (mzData, intensityData) -> PeakArray.zip mzData intensityData
                    Profile.add ticks 5 t
                    let scanRange =
                        let floorToClosest10 x =
                            Math.Floor(x / 10.) * 10.
                        let ceilToClosest10 x =
                            Math.Ceiling(x / 10.) * 10.
                        let low = Math.Max(0., System.Math.Round((recSpec |> Array.minBy (fun x -> x.Mz)).Mz, 0) |> floorToClosest10)
                        let top = Math.Round((recSpec |> Array.maxBy (fun x -> x.Mz)).Mz, 0) |> ceilToClosest10
                        low, top

                    let t = Diagnostics.Stopwatch.GetTimestamp()
                    let windowOf precMz =
                        let lowerMass, upperMass =
                            let massWithH2O = BioFSharp.Mass.ofMZ precMz (ch |> float)
                            Mass.rangePpm processParams.LookUpPPM massWithH2O
                        rowRange table lowerMass upperMass
                    let mzMinusOne = assCh.PrecursorMZ - (Mass.Table.PMassInU / (float ch))
                    let windowMinusOne = windowOf mzMinusOne
                    let window = windowOf assCh.PrecursorMZ
                    candidates.DropBelow (min (fst windowMinusOne) (fst window))
                    let theoSpecs =
                        (candidates.Window windowMinusOne) @ (candidates.Window window)
                        |> List.distinctBy (fun (lookUpResult, _) -> lookUpResult.ModSequenceID)
                    Profile.add ticks 6 t

                    let t = Diagnostics.Stopwatch.GetTimestamp()
                    let sequestTheoreticalSpecs = SequestLike.getTheoSpecs scanRange assCh.PrecCharge theoSpecs
                    let sequestLikeScored =
                        SequestLike.calcSequestScore scanRange recSpec scanTime assCh.PrecCharge
                            assCh.PrecursorMZ sequestTheoreticalSpecs ms2Id
                    Profile.add ticks 7 t

                    let bestTargetSequest =
                        sequestLikeScored
                        |> List.filter (fun (x:SearchEngineResult.SearchEngineResult<float>) -> x.IsTarget)
                        |> List.truncate 10
                        |> List.map (fun x -> (x.ModSequenceID, x.GlobalMod), x)
                        |> Map.ofList

                    let bestDecoySequest =
                        sequestLikeScored
                        |> List.filter (fun x -> not x.IsTarget)
                        |> List.truncate 10
                        |> List.map (fun x -> (x.ModSequenceID, x.GlobalMod), x)
                        |> Map.ofList

                    let t = Diagnostics.Stopwatch.GetTimestamp()
                    let andromedaTheorticalSpecs =
                        theoSpecs
                        |> List.filter (fun (lookUpResult, fragments) ->
                            bestTargetSequest |> Map.containsKey (lookUpResult.ModSequenceID, lookUpResult.GlobalMod) ||
                            bestDecoySequest  |> Map.containsKey (lookUpResult.ModSequenceID, lookUpResult.GlobalMod))
                        |> XScoring.getTheoSpecs scanRange assCh.PrecCharge

                    let andromedaLikeScored, xtandemScored =
                        XScoring.calcAndromedaAndXTandemScore processParams.AndromedaParams.PMinPMax scanRange processParams.AndromedaParams.MatchingIonTolerancePPM
                            recSpec scanTime assCh.PrecCharge assCh.PrecursorMZ andromedaTheorticalSpecs ms2Id
                    // Both lists are sorted by their own score, so the X!Tandem-like result of a
                    // peptide is found by its key, not by its position.
                    let xtandemOf =
                        xtandemScored
                        |> List.map (fun (x: SearchEngineResult.SearchEngineResult<float>) -> (x.ModSequenceID, x.GlobalMod, x.IsTarget), x)
                        |> dict
                    Profile.add ticks 8 t

                    let t = Diagnostics.Stopwatch.GetTimestamp()
                    let result =
                        andromedaLikeScored
                        |> List.map (fun (androRes:SearchEngineResult.SearchEngineResult<float>) ->
                            let xTandemRes = xtandemOf.[(androRes.ModSequenceID, androRes.GlobalMod, androRes.IsTarget)]
                            let label = if androRes.IsTarget then 1 else -1
                            let scanNr = ascendingID
                            let absDeltaMass = (androRes.TheoMass-androRes.MeasuredMass) |> abs
                            let bestSequest = if androRes.IsTarget then bestTargetSequest else bestDecoySequest
                            match Map.tryFind (androRes.ModSequenceID, androRes.GlobalMod) bestSequest with
                            | Some sequestScore ->
                                let res i : Dto.PeptideSpectrumMatchingResult =
                                    {
                                        PSMId                        = androRes.SpectrumID.Replace(' ', '-') + "_" + ascendingID.ToString() + "_" + ch.ToString() + "_" + i.ToString()
                                        GlobalMod                    = androRes.GlobalMod
                                        PepSequenceID                = androRes.PepSequenceID
                                        ModSequenceID                = androRes.ModSequenceID
                                        Label                        = label
                                        ScanNr                       = scanNr
                                        ScanTime                     = androRes.ScanTime
                                        Charge                       = ch
                                        PrecursorMZ                  = androRes.PrecursorMZ
                                        TheoMass                     = androRes.TheoMass
                                        AbsDeltaMass                 = absDeltaMass
                                        PeptideLength                = androRes.PeptideLength
                                        MissCleavages                = -1
                                        SequestScore                 = sequestScore.Score
                                        SequestNormDeltaBestToRest   = sequestScore.NormDeltaBestToRest
                                        SequestNormDeltaNext         = sequestScore.NormDeltaNext
                                        AndroScore                   = androRes.Score
                                        AndroNormDeltaBestToRest     = androRes.NormDeltaBestToRest
                                        AndroNormDeltaNext           = androRes.NormDeltaNext
                                        XtandemScore                 = xTandemRes.Score
                                        XtandemNormDeltaBestToRest   = xTandemRes.NormDeltaBestToRest
                                        XtandemNormDeltaNext         = xTandemRes.NormDeltaNext
                                        StringSequence               = androRes.StringSequence
                                        IonMobility                  = nan
                                        Hyperscore                   = nan
                                        Expectscore                  = nan
                                        MatchedIons                  = 0
                                        TotalIons                    = 0
                                    }
                                Some (label, res)
                            | None -> None)
                        |> List.choose id
                        |> List.groupBy (fun x -> fst x)
                        |> List.sortBy fst
                        |> List.map (fun (x, y) ->
                            y
                            |> List.map snd
                            |> List.mapi (fun i x -> x i))
                        |> List.concat
                    result |> List.iter (writeRow resultWriter)
                    Profile.add ticks 9 t
                with
                | _ as ex ->
                    failed <- failed + 1
                    logger.Warn (sprintf "spec with id: %s at charge %i fails with: %A" ms2Id ch ex)))
        attempts, failed

    let scoreSpectra (processParams:PeptideSpectrumMatchingParams) (outputDir:string) (sdbParams: SearchDbParams) (table: PeptideTable) (instrumentOutput:string) =

        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension instrumentOutput)
        let stopwatch = Diagnostics.Stopwatch.StartNew()
        let ticks = Profile.create ()

        logger.Trace (sprintf "Input file: %s" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)

        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".psm"
            Path.Combine [|outputDir;fileName|]
        logger.Trace (sprintf "Result file path: %s" outFilePath)

        logger.Trace "Prepare processing functions."
        let chargeParams = processParams.ChargeStateDeterminationParams
        logger.Trace (sprintf "Charge parameters: %A" chargeParams)
        logger.Trace (sprintf "DB parameters: %A" sdbParams)
        let calcIonSeries aal =
            Fragmentation.Series.fragmentMassesPseudoReversed Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList sdbParams.MassFunction aal
        let parse = SearchDB.initOfModAminoAcidString sdbParams.IsotopicMod (sdbParams.FixedMods@sdbParams.VariableMods)
        let candidates = CandidateCache(table, parse, calcIonSeries)
        let rnd = new System.Random()
        logger.Trace "Finished preparing processing functions."

        logger.Trace "Init connection to input data base."
        let inReader = Core.MzIO.Reader.getReader instrumentOutput
        try
            Core.MzIO.Reader.openConnection inReader
            let inRunID = Core.MzIO.Reader.getDefaultRunID inReader
            logger.Trace (sprintf "Run ID: %s" inRunID)

            let t = Diagnostics.Stopwatch.GetTimestamp()
            let ms1s, ms2s = SpectrumIndex.readMs1AndMs2 inReader inRunID (fun msg -> logger.Warn msg) (fun msg -> logger.Warn msg)
            Profile.add ticks 0 t
            logger.Trace (sprintf "%i MS1 and %i MS2 spectrum headers read, %.1f s." ms1s.Length ms2s.Length stopwatch.Elapsed.TotalSeconds)
            if ms1s.Length = 0 then
                logger.Warn "The run has no MS1 spectrum with an ID and a finite scan time, so no charge state is determined and no spectrum is searched."
            let scanTimeOf =
                let byId = Dictionary<string, float>()
                for h in ms2s do byId.[h.Id] <- h.ScanTime
                fun (id: string) -> byId.[id]

            logger.Trace "Starting charge state determination."
            let ms2sAndAssignedCharges, peakPosStdDev = getPrecursorCharge chargeParams rnd ticks ms1s ms2s inReader
            logger.Trace (sprintf "Finished charge state determination, %.1f s." stopwatch.Elapsed.TotalSeconds)
            logger.Trace (sprintf "peak position standard deviation: %f" peakPosStdDev)

            logger.Trace "Starting peptide spectrum matching."
            let attempts, failed = psm processParams ticks table candidates scanTimeOf inReader outFilePath ms2sAndAssignedCharges
            logger.Trace (sprintf "Finished peptide spectrum matching, %.1f s." stopwatch.Elapsed.TotalSeconds)
            if failed > 0 then logger.Warn (sprintf "%i of %i spectrum and charge attempts failed and have no rows." failed attempts)
            logger.Trace (sprintf "Time by phase: %s." (Profile.summary ticks))
        finally
            inReader.Dispose()
        logger.Trace "Done."
