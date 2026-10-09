namespace ProteomIQon

open System
open System.IO
open System.Collections.Generic
open System.Globalization
open BioFSharp
open BioFSharp.Mz
open MzIO
open MzIO.IO
open MzIO.MzSQL
open MzIO.Model
open MzIO.Model.CvParam
open MzIO.MetaData.PSIMSExtension
open MzIO.Processing
open FSharpAux.IO
open ProteomIQon.Domain
open ProteomIQon.Core
open PeptideIndex
open SpectrumSearch

/// Peptide spectrum matching for ion mobility runs. The peptide database is turned into a
/// fragment ion index once per invocation. Every MS2 spectrum of a run is searched against it
/// and the best candidates then go through the scoring functions of PeptideSpectrumMatching
/// (SEQUEST-like and Andromeda-like, which also yields the X!Tandem-like score), so the result
/// has the layout of PeptideSpectrumMatching and PSMStatistics reads it directly. A run is
/// searched on one thread, parallelization happens on the file level.
module PeptideSpectrumMatchingTIMs =

    /// Everything the search needs to know about one MS2 spectrum.
    type Ms2Spectrum =
        {
            Id          : string
            ScanNr      : int
            ScanTime    : float
            PrecursorMz : float
            /// Charge state from the file, 0 when unknown.
            Charge      : int
            IonMobility : float
            Mz          : float[]
            Intensity   : float[]
        }

    /// The BioFSharp.Mz scoring functions bound to the database parameters.
    type ClassicScoring =
        {
            CalcIonSeries         : AminoAcids.AminoAcid list -> Fragmentation.FragmentMasses
            ParseSequence         : int -> string -> AminoAcids.AminoAcid list
            ParseRank             : int -> AminoAcids.AminoAcid list
            AndromedaPMinPMax     : int * int
            AndromedaTolerancePpm : float
        }

    /// What a search of one spectrum at one charge state produced.
    type SearchOutcome =
        | Rows of int
        | TooFewPeaks
        | NoCandidates
        | NoHits
        | Failed

    let toSearchSettings (processParams: PeptideSpectrumMatchingTIMsParams) : SearchSettings =
        {
            PrecursorTolerancePpm = processParams.PrecursorTolerancePPM
            FragmentTolerancePpm  = processParams.FragmentTolerancePPM
            IsotopeErrors         = Array.ofList processParams.IsotopeErrors
            MaxFragmentCharge     = processParams.MaxFragmentCharge
            TopNPeaks             = processParams.TopNPeaks
            MinimumRatio          = processParams.MinimumPeakRatio
            RemovePrecursorRange  = processParams.RemovePrecursorRange
            Deisotope             = processParams.Deisotope
            MinimumPeaks          = processParams.MinimumPeaks
            MinMatchedFragments   = processParams.MinMatchedFragments
            MinFragmentsModelling = processParams.MinFragmentsModelling
            ReportedHitsPerLabel  = processParams.ReportedHitsPerLabel
        }

    let toClassicScoring (processParams: PeptideSpectrumMatchingTIMsParams) (sdbParams: SearchDB.SearchDbParams) (table: PeptideTable) : ClassicScoring =
        let parse = SearchDB.initOfModAminoAcidString sdbParams.IsotopicMod (sdbParams.FixedMods @ sdbParams.VariableMods)
        let residues =
            table.Codes |> Array.map (fun code ->
                code.Masses |> Array.mapi (fun globalMod _ ->
                    match parse globalMod code.Token with
                    | [aa] -> aa
                    | _ -> failwithf "Invalid residue token %s" code.Token))
        // MassFunction is pure for one database's chemistry. Preserve its exact results,
        // including labelled neutral-loss objects, without recalculating their formulas.
        let masses = System.Collections.Concurrent.ConcurrentDictionary<IBioItem, float>()
        let compute = Func<IBioItem, float>(sdbParams.MassFunction)
        let massOf item = masses.GetOrAdd(item, compute)
        {
            CalcIonSeries = Fragmentation.Series.fragmentMasses processParams.nTerminalSeries processParams.cTerminalSeries massOf
            ParseSequence = parse
            ParseRank = fun rank ->
                let start = table.Start rank
                let globalMod = int (table.GlobalMod rank)
                List.init (table.Length rank) (fun i -> residues.[int table.Residues.[start + i]].[globalMod])
            AndromedaPMinPMax = processParams.AndromedaParams.PMinPMax
            AndromedaTolerancePpm = processParams.AndromedaParams.MatchingIonTolerancePPM
        }

    /// Value of a cv param as float, parsed culture independent.
    let private cvDouble (value: obj option) =
        match value with
        | Some (:? IParamBase<IConvertible> as p) ->
            tryGetValue p |> Option.map (fun v -> Convert.ToDouble(v, CultureInfo.InvariantCulture))
        | _ -> None

    /// Precursor charge of the first selected ion, 0 when the file does not carry one.
    let precursorCharge (ms: MassSpectrum) =
        ms.Precursors.GetProperties false
        |> Seq.tryPick (fun kv ->
            match kv.Value with
            | :? Precursor as p ->
                p.SelectedIons.GetProperties false
                |> Seq.tryPick (fun si ->
                    match si.Value with
                    | :? SelectedIon as ion -> cvDouble (ion.TryGetValue PSIMS_Precursor.ChargeState)
                    | _ -> None)
            | _ -> None)
        |> Option.map int
        |> Option.defaultValue 0

    /// Inverse reduced ion mobility (MS:1002815) of the scan, NaN when absent.
    let ionMobility (ms: MassSpectrum) =
        ms.Scans.GetProperties false
        |> Seq.tryPick (fun kv ->
            match kv.Value with
            | :? Scan as scan -> cvDouble (scan.TryGetValue "MS:1002815")
            | _ -> None)
        |> Option.defaultValue nan

    /// Time spent per phase of a run, for finding out where the time goes. Every run has its own
    /// counters. A run is searched on one thread at a time, which counts into the counters of that
    /// run, so plain accumulators are enough.
    module Profile =
        let names = [| "headers"; "peaks"; "preprocess"; "windows"; "scatter"; "collect"; "expectation"; "classic"; "rescore"; "write"; "lookups"; "sequest"; "xscoring"; "parse"; "ionseries" |]
        let create () : int64[] = Array.zeroCreate names.Length
        /// The counters of the run this thread works on.
        let private current = new Threading.ThreadLocal<int64[]>(fun () -> create ())
        /// Runs f with the given counters as the counters of this thread.
        let countInto (ticks: int64[]) (f: unit -> 'T) =
            let previous = current.Value
            current.Value <- ticks
            try f () finally current.Value <- previous
        let add (phase: int) (start: int64) =
            let ticks = current.Value
            ticks.[phase] <- ticks.[phase] + (Diagnostics.Stopwatch.GetTimestamp() - start)
        let summary (ticks: int64[]) =
            let seconds (t: int64) = float t / float Diagnostics.Stopwatch.Frequency
            names
            |> Array.mapi (fun i name -> sprintf "%s %.1f s" name (seconds ticks.[i]))
            |> String.concat ", "

    /// The few values the search needs from a spectrum description. Only these are kept, one
    /// description at a time is read, and Position is the place of the spectrum in the file, which
    /// orders spectra that share a scan time.
    type Ms2Header =
        {
            Id          : string
            Position    : int
            ScanTime    : float
            PrecursorMz : float
            Charge      : int
            IonMobility : float
        }

    /// The MS2 headers of a run in scan time order, equal scan times in file order, which the
    /// reader delivers by row id. The descriptions are deserialized one at a time, so the metadata
    /// of the run never exists as a whole. A description that cannot be read is logged and skipped.
    let private readMs2HeadersOf (reader: IMzIODataReader) (runId: string) (log: string -> unit) =
        let started = Diagnostics.Stopwatch.GetTimestamp()
        let descriptions =
            match reader with
            | :? MzSQL as sql ->
                seq {
                    use cmd = new System.Data.SQLite.SQLiteCommand("SELECT Description FROM Spectrum WHERE RunID = @runID ORDER BY rowid", sql.Connection)
                    cmd.Parameters.AddWithValue("@runID", runId) |> ignore
                    use rows = cmd.ExecuteReader()
                    while rows.Read() do yield SpectrumDescription.tryParse (rows.GetString 0)
                }
            | _ -> failwith "PeptideSpectrumMatchingTIMs reads mzlite files, which MzIO serves one spectrum at a time."
        let headers = ResizeArray<Ms2Header>()
        let mutable position = 0
        for h in descriptions do
            if not (isNull h.Error) then
                log (sprintf "spectrum at position %i cannot be read: %s" position h.Error)
            elif h.MsLevel = Some 2 &&
                 (String.IsNullOrWhiteSpace h.ID ||
                  not (h.ScanTime |> Option.exists (fun x -> Double.IsFinite x && x >= 0.)) ||
                  not (h.PrecursorMz |> Option.exists (fun x -> Double.IsFinite x && x > protonMass)) ||
                  (h.ChargeState |> Option.exists (fun x -> x < 0))) then
                log (sprintf "spectrum at position %i has missing or invalid search metadata; skipped." position)
            elif h.MsLevel = Some 2 then
                headers.Add
                    {
                        Id = h.ID
                        Position = position
                        ScanTime = defaultArg h.ScanTime -1.
                        PrecursorMz = defaultArg h.PrecursorMz -1.
                        Charge = defaultArg h.ChargeState 0
                        IonMobility = defaultArg h.IonMobility nan
                    }
            position <- position + 1
        let ordered = headers.ToArray()
        Array.sortInPlaceWith (fun a b -> compare (a.ScanTime, a.Position) (b.ScanTime, b.Position)) ordered
        Profile.add 0 started
        ordered

    /// MS2 spectra of the run in scan time order. ScanNr is the position in that order. The peaks
    /// of a spectrum are read when the sequence reaches it and are dropped again afterwards. A
    /// spectrum whose peaks cannot be read is logged and skipped.
    let readMs2Headers (reader: IMzIODataReader) (runId: string) (log: string -> unit) = readMs2HeadersOf reader runId log

    let validateMobilityMergePpm ppm =
        if not (Double.IsFinite ppm) || ppm < 0. then
            invalidArg "MobilityMergePpm" "Must be finite and nonnegative."

    /// Optional m/z aggregation after exact-coordinate summation. Each cluster spans at
    /// most ppm relative to its lowest m/z, preventing a chain of neighbors from bridging
    /// distant peaks. Return the summed intensity and intensity-weighted m/z.
    let mergeNearbyPeaks ppm (mz: float[]) (intensity: float[]) =
        validateMobilityMergePpm ppm
        if ppm = 0. then mz, intensity
        else
            let result = ResizeArray<float * float>()
            let mutable i = 0
            while i < mz.Length do
                let anchor = mz.[i]
                let upper = anchor + anchor * ppm * 1e-6
                let mutable sum = 0.
                let mutable weightedOffset = 0.
                while i < mz.Length && mz.[i] <= upper do
                    sum <- sum + intensity.[i]
                    weightedOffset <- weightedOffset + (mz.[i] - anchor) * intensity.[i]
                    i <- i + 1
                result.Add(anchor + weightedOffset / sum, sum)
            result.ToArray() |> Array.unzip

    /// Integrate the mobility scans already contained in one stored spectrum. As in TIMs
    /// quantification, equal m/z coordinates are summed. A positive ppm optionally combines
    /// nearby coordinates too; separate spectra are never combined. Runs before all scoring.
    let searchPeakArraysWithTolerance ppm (peaks: global.MzIO.Commons.Arrays.IMzIOArray<global.MzIO.Binary.Peak1D>) =
        validateMobilityMergePpm ppm
        if peaks |> Seq.exists (fun peak -> peak.IonMobility.IsSome) then
            let sums = Dictionary<float, float>()
            for peak in peaks do
                if Double.IsFinite peak.Mz && peak.Mz > 0. && Double.IsFinite peak.Intensity && peak.Intensity > 0. then
                    match sums.TryGetValue peak.Mz with
                    | true, intensity -> sums.[peak.Mz] <- intensity + peak.Intensity
                    | _ -> sums.Add(peak.Mz, peak.Intensity)
            sums
            |> Seq.map (fun pair -> pair.Key, pair.Value)
            |> Seq.sortBy fst
            |> Seq.toArray
            |> Array.unzip
            |> fun (mz, intensity) -> mergeNearbyPeaks ppm mz intensity
        else MzIO.Peaks.unzipIMzliteArray peaks

    /// Standalone-reader default follows the library's standard MS2 fragment tolerance.
    [<Literal>]
    let DefaultMobilityMergePpm = 20.

    let searchPeakArrays peaks = searchPeakArraysWithTolerance DefaultMobilityMergePpm peaks

    /// A failed read keeps its header position so every charge attempt can be failed.
    let private readMs2SpectrumUsing readPeaks (reader: IMzIODataReader) (headers: Ms2Header[]) (i: int) =
        let header = headers.[i]
        try
            let started = Diagnostics.Stopwatch.GetTimestamp()
            let mz, intensity = readPeaks (reader.ReadSpectrumPeaks(header.Id).Peaks)
            Profile.add 1 started
            Ok { Id = header.Id; ScanNr = i; ScanTime = header.ScanTime; PrecursorMz = header.PrecursorMz
                 Charge = header.Charge; IonMobility = header.IonMobility; Mz = mz; Intensity = intensity }
        with ex -> Error ex

    let readMs2SpectrumWithTolerance ppm reader headers i =
        readMs2SpectrumUsing (searchPeakArraysWithTolerance ppm) reader headers i

    /// A 1D search must preserve intensity when the same m/z is observed in multiple
    /// mobility scans of the stored spectrum. Drivers supply their fragment tolerance;
    /// this standalone reader uses the standard 20 ppm default.
    let readMs2Spectrum reader headers i = readMs2SpectrumUsing searchPeakArrays reader headers i

    /// Scan range as PeptideSpectrumMatching defines it: the spectrum limits rounded outwards to tens.
    let private scanRangeOf (mz: float[]) =
        let floorToClosest10 x = Math.Floor(x / 10.) * 10.
        let ceilToClosest10 x = Math.Ceiling(x / 10.) * 10.
        Math.Max(0., Math.Round(Array.min mz, 0) |> floorToClosest10), Math.Round(Array.max mz, 0) |> ceilToClosest10

    /// Scores the given peptides with the BioFSharp.Mz functions that PeptideSpectrumMatching
    /// uses. Every peptide yields a target and a reversed decoy result, keyed by
    /// (ModSequenceID, GlobalMod, isTarget).
    let classicScores (classic: ClassicScoring) (table: PeptideTable) (spectrum: Ms2Spectrum) (charge: int) (ranks: int[]) =
        let t = Diagnostics.Stopwatch.GetTimestamp()
        let theoSpecs =
            ranks
            |> List.ofArray
            |> List.map (fun rank ->
                let sequence = table.SequenceString rank
                let globalMod = int (table.GlobalMod rank)
                let mass = table.Mass.[rank]
                let t1 = Diagnostics.Stopwatch.GetTimestamp()
                let parsed = classic.ParseRank rank
                Profile.add 13 t1
                let lookUp =
                    SearchDB.createLookUpResult (table.ModSequenceId rank) (table.PepSequenceId rank) mass
                        (int64 (Math.Round(mass * 1000000.))) sequence parsed globalMod
                let t2 = Diagnostics.Stopwatch.GetTimestamp()
                let series = classic.CalcIonSeries lookUp.BioSequence
                Profile.add 14 t2
                lookUp, series)
        Profile.add 10 t
        let recSpec = PeakArray.zip spectrum.Mz spectrum.Intensity
        let scanRange = scanRangeOf spectrum.Mz
        let t = Diagnostics.Stopwatch.GetTimestamp()
        let sequest =
            SequestLike.getTheoSpecs scanRange charge theoSpecs
            |> fun t -> SequestLike.calcSequestScore scanRange recSpec spectrum.ScanTime charge spectrum.PrecursorMz t spectrum.Id
        Profile.add 11 t
        let t = Diagnostics.Stopwatch.GetTimestamp()
        let andro, xtandem =
            XScoring.getTheoSpecs scanRange charge theoSpecs
            |> fun t -> XScoring.calcAndromedaAndXTandemScore classic.AndromedaPMinPMax scanRange classic.AndromedaTolerancePpm recSpec spectrum.ScanTime charge spectrum.PrecursorMz t spectrum.Id
        Profile.add 12 t
        let toMap (results: SearchEngineResult.SearchEngineResult<float> list) =
            results |> List.map (fun r -> (r.ModSequenceID, r.GlobalMod, r.IsTarget), r) |> dict
        toMap sequest, toMap andro, toMap xtandem

    /// What one index slice found for one spectrum at one charge: the index pass hyperscores of
    /// the candidates that enter the expectation model and the best candidates of both labels.
    type SliceResult =
        | SliceTooFewPeaks
        | SliceNoCandidates
        | SliceFound of histogram: int[] * targets: Candidate[] * decoys: Candidate[]

    /// The index pass of one slice for one spectrum at one charge.
    let indexPass (settings: SearchSettings) (table: PeptideTable) (index: FragmentIndex) (scratch: Scratch)
                  (maxFragments: int) (precursorMz: float) (processed: ProcessedSpectrum) (charge: int) (rankLo: int) (rankHi: int) =
        if processed.Mz.Length < settings.MinimumPeaks then SliceTooFewPeaks
        else
            let precursorMass = (precursorMz - protonMass) * float charge
            let t = Diagnostics.Stopwatch.GetTimestamp()
            let windows = clipWindows (precursorWindows settings table precursorMass) rankLo rankHi
            Profile.add 3 t
            if windows.Length = 0 then SliceNoCandidates
            else
                let t = Diagnostics.Stopwatch.GetTimestamp()
                scratch.Reset (slotCount windows) maxFragments
                scatter settings index scratch windows charge processed
                Profile.add 4 t
                let t = Diagnostics.Stopwatch.GetTimestamp()
                let refine c =
                    let tr = Diagnostics.Stopwatch.GetTimestamp()
                    let nb, ny, sumB, sumY, total = rescore settings table scratch charge processed c.Rank c.IsDecoy
                    Profile.add 8 tr
                    { c with Exact = hyperscore nb ny sumB sumY; ExactMatched = nb + ny; Total = total }
                let scores, targets, decoys = collect settings scratch windows refine
                Profile.add 5 t
                SliceFound (scores, targets, decoys)

    /// The final scoring of one spectrum at one charge from what the slices kept: the expectation
    /// model from the histogram, the classic scores of the reported candidates, and the rows.
    let finishSpectrum (settings: SearchSettings) (classic: ClassicScoring) (table: PeptideTable)
                       (spectrum: Ms2Spectrum) (charge: int) (histogram: int[])
                       (targets: Candidate[]) (decoys: Candidate[]) =
        let precursorMass = (spectrum.PrecursorMz - protonMass) * float charge
        let t = Diagnostics.Stopwatch.GetTimestamp()
        let expectation = expectationModelOfHistogram histogram
        Profile.add 6 t
        // The classic scores and their deltas are computed over every kept candidate, the ones
        // the MinMatchedFragments filter drops from the rows included, as the one pass search
        // did: the deltas measure the distance to the next peptide of the window, whether or
        // not that peptide gets a row. The kept set is bounded and does not depend on the slicing.
        // candidateOrder already defines a slice-independent input order. Preserve that
        // convention: another rank sort changes which tied classic result gets delta-next.
        let ranks = Array.append targets decoys |> Array.map (fun c -> c.Rank) |> Array.distinct
        // The rows of a label in the order of the exact hyperscore. Equal exact scores keep the
        // index pass order, then the lower rank comes first.
        let select (candidates: Candidate[]) =
            candidates
            |> Array.filter (fun c -> c.ExactMatched >= settings.MinMatchedFragments)
            |> Array.sortWith (fun a b ->
                let byExact = compare b.Exact a.Exact
                if byExact <> 0 then byExact
                else
                    let byIndex = compare b.Hyperscore a.Hyperscore
                    if byIndex <> 0 then byIndex else compare a.Rank b.Rank)
        let targets, decoys = select targets, select decoys
        let t = Diagnostics.Stopwatch.GetTimestamp()
        let sequestMap, androMap, xtandemMap = classicScores classic table spectrum charge ranks
        Profile.add 7 t
        let score (m: IDictionary<_, SearchEngineResult.SearchEngineResult<float>>) key =
            match m.TryGetValue key with
            | true, r -> r.Score, r.NormDeltaBestToRest, r.NormDeltaNext
            | _ -> 0., 0., 0.
        // The exact rescoring decides the reported hyperscore, the rank and the
        // MinMatchedFragments filter. The expectation value belongs to the index pass
        // score, the quantity the survival model was built from.
        let makeRows (candidates: Candidate[]) (label: int) =
            candidates
            |> Array.mapi (fun i c ->
                let key = (table.ModSequenceId c.Rank, int (table.GlobalMod c.Rank), not c.IsDecoy)
                let sq, sqRest, sqNext = score sequestMap key
                let an, anRest, anNext = score androMap key
                let xt, xtRest, xtNext = score xtandemMap key
                let theoMass = table.Mass.[c.Rank]
                // Mass error after the isotope error that fits the candidate best.
                let absDeltaMass =
                    settings.IsotopeErrors
                    |> Array.map (fun k -> abs (theoMass - (precursorMass - float k * isotopeSpacing)))
                    |> Array.min
                let row : Dto.PeptideSpectrumMatchingResult =
                    {
                        PSMId = spectrum.Id.Replace(' ', '-') + "_" + string spectrum.ScanNr + "_" + string charge + "_" + string i
                        GlobalMod = int (table.GlobalMod c.Rank)
                        PepSequenceID = table.PepSequenceId c.Rank
                        ModSequenceID = table.ModSequenceId c.Rank
                        Label = label
                        ScanNr = spectrum.ScanNr
                        ScanTime = spectrum.ScanTime
                        Charge = charge
                        PrecursorMZ = spectrum.PrecursorMz
                        TheoMass = theoMass
                        AbsDeltaMass = absDeltaMass
                        PeptideLength = table.Length c.Rank
                        MissCleavages = -1
                        SequestScore = sq
                        SequestNormDeltaBestToRest = sqRest
                        SequestNormDeltaNext = sqNext
                        AndroScore = an
                        AndroNormDeltaBestToRest = anRest
                        AndroNormDeltaNext = anNext
                        XtandemScore = xt
                        XtandemNormDeltaBestToRest = xtRest
                        XtandemNormDeltaNext = xtNext
                        StringSequence = table.SequenceString c.Rank
                        IonMobility = spectrum.IonMobility
                        Hyperscore = c.Exact
                        Expectscore = expectation c.Hyperscore
                        MatchedIons = c.ExactMatched
                        TotalIons = c.Total
                    }
                row)
        let rows = Array.append (makeRows targets 1) (makeRows decoys -1)
        if rows.Length = 0 then NoHits, [||] else Rows rows.Length, rows

    /// Loads the peptide table once for all runs. The index holds monoisotopic fragment masses,
    /// so a database with average masses is refused.
    let prepareTable (processParams: PeptideSpectrumMatchingTIMsParams) (cn: System.Data.SQLite.SQLiteConnection) (log: string -> unit) =
        let sdbParams = SearchDB.getSDBParamsByCn cn
        match sdbParams.MassMode with
        | SearchDB.MassMode.Monoisotopic -> ()
        | SearchDB.MassMode.Average -> failwith "The peptide data base uses average masses. PeptideSpectrumMatchingTIMs matches monoisotopic fragment masses and needs a monoisotopic data base."
        let table = loadPeptideTable cn sdbParams log
        let deviation = verifyMasses table 0.001
        log (sprintf "Residue masses reproduce the database masses, worst deviation %g Da." deviation)
        sdbParams, table

    /// Writes one result row as the tab separated line the CSV writer of FSharpAux produces, every
    /// value through its own ToString, in the field order of the record.
    let private writeRow (w: IO.TextWriter) (r: Dto.PeptideSpectrumMatchingResult) =
        let sep () = w.Write '\t'
        w.Write r.PSMId; sep ()
        w.Write(r.GlobalMod.ToString()); sep ()
        w.Write(r.PepSequenceID.ToString()); sep ()
        w.Write(r.ModSequenceID.ToString()); sep ()
        w.Write(r.Label.ToString()); sep ()
        w.Write(r.ScanNr.ToString()); sep ()
        w.Write(r.ScanTime.ToString()); sep ()
        w.Write(r.Charge.ToString()); sep ()
        w.Write(r.PrecursorMZ.ToString()); sep ()
        w.Write(r.TheoMass.ToString()); sep ()
        w.Write(r.AbsDeltaMass.ToString()); sep ()
        w.Write(r.PeptideLength.ToString()); sep ()
        w.Write(r.MissCleavages.ToString()); sep ()
        w.Write(r.SequestScore.ToString()); sep ()
        w.Write(r.SequestNormDeltaBestToRest.ToString()); sep ()
        w.Write(r.SequestNormDeltaNext.ToString()); sep ()
        w.Write(r.AndroScore.ToString()); sep ()
        w.Write(r.AndroNormDeltaBestToRest.ToString()); sep ()
        w.Write(r.AndroNormDeltaNext.ToString()); sep ()
        w.Write(r.XtandemScore.ToString()); sep ()
        w.Write(r.XtandemNormDeltaBestToRest.ToString()); sep ()
        w.Write(r.XtandemNormDeltaNext.ToString()); sep ()
        w.Write r.StringSequence; sep ()
        w.Write(r.IonMobility.ToString()); sep ()
        w.Write(r.Hyperscore.ToString()); sep ()
        w.Write(r.Expectscore.ToString()); sep ()
        w.Write(r.MatchedIons.ToString()); sep ()
        w.Write(r.TotalIons.ToString())
        w.WriteLine()

    /// What the index passes keep per spectrum and charge across the slices.
    type Attempt =
        {
            mutable TooFewPeaks : bool
            /// A slice failed on this attempt; it gets no rows and counts as failed.
            mutable Failed      : bool
            mutable Searched    : bool
            mutable Histogram   : int[]
            mutable Targets     : Candidate[]
            mutable Decoys      : Candidate[]
        }

    /// The best candidates of a label after one more slice: the kept ones and the new ones
    /// together in candidate order, cut to the reported number.
    let mergeTop (kept: Candidate[]) (fresh: Candidate[]) (limit: int) =
        if limit <= 0 then [||]
        elif kept.Length = 0 then fresh |> Array.sortWith (fun a b -> candidateOrder.Compare(a, b)) |> Array.truncate limit
        elif fresh.Length = 0 then kept |> Array.sortWith (fun a b -> candidateOrder.Compare(a, b)) |> Array.truncate limit
        else
            let all = Array.append kept fresh
            Array.Sort(all, candidateOrder)
            Array.sub all 0 (min all.Length limit)

    /// Combine histograms without materializing the candidate scores.
    let addHistogram (attempt: Attempt) (histogram: int[]) =
        if histogram.Length > attempt.Histogram.Length then
            let grown = Array.zeroCreate (max histogram.Length (attempt.Histogram.Length * 2))
            Array.blit attempt.Histogram 0 grown 0 attempt.Histogram.Length
            attempt.Histogram <- grown
        for bin = 0 to histogram.Length - 1 do
            attempt.Histogram.[bin] <- attempt.Histogram.[bin] + histogram.[bin]

    /// The histogram of an attempt cut to the last bin with a count plus one, which is the
    /// length the one pass model built.
    let private trimmedHistogram (attempt: Attempt) =
        let last = attempt.Histogram |> Array.tryFindIndexBack (fun c -> c > 0)
        match last with
        | Some b -> Array.sub attempt.Histogram 0 (b + 2)
        | None -> Array.zeroCreate 1

    /// One run between the passes. Every run is searched on its own thread, so nothing in here
    /// is shared with another run.
    type Run =
        {
            Path         : string
            Log          : string -> unit
            Reader       : IMzIODataReader
            Headers      : Ms2Header[]
            /// The first attempt of every spectrum; one attempt per spectrum and charge, laid out
            /// one after the other.
            AttemptStart : int[]
            Attempts     : Attempt[]
            /// The preprocessed spectrum of every attempt from the first slice on, null when the
            /// spectra are not kept. Released before the final pass.
            mutable Kept : ProcessedSpectrum[]
            Scratch      : Scratch
            /// The result file, opened when every run of the group has opened, null before.
            mutable Writer : StreamWriter
            Stopwatch    : Diagnostics.Stopwatch
            /// Time spent per phase of this run, see Profile.
            PhaseTicks   : int64[]
        }

    let private chargesOf (fallbackCharges: int[]) (header: Ms2Header) =
        if header.Charge > 0 then [| header.Charge |] else fallbackCharges

    /// Opens one run: the reader, the headers and the attempts. The run's stopwatch starts before
    /// the headers are read. The result file is not touched yet.
    let private openRun (processParams: PeptideSpectrumMatchingTIMsParams) (outputDir: string) (fallbackCharges: int[]) (keepSpectra: bool) (path: string) =
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension path)
        let log (msg: string) = logger.Trace msg
        log (sprintf "Input file: %s" path)
        log (sprintf "Output directory: %s" outputDir)
        log (sprintf "Parameters: %A" processParams)
        let outFilePath = Path.Combine(outputDir, Path.GetFileNameWithoutExtension path + ".psm")
        log (sprintf "Result file path: %s" outFilePath)
        log "Init connection to input data base."
        let reader = MzIO.Reader.getReader path
        try
            MzIO.Reader.openConnection reader
            let runId = MzIO.Reader.getDefaultRunID reader
            log (sprintf "Run ID: %s" runId)
            log "Starting peptide spectrum matching."
            let stopwatch = Diagnostics.Stopwatch.StartNew()
            let phaseTicks = Profile.create ()
            let headers = Profile.countInto phaseTicks (fun () -> readMs2Headers reader runId log)
            let attemptStart = headers |> Array.scan (fun start header -> start + (chargesOf fallbackCharges header).Length) 0
            let attempts = Array.init attemptStart.[headers.Length] (fun _ -> { TooFewPeaks = false; Failed = false; Searched = false; Histogram = Array.zeroCreate 64; Targets = [||]; Decoys = [||] })
            log (sprintf "%i MS2 spectra, %i spectrum and charge attempts." headers.Length attempts.Length)
            {
                Path = path
                Log = log
                Reader = reader
                Headers = headers
                AttemptStart = attemptStart
                Attempts = attempts
                Kept = if keepSpectra then Array.zeroCreate attempts.Length else null
                Scratch = Scratch()
                Writer = null
                Stopwatch = stopwatch
                PhaseTicks = phaseTicks
            }
        with _ ->
            reader.Dispose()
            reraise()

    /// Creates the result file of a run with its header line.
    let private openWriter (outputDir: string) (run: Run) =
        let writer = new StreamWriter(Path.Combine(outputDir, Path.GetFileNameWithoutExtension run.Path + ".psm"), false)
        try
            Reflection.FSharpType.GetRecordFields(typeof<Dto.PeptideSpectrumMatchingResult>)
            |> Array.map (fun field -> field.Name)
            |> String.concat "\t"
            |> writer.WriteLine
            run.Writer <- writer
        with _ ->
            writer.Dispose()
            reraise()

    /// Disposal of the reader is attempted even when the writer cannot flush.
    let private closeRun (run: Run) =
        try
            if not (isNull run.Writer) then run.Writer.Dispose()
        finally run.Reader.Dispose()

    let markReadFailure (run: Run) (i: int) (ex: exn) =
        for a = run.AttemptStart.[i] to run.AttemptStart.[i + 1] - 1 do
            run.Attempts.[a].Failed <- true
        run.Log (sprintf "spec with id: %s cannot be read: %A" run.Headers.[i].Id ex)

    /// The highest peptide mass any precursor window of the run reaches. Fragments above it
    /// cannot be matched by any spectrum of the run.
    let private runCeiling (settings: SearchSettings) (fallbackCharges: int[]) (run: Run) =
        run.Headers
        |> Array.fold (fun acc header ->
            chargesOf fallbackCharges header
            |> Array.fold (fun acc charge -> max acc (windowCeiling settings ((header.PrecursorMz - protonMass) * float charge))) acc) 0.

    /// The index pass of one slice over one run. The first slice reads the peaks of every
    /// spectrum. With kept spectra the later slices search those and read nothing.
    let searchSliceWith (read: IMzIODataReader -> Ms2Header[] -> int -> Result<Ms2Spectrum, exn>) (settings: SearchSettings) (table: PeptideTable) (fallbackCharges: int[]) (maxFragments: int)
                            (index: FragmentIndex) (rankLo: int) (rankHi: int) (slice: int) (run: Run) =
        Profile.countInto run.PhaseTicks (fun () ->
            let log = run.Log
            let searchAttempt (spectrumId: string) (precursorMz: float) (charge: int) (attempt: Attempt) (processed: ProcessedSpectrum) =
                try
                    match indexPass settings table index run.Scratch maxFragments precursorMz processed charge rankLo rankHi with
                    | SliceTooFewPeaks -> attempt.TooFewPeaks <- true
                    | SliceNoCandidates -> ()
                    | SliceFound (scores, targets, decoys) ->
                        attempt.Searched <- true
                        addHistogram attempt scores
                        attempt.Targets <- mergeTop attempt.Targets targets settings.ReportedHitsPerLabel
                        attempt.Decoys <- mergeTop attempt.Decoys decoys settings.ReportedHitsPerLabel
                with ex ->
                    attempt.Failed <- true
                    log (sprintf "spec with id: %s at charge %i fails in slice %i with: %A" spectrumId charge (slice + 1) ex)
            let progress (i: int) =
                if i % 50000 = 0 && i > 0 then
                    log (sprintf "slice %i: %i spectra searched, %.1f s." (slice + 1) i run.Stopwatch.Elapsed.TotalSeconds)
            if slice > 0 && not (isNull run.Kept) then
                run.Headers
                |> Array.iteri (fun i header ->
                    progress i
                    chargesOf fallbackCharges header
                    |> Array.iteri (fun c charge ->
                        let a = run.AttemptStart.[i] + c
                        if not run.Attempts.[a].Failed && not (isNull (box run.Kept.[a])) then
                            searchAttempt header.Id header.PrecursorMz charge run.Attempts.[a] run.Kept.[a]))
            else
                run.Headers |> Array.iteri (fun i header ->
                    progress i
                    match read run.Reader run.Headers i with
                    | Error ex -> markReadFailure run i ex
                    | Ok spectrum ->
                        chargesOf fallbackCharges header
                        |> Array.iteri (fun c charge ->
                            let a = run.AttemptStart.[i] + c
                            if not run.Attempts.[a].Failed then
                                try
                                    let t = Diagnostics.Stopwatch.GetTimestamp()
                                    let processed = preprocess settings spectrum.PrecursorMz charge spectrum.Mz spectrum.Intensity
                                    Profile.add 2 t
                                    if not (isNull run.Kept) then run.Kept.[a] <- processed
                                    searchAttempt spectrum.Id spectrum.PrecursorMz charge run.Attempts.[a] processed
                                with ex ->
                                    run.Attempts.[a].Failed <- true
                                    log (sprintf "spec with id: %s fails preprocessing at charge %i: %A" spectrum.Id charge ex)))
            log (sprintf "Slice %i searched, %.1f s." (slice + 1) run.Stopwatch.Elapsed.TotalSeconds))

    /// The final pass over one run: the expectation model from the histogram, the classic scores
    /// of the kept candidates, the rows.
    let finishRunWith (read: IMzIODataReader -> Ms2Header[] -> int -> Result<Ms2Spectrum, exn>) (settings: SearchSettings) (classic: ClassicScoring) (table: PeptideTable) (fallbackCharges: int[]) (run: Run) =
        Profile.countInto run.PhaseTicks (fun () ->
            let log = run.Log
            log "Final scoring of the kept candidates."
            let results = ResizeArray<SearchOutcome>()
            let withoutScoring (attempt: Attempt) =
                if attempt.Failed then Some Failed
                elif attempt.TooFewPeaks then Some TooFewPeaks
                elif not attempt.Searched then Some NoCandidates
                elif not (Array.exists (fun c -> c.ExactMatched >= settings.MinMatchedFragments) attempt.Targets ||
                          Array.exists (fun c -> c.ExactMatched >= settings.MinMatchedFragments) attempt.Decoys) then Some NoHits
                else None
            run.Headers |> Array.iteri (fun i header ->
                let charges = chargesOf fallbackCharges header
                let known = charges |> Array.mapi (fun c _ -> withoutScoring run.Attempts.[run.AttemptStart.[i] + c])
                if Array.forall Option.isSome known then
                    for outcome in known do results.Add outcome.Value
                else
                    match read run.Reader run.Headers i with
                    | Error ex ->
                        markReadFailure run i ex
                        for _ in charges do results.Add Failed
                    | Ok spectrum ->
                        charges |> Array.iteri (fun c charge ->
                            let attempt = run.Attempts.[run.AttemptStart.[i] + c]
                            let outcome =
                                match known.[c] with
                                | Some outcome -> outcome
                                | None ->
                                    try
                                        let outcome, rows = finishSpectrum settings classic table spectrum charge (trimmedHistogram attempt) attempt.Targets attempt.Decoys
                                        let t = Diagnostics.Stopwatch.GetTimestamp()
                                        rows |> Array.iter (writeRow run.Writer)
                                        Profile.add 9 t
                                        outcome
                                    with ex ->
                                        attempt.Failed <- true
                                        log (sprintf "spec with id: %s at charge %i fails with: %A" spectrum.Id charge ex)
                                        Failed
                            results.Add outcome))
            let outcomes = results |> Seq.countBy id |> Map.ofSeq
            let count outcome = outcomes |> Map.tryFind outcome |> Option.defaultValue 0
            let rows = outcomes |> Map.toSeq |> Seq.sumBy (fun (outcome, n) -> match outcome with Rows r -> r * n | _ -> 0)
            log (sprintf "Finished peptide spectrum matching: %i rows, %.1f s." rows run.Stopwatch.Elapsed.TotalSeconds)
            log (sprintf "Time by phase: %s." (Profile.summary run.PhaseTicks))
            log (sprintf "Spectrum and charge attempts without rows: %i with too few peaks, %i without candidate peptides, %i without matching candidates, %i failed." (count TooFewPeaks) (count NoCandidates) (count NoHits) (count Failed))
            run.Writer.Flush()
            log "Done.")

    /// Reject invalid settings before input or output resources are acquired.
    let validateParameters (p: PeptideSpectrumMatchingTIMsParams) (ceiling: float) (slices: int) =
        let nonnegative name x =
            if not (Double.IsFinite x) || x < 0. then invalidArg name "Must be finite and nonnegative."
        nonnegative "PrecursorTolerancePPM" p.PrecursorTolerancePPM
        nonnegative "FragmentTolerancePPM" p.FragmentTolerancePPM
        nonnegative "RemovePrecursorRange" p.RemovePrecursorRange
        nonnegative "MinimumPeakRatio" p.MinimumPeakRatio
        nonnegative "MaxFragmentMass" ceiling
        nonnegative "AndromedaTolerance" p.AndromedaParams.MatchingIonTolerancePPM
        if p.MinimumPeakRatio > 1. then invalidArg "MinimumPeakRatio" "Must not exceed one."
        if not (Double.IsFinite p.FragmentIndexBinWidth) || p.FragmentIndexBinWidth <= 0. then invalidArg "FragmentIndexBinWidth" "Must be finite and positive."
        if slices <= 0 then invalidArg "Slices" "Must be positive."
        if p.ReportedHitsPerLabel <= 0 || p.TopNPeaks <= 0 || p.MinimumPeaks <= 0 || p.MaxFragmentCharge <= 0 then
            invalidArg "Parameters" "Hit, peak and charge limits must be positive."
        if p.MinMatchedFragments < 0 || p.MinFragmentsModelling < 0 then invalidArg "Parameters" "Fragment thresholds must be nonnegative."
        if List.isEmpty p.IsotopeErrors then invalidArg "IsotopeErrors" "At least one isotope error is required."
        if List.isEmpty p.FallbackChargeStates then invalidArg "FallbackChargeStates" "At least one fallback charge is required."
        if p.FallbackChargeStates |> List.exists (fun c -> c <= 0) then invalidArg "FallbackChargeStates" "Charges must be positive."
        let pmin, pmax = p.AndromedaParams.PMinPMax
        if pmin <= 0 || pmax < pmin || pmax > 100 then invalidArg "Andromeda" "Require 1 <= PMin <= PMax <= 100."

    let validateOutputPaths (outputDir: string) (paths: string[]) =
        let seen = HashSet<string>(if OperatingSystem.IsWindows() then StringComparer.OrdinalIgnoreCase else StringComparer.Ordinal)
        for path in paths do
            let output = Path.GetFullPath(Path.Combine(outputDir, Path.GetFileNameWithoutExtension path + ".psm"))
            if not (seen.Add output) then invalidArg "paths" (sprintf "Inputs would overwrite the same output: %s" output)

    /// Searches the given runs together and writes <run>.psm for each of them into the output
    /// directory. The fragment index is built in slices, every slice once, and the runs are
    /// searched against it side by side, each on one thread. The memory of a call is one slice
    /// plus what the runs keep between the slices. The ceiling of the index is the highest
    /// peptide mass any precursor window of the runs reaches unless one is given, in which case
    /// a value below that bound is warned about.
    let private scoreRunsWithPeakReader read (processParams: PeptideSpectrumMatchingTIMsParams) (outputDir: string) (log: string -> unit) (maxFragmentMass: float) (slices: int) (keepSpectra: bool)
                  (sdbParams: SearchDB.SearchDbParams) (table: PeptideTable) (paths: string[]) =
        validateParameters processParams maxFragmentMass slices
        validateOutputPaths outputDir paths
        let settings = toSearchSettings processParams
        let classic = toClassicScoring processParams sdbParams table
        let fallbackCharges = Array.ofList processParams.FallbackChargeStates
        let maxFragments = 2 * table.MaxLength
        let workers = max 1 paths.Length
        // every run that opened is closed again when anything in the group fails
        let opened : Run[] = Array.zeroCreate paths.Length
        try
            Threading.Tasks.Parallel.For(0, paths.Length, fun i -> opened.[i] <- openRun processParams outputDir fallbackCharges keepSpectra paths.[i]) |> ignore
            let runs = opened
            let reach = runs |> Array.map (runCeiling settings fallbackCharges) |> Array.fold max 0.
            let ceiling =
                let bound =
                    if maxFragmentMass > 0. then
                        if maxFragmentMass < reach then
                            log (sprintf "The fragment mass ceiling %.1f Da is below the reach of the precursor windows, %.1f Da. Fragments between the two cannot be matched." maxFragmentMass reach)
                        maxFragmentMass
                    else reach + 1.
                min bound (table.Mass.[table.Count - 1] + 1.)
            if not (Double.IsFinite ceiling) || ceiling <= 0. then failwithf "The fragment mass ceiling %g Da is not a mass." ceiling
            log (sprintf "Fragment mass ceiling %.1f Da, the precursor windows reach %.1f Da." ceiling reach)
            runs |> Array.iter (fun run -> run.Log (sprintf "Fragment mass ceiling %.1f Da, the precursor windows of the runs searched together reach %.1f Da." ceiling reach))
            if ceiling / processParams.FragmentIndexBinWidth > float (Int32.MaxValue - 2) then
                invalidArg "FragmentIndexBinWidth" "Too many fragment bins."
            runs |> Array.iter (openWriter outputDir)
            let precursorMasses =
                seq {
                    for run in runs do
                        for header in run.Headers do
                            for charge in chargesOf fallbackCharges header do
                                yield (header.PrecursorMz - protonMass) * float charge
                }
            let eligible = reachableRanks settings table precursorMasses
            log (sprintf "Indexing %i of %i peptide ranks reached by the input precursor windows." (eligible |> Array.sumBy (fun x -> if x then 1 else 0)) table.Count)
            let boundaries = sliceBoundariesFor eligible table slices processParams.FragmentIndexBinWidth ceiling
            for slice = 0 to slices - 1 do
                let rankLo, rankHi = boundaries.[slice], boundaries.[slice + 1]
                log (sprintf "Slice %i of %i: ranks %i to %i." (slice + 1) slices rankLo rankHi)
                let index = buildFragmentIndexFor eligible table processParams.FragmentIndexBinWidth ceiling workers rankLo rankHi log
                runs |> Array.Parallel.iter (searchSliceWith read settings table fallbackCharges maxFragments index rankLo rankHi slice)
                log (sprintf "Slice %i of %i searched." (slice + 1) slices)
                GC.Collect()
            runs |> Array.iter (fun run -> run.Kept <- null)
            runs |> Array.Parallel.iter (finishRunWith read settings classic table fallbackCharges)
        finally
            for run in opened do
                if not (isNull (box run)) then
                    try closeRun run
                    with ex -> log (sprintf "Failed to close run %s: %A" run.Path ex)

    /// Optional additional centroid merging, independent of the matching tolerance.
    let scoreRunsWithMobilityTolerance mobilityMergePpm processParams outputDir log maxFragmentMass slices keepSpectra sdbParams table paths =
        validateMobilityMergePpm mobilityMergePpm
        log (sprintf "Summing mobility peaks within each stored spectrum; additional m/z merge tolerance %g ppm (0 = exact m/z)." mobilityMergePpm)
        scoreRunsWithPeakReader (readMs2SpectrumWithTolerance mobilityMergePpm) processParams outputDir log maxFragmentMass slices keepSpectra sdbParams table paths

    /// Mobility summation cannot be disabled. Without an explicit override, use the
    /// fragment matching tolerance of this run as the m/z cluster width.
    let scoreRuns processParams outputDir log maxFragmentMass slices keepSpectra sdbParams table paths =
        scoreRunsWithMobilityTolerance processParams.FragmentTolerancePPM processParams outputDir log maxFragmentMass slices keepSpectra sdbParams table paths
