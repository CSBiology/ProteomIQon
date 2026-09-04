namespace ProteomIQon

open System
open System.IO
open System.Collections.Generic
open System.Globalization
open BioFSharp
open BioFSharp.Mz
open MzIO
open MzIO.IO
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

    let toClassicScoring (processParams: PeptideSpectrumMatchingTIMsParams) (sdbParams: SearchDB.SearchDbParams) : ClassicScoring =
        {
            CalcIonSeries         = Fragmentation.Series.fragmentMasses processParams.nTerminalSeries processParams.cTerminalSeries sdbParams.MassFunction
            ParseSequence         = SearchDB.initOfModAminoAcidString sdbParams.IsotopicMod (sdbParams.FixedMods @ sdbParams.VariableMods)
            AndromedaPMinPMax     = processParams.AndromedaParams.PMinPMax
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

    /// The few values the search needs from a spectrum description, so that the descriptions
    /// themselves are not kept for the whole run.
    type private Ms2Header =
        {
            Id          : string
            ScanTime    : float
            PrecursorMz : float
            Charge      : int
            IonMobility : float
        }

    /// MS2 spectra of the run in scan time order. ScanNr is the position in that order. The
    /// peaks of a spectrum are read when the sequence reaches it. A spectrum whose precursor or
    /// peaks cannot be read is logged and skipped.
    let readMs2Spectra (reader: IMzIODataReader) (runId: string) (log: string -> unit) =
        reader.ReadMassSpectra runId
        |> Seq.filter (fun ms -> MassSpectrum.getMsLevel ms = 2)
        |> Seq.choose (fun ms ->
            try
                Some
                    {
                        Id = ms.ID
                        ScanTime = MassSpectrum.getScanTime ms
                        PrecursorMz = MassSpectrum.getPrecursorMZ ms
                        Charge = precursorCharge ms
                        IonMobility = ionMobility ms
                    }
            with ex ->
                log (sprintf "spec with id: %s cannot be read: %A" ms.ID ex)
                None)
        |> Seq.sortBy (fun header -> header.ScanTime)
        |> Seq.mapi (fun i header ->
            try
                let mz, intensity = MzIO.Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(header.Id).Peaks)
                Some
                    {
                        Id = header.Id
                        ScanNr = i
                        ScanTime = header.ScanTime
                        PrecursorMz = header.PrecursorMz
                        Charge = header.Charge
                        IonMobility = header.IonMobility
                        Mz = mz
                        Intensity = intensity
                    }
            with ex ->
                log (sprintf "spec with id: %s cannot be read: %A" header.Id ex)
                None)
        |> Seq.choose id

    /// Scan range as PeptideSpectrumMatching defines it: the spectrum limits rounded outwards to tens.
    let private scanRangeOf (mz: float[]) =
        let floorToClosest10 x = Math.Floor(x / 10.) * 10.
        let ceilToClosest10 x = Math.Ceiling(x / 10.) * 10.
        Math.Max(0., Math.Round(Array.min mz, 0) |> floorToClosest10), Math.Round(Array.max mz, 0) |> ceilToClosest10

    /// Scores the given peptides with the BioFSharp.Mz functions that PeptideSpectrumMatching
    /// uses. Every peptide yields a target and a reversed decoy result, keyed by
    /// (ModSequenceID, GlobalMod, isTarget).
    let classicScores (classic: ClassicScoring) (table: PeptideTable) (spectrum: Ms2Spectrum) (charge: int) (ranks: int[]) =
        let theoSpecs =
            ranks
            |> List.ofArray
            |> List.map (fun rank ->
                let sequence = table.SequenceString rank
                let globalMod = int table.GlobalMod.[rank]
                let mass = table.Mass.[rank]
                let lookUp =
                    SearchDB.createLookUpResult table.ModSequenceId.[rank] table.PepSequenceId.[rank] mass
                        (int64 (Math.Round(mass * 1000000.))) sequence (classic.ParseSequence globalMod sequence) globalMod
                lookUp, classic.CalcIonSeries lookUp.BioSequence)
        let recSpec = PeakArray.zip spectrum.Mz spectrum.Intensity
        let scanRange = scanRangeOf spectrum.Mz
        let sequest =
            SequestLike.getTheoSpecs scanRange charge theoSpecs
            |> fun t -> SequestLike.calcSequestScore scanRange recSpec spectrum.ScanTime charge spectrum.PrecursorMz t spectrum.Id
        let andro, xtandem =
            XScoring.getTheoSpecs scanRange charge theoSpecs
            |> fun t -> XScoring.calcAndromedaAndXTandemScore classic.AndromedaPMinPMax scanRange classic.AndromedaTolerancePpm recSpec spectrum.ScanTime charge spectrum.PrecursorMz t spectrum.Id
        let toMap (results: SearchEngineResult.SearchEngineResult<float> list) =
            results |> List.map (fun r -> (r.ModSequenceID, r.GlobalMod, r.IsTarget), r) |> dict
        toMap sequest, toMap andro, toMap xtandem

    /// Searches one spectrum at one charge state and returns the result rows, best hyperscore
    /// first within the targets and within the decoys.
    let searchSpectrum (settings: SearchSettings) (classic: ClassicScoring) (table: PeptideTable) (index: FragmentIndex)
                       (scratch: Scratch) (maxFragments: int) (spectrum: Ms2Spectrum) (charge: int) =
        let processed = preprocess settings spectrum.PrecursorMz charge spectrum.Mz spectrum.Intensity
        if processed.Mz.Length < settings.MinimumPeaks then TooFewPeaks, [||]
        else
            let precursorMass = (spectrum.PrecursorMz - protonMass) * float charge
            let windows = precursorWindows settings table precursorMass
            if windows.Length = 0 then NoCandidates, [||]
            else
                scratch.Reset (slotCount windows) maxFragments
                scatter settings index scratch windows charge processed
                let scores, targets, decoys = collect settings scratch windows
                if targets.Length = 0 && decoys.Length = 0 then NoHits, [||]
                else
                    let expectation = expectationModel scores
                    let ranks = Array.append targets decoys |> Array.map (fun c -> c.Rank) |> Array.distinct
                    let sequestMap, androMap, xtandemMap = classicScores classic table spectrum charge ranks
                    let score (m: IDictionary<_, SearchEngineResult.SearchEngineResult<float>>) key =
                        match m.TryGetValue key with
                        | true, r -> r.Score, r.NormDeltaBestToRest, r.NormDeltaNext
                        | _ -> 0., 0., 0.
                    // The exact rescoring decides the reported hyperscore, the rank and the
                    // MinMatchedFragments filter. The expectation value belongs to the index pass
                    // score, the quantity the survival model was built from.
                    let makeRows (candidates: Candidate[]) (label: int) =
                        candidates
                        |> Array.map (fun c ->
                            let nb, ny, sumB, sumY, total = rescore settings table scratch charge processed c.Rank c.IsDecoy
                            c, hyperscore nb ny sumB sumY, nb + ny, total)
                        |> Array.filter (fun (_, _, matched, _) -> matched >= settings.MinMatchedFragments)
                        |> Array.sortByDescending (fun (_, hs, _, _) -> hs)
                        |> Array.mapi (fun i (c, hs, matched, total) ->
                            let key = (table.ModSequenceId.[c.Rank], int table.GlobalMod.[c.Rank], not c.IsDecoy)
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
                                    GlobalMod = int table.GlobalMod.[c.Rank]
                                    PepSequenceID = table.PepSequenceId.[c.Rank]
                                    ModSequenceID = table.ModSequenceId.[c.Rank]
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
                                    Hyperscore = hs
                                    Expectscore = expectation c.Hyperscore
                                    MatchedIons = matched
                                    TotalIons = total
                                }
                            row)
                    let rows = Array.append (makeRows targets 1) (makeRows decoys -1)
                    if rows.Length = 0 then NoHits, [||] else Rows rows.Length, rows

    /// Loads the peptide database and builds the fragment index once for all runs, with as
    /// many workers as runs are processed at the same time. The index holds monoisotopic
    /// fragment masses, so a database with average masses is refused.
    let prepareIndex (processParams: PeptideSpectrumMatchingTIMsParams) (cn: System.Data.SQLite.SQLiteConnection) (workers: int) (log: string -> unit) =
        let sdbParams = SearchDB.getSDBParamsByCn cn
        match sdbParams.MassMode with
        | SearchDB.MassMode.Monoisotopic -> ()
        | SearchDB.MassMode.Average -> failwith "The peptide data base uses average masses. PeptideSpectrumMatchingTIMs matches monoisotopic fragment masses and needs a monoisotopic data base."
        let table = loadPeptideTable cn sdbParams log
        let deviation = verifyMasses table 0.001
        log (sprintf "Residue masses reproduce the database masses, worst deviation %g Da." deviation)
        let index = buildFragmentIndex table processParams.FragmentIndexBinWidth (table.Mass.[table.Count - 1] + 10.) workers log
        sdbParams, table, index

    /// Searches one mzlite file and writes <run>.psm into the output directory.
    let scoreSpectra (processParams: PeptideSpectrumMatchingTIMsParams) (outputDir: string)
                     (sdbParams: SearchDB.SearchDbParams) (table: PeptideTable) (index: FragmentIndex) (instrumentOutput: string) =
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension instrumentOutput)
        let log (msg: string) = logger.Trace msg
        log (sprintf "Input file: %s" instrumentOutput)
        log (sprintf "Output directory: %s" outputDir)
        log (sprintf "Parameters: %A" processParams)
        let outFilePath = Path.Combine(outputDir, Path.GetFileNameWithoutExtension instrumentOutput + ".psm")
        log (sprintf "Result file path: %s" outFilePath)
        let settings = toSearchSettings processParams
        let classic = toClassicScoring processParams sdbParams
        let fallbackCharges = Array.ofList processParams.FallbackChargeStates
        let maxFragments = 2 * table.MaxLength
        let scratch = Scratch()
        log "Init connection to input data base."
        let inReader = MzIO.Reader.getReader instrumentOutput
        MzIO.Reader.openConnection inReader
        let inRunID = MzIO.Reader.getDefaultRunID inReader
        log (sprintf "Run ID: %s" inRunID)
        let inTr = inReader.BeginTransaction()
        use resultWriter = new StreamWriter(outFilePath, false)
        Reflection.FSharpType.GetRecordFields(typeof<Dto.PeptideSpectrumMatchingResult>)
        |> Array.map (fun field -> field.Name)
        |> String.concat "\t"
        |> resultWriter.WriteLine
        log "Starting peptide spectrum matching."
        let stopwatch = Diagnostics.Stopwatch.StartNew()
        let outcomes =
            readMs2Spectra inReader inRunID log
            |> Seq.collect (fun spectrum ->
                if spectrum.ScanNr % 50000 = 0 && spectrum.ScanNr > 0 then
                    log (sprintf "%i spectra searched, %.1f s." spectrum.ScanNr stopwatch.Elapsed.TotalSeconds)
                let charges = if spectrum.Charge > 0 then [| spectrum.Charge |] else fallbackCharges
                charges
                |> Array.map (fun charge ->
                    try
                        let outcome, rows = searchSpectrum settings classic table index scratch maxFragments spectrum charge
                        rows |> SeqIO'.csv "\t" false false |> Seq.iter resultWriter.WriteLine
                        outcome
                    with ex ->
                        log (sprintf "spec with id: %s at charge %i fails with: %A" spectrum.Id charge ex)
                        Failed))
            |> Seq.countBy id
            |> Map.ofSeq
        let count outcome = outcomes |> Map.tryFind outcome |> Option.defaultValue 0
        let rows = outcomes |> Map.toSeq |> Seq.sumBy (fun (outcome, n) -> match outcome with Rows r -> r * n | _ -> 0)
        log (sprintf "Finished peptide spectrum matching: %i rows, %.1f s." rows stopwatch.Elapsed.TotalSeconds)
        log (sprintf "Spectrum and charge attempts without rows: %i with too few peaks, %i without candidate peptides, %i without matching candidates, %i failed." (count TooFewPeaks) (count NoCandidates) (count NoHits) (count Failed))
        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()
        log "Done."
