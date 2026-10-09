namespace ProteomIQon

open System.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute
open Deedle
open System.Data
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open System.Data.SQLite
open ProteomIQon
open MzIO.Processing
open ProteomIQon.Dto
open ProteomIQon.Domain


module MsFraggerToPSM =

    let readFragegrPsms (path: string) =
        Frame.ReadCsv(path, true, separators = "\t")
        |> Frame.mapRows (fun k row ->
            let assignedModifications =
                row.GetAs<string>("Assigned Modifications").Split(',')
                |> Array.sort
                |> Array.filter (fun s -> s <> "")
            let stringSequence = 
                row.GetAs<string>("Peptide").ToCharArray()
                |> Array.map string
            let modStringSequence =
                assignedModifications
                |> Array.fold (fun (acc: string []) (modString: string) ->
                    let modInfo = modString.Split('(') |> Array.head
                    let modPosition, modification =
                        if modInfo.StartsWith("N-term") then 0, "nTermAcetylation"
                        else
                            if modInfo.Length = 2 then
                                int modInfo.[0..0] - 1, modInfo.[1..]
                            else
                                int modInfo.[0..1] - 1, modInfo.[2..]
                    match modification with
                    | "M" ->
                        acc.[modPosition] <- "[ox]M"
                        acc
                    | "C" ->
                        acc.[modPosition] <- "[ca]C"
                        acc
                    | "nTermAcetylation" ->
                        acc.[0] <- "[ac]" + acc.[0]
                        acc
                    | m ->
                        printfn "Modification %s is not yet implemented" m
                        Array.append acc [|"UnknownMod"|]
                ) stringSequence
                |> String.concat ""
            {|
                ScanTime = row.GetAs<float>("Retention") / 60.
                Charge = row.GetAs<int>("Charge")
                PrecursorMZ = row.GetAs<float>("Calibrated Observed M/Z")
                TheoMass = row.GetAs<float>("Calculated Peptide Mass")
                AbsDeltaMass = row.GetAs<float>("Delta Mass") |> abs
                // NaN when the column or the value is missing
                IonMobility =
                    let ionMobility = row.TryGetAs<float>("Ion Mobility")
                    if ionMobility.HasValue then ionMobility.Value else nan
                PeptideLength = row.GetAs<int>("Peptide Length")
                MissCleavages = row.GetAs<int>("Number of Missed Cleavages")
                Hyperscore = row.GetAs<float>("Hyperscore")
                Probability = row.GetAs<float>("Probability")
                StringSequence = modStringSequence
                ProteinNames = row.GetAs<string>("Protein")
            |}
        )
        |> Series.values
        |> Array.ofSeq

    /// q-values from posterior error probabilities: the q-value of a PSM is the mean PEP of all
    /// PSMs whose PEP is not larger, the expected false discovery rate when accepting down to it.
    /// A PSM without PEP (NaN) gets the q-value NaN and does not count for the others.
    let qValuesOfPeps (peps: float []) =
        let order =
            Array.init peps.Length id
            |> Array.filter (fun i -> not (System.Double.IsNaN peps.[i]))
            |> Array.sortBy (fun i -> peps.[i])
        let qValues = Array.create peps.Length nan
        let mutable sum = 0.
        let mutable k = 0
        while k < order.Length do
            // PSMs with the same PEP share the q-value of the last of them
            let mutable last = k
            while last + 1 < order.Length && peps.[order.[last + 1]] = peps.[order.[k]] do
                last <- last + 1
            for j in k .. last do
                sum <- sum + peps.[order.[j]]
            let q = sum / float (last + 1)
            for j in k .. last do
                qValues.[order.[j]] <- q
            k <- last + 1
        qValues

    /// Prepares statement to select a ModSequence entry by Sequence
    let prepareSelectModsequenceBySequence (cn:SQLiteConnection) (tr) =
        let querystring = "SELECT * FROM ModSequence WHERE Sequence=@sequence AND GlobalMod=@globalMod"
        let cmd = new SQLiteCommand(querystring, cn, tr) 
        cmd.Parameters.Add("@sequence", System.Data.DbType.String) |> ignore
        cmd.Parameters.Add("@globalMod", System.Data.DbType.Int32) |> ignore
        (fun (sequence:string) (globalMod:int) ->        
            cmd.Parameters.["@sequence"].Value <- sequence
            cmd.Parameters.["@globalMod"].Value <- globalMod
            use reader = cmd.ExecuteReader()            
            match reader.Read() with
            | true -> 
                {|
                    ID = reader.GetInt32(0);
                    PepSequenceID = reader.GetInt32(1);
                    RealMass = reader.GetDouble(2);
                    RoundedMass = reader.GetInt64(3);
                    Sequence = reader.GetString(4);
                    GlobalMod = reader.GetInt32(5)
                |}
                |> Some
            | false -> None
        )

    let initModSeqLookup (cn:SQLiteConnection) tr =
        let selectModsequenceBySequence = prepareSelectModsequenceBySequence cn tr
        selectModsequenceBySequence

    /// An MS2 spectrum of the mzlite. ScanNr is its position among the MS2 spectra in scan time order, equal scan
    /// times in file order, as in PeptideSpectrumMatchingTIMs.
    type Ms2Spectrum =
        {
            ID          : string
            ScanNr      : int
            ScanTime    : float
            PrecursorMZ : float
            Charge      : int option
            IonMobility : float
        }

    /// The MS2 spectra of the run from one pass over the spectrum descriptions. Like PeptideSpectrumMatchingTIMs,
    /// spectra without ID, scan time or precursor m/z are left out and do not count for ScanNr.
    let readMs2Spectra (reader: MzIO.MzSQL.MzSQL) (runID: string) (logger: NLog.Logger) =
        let spectra = ResizeArray<int*SpectrumDescription.SpectrumHeader>()
        do
            use cmd = new SQLiteCommand("SELECT Description FROM Spectrum WHERE RunID = @runID ORDER BY rowid", reader.Connection)
            cmd.Parameters.AddWithValue("@runID", runID) |> ignore
            use rows = cmd.ExecuteReader()
            let mutable position = 0
            while rows.Read() do
                let h = SpectrumDescription.tryParse (rows.GetString 0)
                if not (isNull h.Error) then
                    logger.Warn (sprintf "spectrum at position %i cannot be read: %s" position h.Error)
                elif h.MsLevel = Some 2 &&
                     (System.String.IsNullOrWhiteSpace h.ID ||
                      not (h.ScanTime |> Option.exists (fun x -> System.Double.IsFinite x && x >= 0.)) ||
                      not (h.PrecursorMz |> Option.exists (fun x -> System.Double.IsFinite x && x > BioFSharp.Mass.Table.PMassInU)) ||
                      (h.ChargeState |> Option.exists (fun x -> x < 0))) then
                    logger.Warn (sprintf "spectrum at position %i has missing or invalid metadata; skipped." position)
                elif h.MsLevel = Some 2 then
                    spectra.Add (position, h)
                position <- position + 1
        spectra.ToArray()
        |> Array.sortBy (fun (position, h) -> h.ScanTime.Value, position)
        |> Array.mapi (fun i (_, h) ->
            {
                ID          = h.ID
                ScanNr      = i
                ScanTime    = h.ScanTime.Value
                PrecursorMZ = h.PrecursorMz.Value
                // charge 0 means unknown
                Charge      = h.ChargeState |> Option.filter (fun x -> x > 0)
                IonMobility = h.IonMobility |> Option.defaultValue nan
            }
        )

    /// Ion mobility (1/K0) within which a spectrum is a candidate for a PSM of the same ion mobility. FragPipe writes the
    /// ion mobility of the precursor with 4 decimals, and the spectra of that precursor carry the same value.
    let ionMobilityTolerance = 0.001

    /// Precursor m/z deviation (ppm) of the matched spectrum beyond which a warning is logged.
    let precursorMzTolerancePpm = 20.

    /// The spectrum a PSM of the tsv was identified from. FragPipe names it by its own scan numbers, which do
    /// not map to the mzlite, so it is found by charge, ion mobility, scan time and precursor m/z: candidates
    /// are the spectra of the PSM charge or without charge. Among the candidates within ionMobilityTolerance of
    /// the PSM ion mobility the one with the smallest |dScanTime| + |dPrecursorMZ| is taken, among all candidates
    /// when none is that close or the PSM has no ion mobility, the first in scan time order on a tie. None when
    /// there is no candidate. The spectra are sorted by scan time, so the search runs outward from the scan time
    /// of the PSM and stops where the scan time difference alone exceeds the best distance.
    let findSpectrum (scanTime: float, precursorMz: float, charge: int, ionMobility: float) (spectra: Ms2Spectrum []) =
        let distance (s: Ms2Spectrum) = abs (s.ScanTime - scanTime) + abs (s.PrecursorMZ - precursorMz)
        let mutable bestWithin = -1
        let mutable bestOther = -1
        let better i best =
            best < 0 || (let d, b = distance spectra.[i], distance spectra.[best] in d < b || (d = b && i < best))
        let consider i =
            let s = spectra.[i]
            if s.Charge |> Option.forall ((=) charge) then
                if abs (s.IonMobility - ionMobility) <= ionMobilityTolerance then
                    if better i bestWithin then bestWithin <- i
                elif better i bestOther then bestOther <- i
        // a candidate within the ion mobility tolerance wins at any distance, so without one the bound is open
        let bound () =
            if bestWithin >= 0 then distance spectra.[bestWithin]
            elif System.Double.IsNaN ionMobility && bestOther >= 0 then distance spectra.[bestOther]
            else infinity
        // first spectrum not before scanTime
        let mutable lo = 0
        let mutable hi = spectra.Length
        while lo < hi do
            let mid = (lo + hi) / 2
            if spectra.[mid].ScanTime < scanTime then lo <- mid + 1 else hi <- mid
        let mutable left = lo - 1
        let mutable right = lo
        let mutable searching = true
        while searching && (left >= 0 || right < spectra.Length) do
            let i =
                if left < 0 then right
                elif right >= spectra.Length then left
                elif scanTime - spectra.[left].ScanTime <= spectra.[right].ScanTime - scanTime then left
                else right
            if abs (spectra.[i].ScanTime - scanTime) > bound () then
                searching <- false
            else
                consider i
                if i = left then left <- left - 1 else right <- right + 1
        if bestWithin >= 0 then Some spectra.[bestWithin]
        elif bestOther >= 0 then Some spectra.[bestOther]
        else None

    let convertToPSM (inputFileMsFragger: string) (inputFileMzLite: string) (outputDir: string) (cn: SQLite.SQLiteConnection) =

        let logger = Logging.createLogger "MSFraggerToPSM"

        let psms =
            let rows = readFragegrPsms inputFileMsFragger
            // FragPipe reports the probability that a PSM is correct, so PEP = 1 - Probability. The q-values
            // are computed over all rows of the tsv, the PSMs FragPipe accepted.
            let peps = rows |> Array.map (fun x -> 1. - x.Probability)
            let qValues = qValuesOfPeps peps
            Array.zip3 rows peps qValues
            |> Array.filter (fun (x, _, _) -> x.StringSequence.Contains "UnknownMod" |> not)

        let outputPath = Path.Combine(outputDir, Path.GetFileNameWithoutExtension(inputFileMzLite) + ".qpsm")

        let memoryDB = SearchDB.copyDBIntoMemory cn
        let pepDBTr = memoryDB.BeginTransaction()

        let modSeqLookup = initModSeqLookup memoryDB pepDBTr

        let idPSMs =
            psms
            |> Array.choose (fun (v, pep, qValue) ->
                modSeqLookup v.StringSequence 0
                |> Option.map (fun res -> v, pep, qValue, res.ID, res.PepSequenceID)
            )

        let inReader = Core.MzIO.Reader.getReader inputFileMzLite
        Core.MzIO.Reader.openConnection inReader
        let inTr = inReader.BeginTransaction()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader


        let spectra =
            match inReader with
            | :? MzIO.MzSQL.MzSQL as sql -> readMs2Spectra sql inRunID logger
            | _ -> failwith "MsFraggerToPSM reads mzlite files."

        let mutable withoutIonMobility = 0
        let psmsWithSpecID : PSMStatisticsResult [] =
            idPSMs
            |> Array.choose (fun (psm, pep, qValue, modSequenceID, pepSequenceID) ->
                match findSpectrum (psm.ScanTime, psm.PrecursorMZ, psm.Charge, psm.IonMobility) spectra with
                | None ->
                    logger.Warn (sprintf "%s at scan time %f, m/z %f, charge %i: no MS2 spectrum of this charge; dropped." psm.StringSequence psm.ScanTime psm.PrecursorMZ psm.Charge)
                    None
                | Some spectrum ->
                    let mzDeviationPpm = abs (spectrum.PrecursorMZ - psm.PrecursorMZ) / psm.PrecursorMZ * 1e6
                    if System.Double.IsNaN spectrum.IonMobility then withoutIonMobility <- withoutIonMobility + 1
                    if abs (spectrum.IonMobility - psm.IonMobility) > ionMobilityTolerance || mzDeviationPpm > precursorMzTolerancePpm then
                        logger.Warn (
                            sprintf "%s at scan time %f, m/z %f, charge %i, ion mobility %f matched to %s at scan time %f, m/z %f (%.1f ppm), ion mobility %f."
                                psm.StringSequence psm.ScanTime psm.PrecursorMZ psm.Charge psm.IonMobility spectrum.ID spectrum.ScanTime spectrum.PrecursorMZ mzDeviationPpm spectrum.IonMobility)
                    {
                        PSMId                       = spectrum.ID
                        GlobalMod                   = 0
                        PepSequenceID               = pepSequenceID
                        ModSequenceID               = modSequenceID
                        // psm.tsv lists only target PSMs
                        Label                       = 1
                        ScanNr                      = spectrum.ScanNr
                        ScanTime                    = psm.ScanTime
                        Charge                      = psm.Charge
                        PrecursorMZ                 = psm.PrecursorMZ
                        TheoMass                    = psm.TheoMass
                        AbsDeltaMass                = psm.AbsDeltaMass
                        PeptideLength               = psm.PeptideLength
                        MissCleavages               = psm.MissCleavages
                        SequestScore                = nan
                        SequestNormDeltaBestToRest  = nan
                        SequestNormDeltaNext        = nan
                        AndroScore                  = nan
                        AndroNormDeltaBestToRest    = nan
                        AndroNormDeltaNext          = nan
                        // the MSFragger hyperscore is the X!Tandem hyperscore
                        XtandemScore                = psm.Hyperscore
                        XtandemNormDeltaBestToRest  = nan
                        XtandemNormDeltaNext        = nan
                        ModelScore                  = psm.Probability
                        QValue                      = qValue
                        PEPValue                    = pep
                        StringSequence              = psm.StringSequence
                        ProteinNames                = psm.ProteinNames
                    }
                    |> Some
            )

        if withoutIonMobility > 0 then
            logger.Warn (sprintf "%i of %i PSMs were matched to a spectrum without ion mobility." withoutIonMobility psmsWithSpecID.Length)

        psmsWithSpecID
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outputPath)
