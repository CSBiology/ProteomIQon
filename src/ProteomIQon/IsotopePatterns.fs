namespace ProteomIQon

open System
open System.Data
open System.Data.SQLite
open System.Globalization
open System.Reflection
open BioFSharp
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB

/// The theoretical isotope patterns of the precursor ions, computed in one place and stored in the peptide data base by PSMStatistics.
module IsotopePatterns =

    /// Resolution in Da of the isotope patterns.
    let resolution = 0.01

    /// Peaks of the isotope patterns with a lower probability are dropped.
    let minProbability = 0.001

    ///Predicts an isotopic distribution of the given peptide at the given charge, normalized by the sum of probabilities, using the MIDAs algorithm
    let compute (charge:int) (seq:AminoAcids.AminoAcid list) =
        seq
        |> BioList.toFormula
        |> Formula.add Formula.Table.H2O
        |> IsotopicDistribution.MIDA.ofFormula IsotopicDistribution.MIDA.normalizeByProbSum resolution minProbability charge
        |> Array.ofList

    let private versionOf (assembly: Assembly) =
        match assembly.GetCustomAttribute<AssemblyInformationalVersionAttribute>() with
        | null -> assembly.GetName().Version.ToString()
        | attribute -> attribute.InformationalVersion

    /// Everything the patterns depend on besides the peptide and the charge: the versions of BioFSharp and BioFSharp.Mz (the latter builds the sequence of a
    /// ModSequence row), the function with its parameters, the formula and the encoding. Patterns are only used if the data base holds them under this text.
    let settings =
        sprintf "BioFSharp %s; BioFSharp.Mz %s; IsotopicDistribution.MIDA.ofFormula normalizeByProbSum resolution=%s minProbability=%s; formula=BioList.toFormula+H2O; peaks=float64 little-endian (m/z, probability) pairs"
            (versionOf typeof<AminoAcids.AminoAcid>.Assembly)
            (versionOf typeof<SearchDbParams>.Assembly)
            (resolution.ToString("R", CultureInfo.InvariantCulture))
            (minProbability.ToString("R", CultureInfo.InvariantCulture))

    /// The peaks as float64 little-endian (m/z, probability) pairs, 16 bytes per peak, with the exact bits of the values.
    let encode (peaks: (float*float)[]) : byte[] =
        let bytes = Array.zeroCreate<byte> (peaks.Length * 16)
        let write (offset: int) (value: float) =
            let bits = BitConverter.DoubleToInt64Bits value
            for i = 0 to 7 do
                bytes.[offset + i] <- byte (bits >>> (8 * i))
        peaks
        |> Array.iteri (fun i (mz, probability) ->
            write (i * 16) mz
            write (i * 16 + 8) probability)
        bytes

    /// The peaks of encode, None if the length of the bytes is not a multiple of 16.
    let decode (bytes: byte[]) : (float*float)[] option =
        if isNull bytes || bytes.Length % 16 <> 0 then None
        else
            let read (offset: int) =
                let mutable bits = 0L
                for i = 7 downto 0 do
                    bits <- (bits <<< 8) ||| int64 bytes.[offset + i]
                BitConverter.Int64BitsToDouble bits
            Some (Array.init (bytes.Length / 16) (fun i -> read (i * 16), read (i * 16 + 8)))

    /// Returns the isotope pattern of a ModSequence row at the given charge from the table IsotopePattern of the data base of the connection and computes it if the
    /// table, the settings row of this version or the row of the pattern is missing, or on any error. Never writes. Not thread safe, like the peptide look up of BioFSharp.Mz.
    type IsotopePatternLookUp(logger: NLog.Logger, cn: SQLiteConnection) =
        let mutable read = 0
        let mutable computed = 0
        let mutable errorLogged = false
        let settingsId =
            try
                use cmd = new SQLiteCommand("SELECT ID FROM IsotopePatternSettings WHERE Settings=@settings", cn)
                cmd.Parameters.AddWithValue("@settings", settings) |> ignore
                use reader = cmd.ExecuteReader()
                if reader.Read() then Some (reader.GetInt32 0) else None
            with
            | :? SQLiteException -> None
        let select =
            match settingsId with
            | Some id ->
                let cmd = new SQLiteCommand("SELECT Peaks FROM IsotopePattern WHERE SettingsID=@settingsId AND ModSequenceID=@modSequenceId AND Charge=@charge", cn)
                cmd.Parameters.AddWithValue("@settingsId", id) |> ignore
                cmd.Parameters.Add("@modSequenceId", DbType.Int32) |> ignore
                cmd.Parameters.Add("@charge", DbType.Int32) |> ignore
                Some cmd
            | None ->
                logger.Info "The peptide data base holds no isotope patterns of this version, they are computed."
                None
        let tryFind (lookUpResult: LookUpResult<AminoAcids.AminoAcid>) (charge: int) =
            match select with
            | Some cmd when lookUpResult.ModSequenceID >= 0 ->
                try
                    cmd.Parameters.["@modSequenceId"].Value <- lookUpResult.ModSequenceID
                    cmd.Parameters.["@charge"].Value <- charge
                    use reader = cmd.ExecuteReader()
                    if reader.Read() then decode (reader.["Peaks"] :?> byte[]) else None
                with ex ->
                    if not errorLogged then
                        errorLogged <- true
                        logger.Warn (sprintf "Reading an isotope pattern failed, the patterns are computed: %s" ex.Message)
                    None
            | _ -> None

        member this.Get (lookUpResult: LookUpResult<AminoAcids.AminoAcid>, charge: int) =
            match tryFind lookUpResult charge with
            | Some peaks ->
                read <- read + 1
                peaks
            | None ->
                computed <- computed + 1
                compute charge lookUpResult.BioSequence

        /// Number of patterns read from the data base.
        member this.Read = read

        /// Number of patterns computed.
        member this.Computed = computed

        interface IDisposable with
            member this.Dispose() =
                select |> Option.iter (fun cmd -> cmd.Dispose())

    /// Writes the isotope patterns that the table IsotopePattern of the peptide data base in the file dbPath lacks for the (sequence, global modification, charge) ions: for the global
    /// modification of the ion, or for both if the data base has isotopic modifications (ions that differ only in the global modification are computed once). The patterns are looked up
    /// and computed without a write lock and inserted in one transaction. If none is missing nothing is written, except the index on (Sequence, GlobalMod) if the data base lacks it.
    /// Ions whose look up or pattern fails are skipped and counted in the log. A failure to read or write the data base is a warning; the tools reading it compute what is missing.
    let store (logger: NLog.Logger) (dbPath: string) (ions: (string*int*int)[]) =
        try
            use cn = new SQLiteConnection(sprintf "Data Source=%s;Version=3;FailIfMissing=True;Default Timeout=600;BusyTimeout=600000" dbPath)
            cn.Open()
            let hasIndex =
                use cmd = new SQLiteCommand("SELECT COUNT(*) FROM sqlite_master WHERE type='index' AND name='SequenceAndGlobalModIndex'", cn)
                Convert.ToInt32 (cmd.ExecuteScalar()) > 0
            if not hasIndex then
                Db.SQLiteQuery.setIndexOnModSequenceAndGlobalMod cn |> ignore
            let sdbParams = getSDBParamsByCn cn
            let lookUp = getThreadSafePeptideLookUpFromFileBySequenceAndGMod cn sdbParams
            let labelled = not sdbParams.IsotopicMod.IsEmpty
            let existing =
                try
                    use cmd = new SQLiteCommand("SELECT p.ModSequenceID, p.Charge FROM IsotopePattern p JOIN IsotopePatternSettings s ON s.ID = p.SettingsID WHERE s.Settings=@settings", cn)
                    cmd.Parameters.AddWithValue("@settings", settings) |> ignore
                    use reader = cmd.ExecuteReader()
                    let keys = Collections.Generic.HashSet<struct (int*int)>()
                    while reader.Read() do
                        keys.Add (struct (reader.GetInt32 0, reader.GetInt32 1)) |> ignore
                    keys
                with
                | :? SQLiteException -> Collections.Generic.HashSet<struct (int*int)>()
            let lookUpFailed = ref 0
            let notFound = ref 0
            let patternFailed = ref 0
            let tryLookUp sequence gMod =
                try
                    Some (lookUp sequence gMod)
                with ex ->
                    if lookUpFailed.Value = 0 then
                        logger.Warn (sprintf "The look up of %s with global modification %i failed, the ion is skipped: %s" sequence gMod ex.Message)
                    lookUpFailed.Value <- lookUpFailed.Value + 1
                    None
            // The ID of the sequence look up is the key, it is the lowest ID of the rows with the sequence and the global modification as it is for the readers.
            let missing =
                (if labelled then ions |> Array.distinctBy (fun (sequence, _, charge) -> sequence, charge) else ions)
                |> Array.collect (fun (sequence, globalMod, charge) ->
                    (if labelled then [|0; 1|] else [|globalMod|])
                    |> Array.choose (fun gMod -> tryLookUp sequence gMod |> Option.map (fun lookUpResult -> lookUpResult, charge)))
                |> Array.filter (fun (lookUpResult, charge) ->
                    if lookUpResult.ModSequenceID < 0 then
                        notFound.Value <- notFound.Value + 1
                        false
                    else not (existing.Contains (struct (lookUpResult.ModSequenceID, charge))))
                |> Array.distinctBy (fun (lookUpResult, charge) -> lookUpResult.ModSequenceID, charge)
            let computedPatterns =
                missing
                |> Array.choose (fun (lookUpResult, charge) ->
                    try
                        Some (lookUpResult.ModSequenceID, charge, encode (compute charge lookUpResult.BioSequence))
                    with ex ->
                        if patternFailed.Value = 0 then
                            logger.Warn (sprintf "The isotope pattern of %s at charge %i failed, the ion is skipped: %s" lookUpResult.StringSequence charge ex.Message)
                        patternFailed.Value <- patternFailed.Value + 1
                        None)
            let skipped = sprintf "%i not in the data base, %i look ups and %i patterns failed" notFound.Value lookUpFailed.Value patternFailed.Value
            if computedPatterns.Length > 0 then
                use tr = cn.BeginTransaction(IsolationLevel.Serializable)
                let execute (query: string) =
                    use cmd = new SQLiteCommand(query, cn, tr)
                    cmd.ExecuteNonQuery() |> ignore
                execute "CREATE TABLE IF NOT EXISTS IsotopePatternSettings (ID INTEGER PRIMARY KEY, Settings TEXT NOT NULL UNIQUE)"
                execute "CREATE TABLE IF NOT EXISTS IsotopePattern (SettingsID INTEGER NOT NULL REFERENCES IsotopePatternSettings(ID), ModSequenceID INTEGER NOT NULL REFERENCES ModSequence(ID), Charge INTEGER NOT NULL, Peaks BLOB NOT NULL, PRIMARY KEY (SettingsID, ModSequenceID, Charge)) WITHOUT ROWID"
                use insertSettings = new SQLiteCommand("INSERT OR IGNORE INTO IsotopePatternSettings (Settings) VALUES (@settings)", cn, tr)
                insertSettings.Parameters.AddWithValue("@settings", settings) |> ignore
                insertSettings.ExecuteNonQuery() |> ignore
                use selectSettings = new SQLiteCommand("SELECT ID FROM IsotopePatternSettings WHERE Settings=@settings", cn, tr)
                selectSettings.Parameters.AddWithValue("@settings", settings) |> ignore
                let settingsId = Convert.ToInt32 (selectSettings.ExecuteScalar())
                use insert = new SQLiteCommand("INSERT OR IGNORE INTO IsotopePattern (SettingsID, ModSequenceID, Charge, Peaks) VALUES (@settingsId, @modSequenceId, @charge, @peaks)", cn, tr)
                insert.Parameters.AddWithValue("@settingsId", settingsId) |> ignore
                insert.Parameters.Add("@modSequenceId", DbType.Int32) |> ignore
                insert.Parameters.Add("@charge", DbType.Int32) |> ignore
                insert.Parameters.Add("@peaks", DbType.Binary) |> ignore
                let inserted =
                    computedPatterns
                    |> Array.sumBy (fun (modSequenceId, charge, peaks) ->
                        insert.Parameters.["@modSequenceId"].Value <- modSequenceId
                        insert.Parameters.["@charge"].Value <- charge
                        insert.Parameters.["@peaks"].Value <- peaks
                        insert.ExecuteNonQuery())
                tr.Commit()
                let bytes = computedPatterns |> Array.sumBy (fun (_, _, peaks) -> peaks.Length)
                logger.Info (sprintf "Isotope patterns in %s: %i ions, %i computed, %i added, %.1f peaks on average, %s" dbPath ions.Length computedPatterns.Length inserted (float bytes / 16. / float computedPatterns.Length) skipped)
            else
                logger.Info (sprintf "Isotope patterns in %s: %i ions, no pattern missing, %s" dbPath ions.Length skipped)
        with ex ->
            logger.Warn (sprintf "The isotope patterns could not be stored in %s, the tools reading it compute them: %s" dbPath ex.Message)
