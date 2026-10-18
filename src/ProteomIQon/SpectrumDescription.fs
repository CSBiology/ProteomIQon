namespace ProteomIQon

// A streaming parser of the spectrum description that reads only the fields the search needs.
module SpectrumDescription =

    open Newtonsoft.Json

    type SpectrumHeader =
        {
            ID          : string
            /// MS:1000511 of the spectrum, None when absent.
            MsLevel     : int option
            /// MS:1000016 of the first scan that carries it.
            ScanTime    : float option
            /// Of the first selected ion that carries one of the two.
            PrecursorMz : float option
            /// MS:1000041 of the first selected ion that carries it.
            ChargeState : int option
            /// MS:1002815 of the first scan that carries it.
            IonMobility : float option
            /// The message of the failure when the description could not be parsed, else null.
            /// Such a header carries no other value.
            Error       : string
        }

    let private parseFloat (text: string) =
        match System.Double.TryParse(text, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture) with
        | true, v -> Some v
        | _ -> None

    /// Reads the first entry of the Values array of the CV parameter object the reader stands in,
    /// and leaves the reader at the end of that object.
    let private firstValue (reader: JsonTextReader) =
        let mutable value : string = null
        let mutable depth = 0
        let mutable started = false
        while (not started || depth > 0) && reader.Read() do
            if reader.TokenType = JsonToken.StartObject then started <- true
            match reader.TokenType with
            | JsonToken.StartObject -> depth <- depth + 1
            | JsonToken.EndObject -> depth <- depth - 1
            | JsonToken.PropertyName when depth = 1 && (reader.Value :?> string) = "Values" ->
                if reader.Read() && reader.TokenType = JsonToken.StartArray then
                    if reader.Read() && reader.TokenType <> JsonToken.EndArray then
                        value <- (match reader.Value with null -> null | v -> System.Convert.ToString(v, System.Globalization.CultureInfo.InvariantCulture))
                        // skip the rest of the array
                        let mutable arrayDepth = 1
                        while arrayDepth > 0 && reader.Read() do
                            match reader.TokenType with
                            | JsonToken.StartArray -> arrayDepth <- arrayDepth + 1
                            | JsonToken.EndArray -> arrayDepth <- arrayDepth - 1
                            | _ -> ()
            | _ -> ()
        value

    /// Parses one description. The path of property names from the root tells where a CV
    /// parameter sits. A parameter of the spectrum itself has the shortest path. A scan and a
    /// selected ion each add their own levels to it.
    let parse (json: string) =
        use reader = new JsonTextReader(new System.IO.StringReader(json))
        reader.DateParseHandling <- DateParseHandling.None
        reader.FloatParseHandling <- FloatParseHandling.Double
        let path = System.Collections.Generic.List<string>()
        let mutable rootStarted = false
        let mutable rootClosed = false
        let mutable pending : string = null
        let mutable id : string = null
        let mutable msLevel = None
        let mutable scanTime = None
        let mutable mobility = None
        let mutable precursorMz = None
        let mutable charge = None
        // per selected ion
        let mutable ionPrecursorMz = None
        let mutable ionSelectedMz = None
        // path.[0] is the root object
        let inSelectedIon () =
            path.Count = 8 && path.[1] = "Precursors" && path.[2] = "properties" && path.[4] = "SelectedIons" && path.[5] = "properties" && path.[7] = "properties"
        let inScan () =
            path.Count = 5 && path.[1] = "Scans" && path.[2] = "properties" && path.[4] = "properties"
        let inSpectrum () = path.Count = 2 && path.[1] = "properties"
        let closeIon () =
            if precursorMz.IsNone then
                precursorMz <- (match ionPrecursorMz with Some v -> Some v | None -> ionSelectedMz)
            ionPrecursorMz <- None
            ionSelectedMz <- None
        while reader.Read() do
            match reader.TokenType with
            | JsonToken.PropertyName ->
                let name = reader.Value :?> string
                pending <- name
                if path.Count = 1 && name = "ID" then
                    if reader.Read() then id <- (reader.Value :?> string)
                elif inSpectrum () && name = "MS:1000511" then
                    let v = firstValue reader
                    if msLevel.IsNone && not (isNull v) then msLevel <- parseFloat v |> Option.map int
                elif inScan () && name = "MS:1000016" then
                    let v = firstValue reader
                    if scanTime.IsNone && not (isNull v) then scanTime <- parseFloat v
                elif inScan () && name = "MS:1002815" then
                    let v = firstValue reader
                    if mobility.IsNone && not (isNull v) then mobility <- parseFloat v
                elif inSelectedIon () && name = "MS:1002234" then
                    let v = firstValue reader
                    if not (isNull v) then ionPrecursorMz <- parseFloat v
                elif inSelectedIon () && name = "MS:1000744" then
                    let v = firstValue reader
                    if not (isNull v) then ionSelectedMz <- parseFloat v
                elif inSelectedIon () && name = "MS:1000041" then
                    let v = firstValue reader
                    if charge.IsNone && not (isNull v) then charge <- parseFloat v |> Option.map int
            | JsonToken.StartObject ->
                if path.Count = 0 then
                    if rootStarted then failwith "More than one JSON root object."
                    rootStarted <- true
                path.Add(if isNull pending then "" else pending)
                pending <- null
            | JsonToken.EndObject ->
                if inSelectedIon () then closeIon ()
                if path.Count > 0 then path.RemoveAt(path.Count - 1)
                if path.Count = 0 then rootClosed <- true
            | JsonToken.StartArray ->
                // arrays outside Values are not expected; skip them whole
                let mutable depth = 1
                while depth > 0 && reader.Read() do
                    match reader.TokenType with
                    | JsonToken.StartArray -> depth <- depth + 1
                    | JsonToken.EndArray -> depth <- depth - 1
                    | _ -> ()
                pending <- null
            | _ -> pending <- null
        if not rootStarted || not rootClosed || path.Count <> 0 then
            failwith "Incomplete spectrum description JSON."
        {
            ID = id
            MsLevel = msLevel
            ScanTime = scanTime
            PrecursorMz = precursorMz
            ChargeState = charge
            IonMobility = mobility
            Error = null
        }

    /// parse, with a failure turned into a header whose Error names it.
    let tryParse (json: string) =
        try parse json
        with ex -> { ID = null; MsLevel = None; ScanTime = None; PrecursorMz = None; ChargeState = None; IonMobility = None; Error = ex.Message }


/// Binary searches on sorted float arrays.
module SortedSearch =

    /// Index of the first value >= v in the sorted array.
    let lowerBound (sorted: float []) (v: float) =
        let mutable lo = 0
        let mutable hi = sorted.Length
        while lo < hi do
            let mid = (lo + hi) / 2
            if sorted.[mid] < v then lo <- mid + 1 else hi <- mid
        lo

    /// Index of the first value > v in the sorted array.
    let upperBound (sorted: float []) (v: float) =
        let mutable lo = 0
        let mutable hi = sorted.Length
        while lo < hi do
            let mid = (lo + hi) / 2
            if sorted.[mid] <= v then lo <- mid + 1 else hi <- mid
        lo

    /// Index of the value of the sorted, non-empty array closest to t: the first of equal values, and the
    /// earlier value when two are equally close (as Array.minBy over the array). A NaN or -infinity target gives
    /// the first index and +infinity the last one, where Array.minBy gives the first.
    let closestIndex (sorted: float []) (t: float) =
        let i = lowerBound sorted t
        let closest =
            if i = 0 then 0
            elif i >= sorted.Length then sorted.Length - 1
            elif sorted.[i] - t < t - sorted.[i - 1] then i
            else i - 1
        let mutable first = closest
        while first > 0 && sorted.[first - 1] = sorted.[closest] do
            first <- first - 1
        first

/// One pass readers of the spectrum descriptions of a mzlite run.
module SpectrumIndex =

    open System
    open System.Data.SQLite
    open SpectrumDescription

    /// Neither NaN nor infinite (netstandard2.0 has no Double.IsFinite).
    let private isFinite (x: float) = not (Double.IsNaN x || Double.IsInfinity x)

    /// The headers of the run in file order. A description that cannot be read is given with its Error.
    let private readHeaders (reader: MzIO.MzSQL.MzSQL) (runID: string) =
        seq {
            use cmd = new SQLiteCommand("SELECT Description FROM Spectrum WHERE RunID = @runID ORDER BY rowid", reader.Connection)
            cmd.Parameters.AddWithValue("@runID", runID) |> ignore
            use rows = cmd.ExecuteReader()
            while rows.Read() do yield tryParse (rows.GetString 0)
        }

    /// Scan time and ID of the MS1 spectra of the run, sorted by scan time (file order on equal scan times), from one
    /// pass over the spectrum descriptions with the streaming parser instead of the object model. A description that
    /// cannot be read and a MS1 without a finite scan time or an ID are not used; one warning gives their number. Fails
    /// if no MS1 is left.
    /// onOther is called with every readable header that is not of a MS1.
    let readMs1Index (reader: MzIO.MzSQL.MzSQL) (runID: string) (warn: string -> unit) (onOther: SpectrumHeader -> unit) =
        let entries = ResizeArray<float*string>()
        let mutable unreadable = 0
        for header in readHeaders reader runID do
            if not (isNull header.Error) then
                unreadable <- unreadable + 1
            elif header.MsLevel = Some 1 then
                match header.ScanTime with
                | Some scanTime when isFinite scanTime && not (String.IsNullOrWhiteSpace header.ID) -> entries.Add (scanTime, header.ID)
                | _ -> unreadable <- unreadable + 1
            else
                onOther header
        if unreadable > 0 then
            warn (sprintf "%i spectrum descriptions cannot be read or lack the ID or scan time; these spectra are not used." unreadable)
        if entries.Count = 0 then
            failwith "The run has no MS1 spectrum with an ID and a finite scan time, which the quantification needs."
        entries.ToArray() |> Array.sortBy fst

    /// The values the search needs from a MS2 spectrum description. Position is the place of the spectrum in the
    /// file, which orders spectra that share a scan time.
    type Ms2Header =
        {
            Id          : string
            Position    : int
            ScanTime    : float
            PrecursorMz : float
            /// 0 when unknown.
            Charge      : int
            /// nan when absent.
            IonMobility : float
        }

    /// The MS2 headers of the run in scan time order, equal scan times in file order, from one pass over the spectrum
    /// descriptions. A description that cannot be read and a MS2 without ID, finite scan time, precursor m/z above the
    /// proton mass or with a negative charge are logged with log and skipped.
    let readMs2Headers (reader: MzIO.MzSQL.MzSQL) (runID: string) (log: string -> unit) =
        let headers = ResizeArray<Ms2Header>()
        let mutable position = 0
        for h in readHeaders reader runID do
            if not (isNull h.Error) then
                log (sprintf "spectrum at position %i cannot be read: %s" position h.Error)
            elif h.MsLevel = Some 2 &&
                 (String.IsNullOrWhiteSpace h.ID ||
                  not (h.ScanTime |> Option.exists isFinite) ||
                  not (h.PrecursorMz |> Option.exists (fun x -> isFinite x && x > BioFSharp.Mass.Table.PMassInU)) ||
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
        ordered

    /// The item of (scan time, item) pairs sorted by scan time whose scan time is closest to a given one; the first
    /// of equal scan times, the earlier item on equal distance.
    let initClosestMs1 (ms1s: (float*'a) []) =
        let scanTimes = ms1s |> Array.map fst
        fun (scanTime: float) -> snd ms1s.[SortedSearch.closestIndex scanTimes scanTime]
