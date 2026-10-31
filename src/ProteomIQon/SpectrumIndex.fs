namespace ProteomIQon

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

/// One pass readers of the spectrum descriptions of a run.
module SpectrumIndex =

    open System
    open MzIO.IO
    open MzIO.Processing
    open MzIO.MzSQL.SpectrumHeader

    /// Neither NaN nor infinite (netstandard2.0 has no Double.IsFinite).
    let private isFinite (x: float) = not (Double.IsNaN x || Double.IsInfinity x)

    /// The headers of the run in file order. A description that cannot be read is given with its Error. An mzlite
    /// run is read from the description text one spectrum at a time; the headers of any other reader come from the
    /// object model and carry no charge state and no ion mobility.
    let readHeaders (reader: IMzIODataReader) (runID: string) =
        match reader with
        | :? MzIO.MzSQL.MzSQL as sql -> sql.ReadSpectrumHeadersLazy runID
        | _ ->
            reader.ReadMassSpectra(runID)
            |> Seq.map (fun ms ->
                let msLevel = MassSpectrum.getMsLevel ms
                {
                    ID = ms.ID
                    MsLevel = Some msLevel
                    // MzIO gives -1 for a spectrum without a scan start time
                    ScanTime = (match MassSpectrum.getScanTime ms with -1. -> None | scanTime -> Some scanTime)
                    PrecursorMz = (if msLevel = 2 then Some (MassSpectrum.getPrecursorMZ ms) else None)
                    ChargeState = None
                    IonMobility = None
                    Error = null
                })

    /// Scan time and ID of a MS1 header, None without an ID or a finite scan time.
    let private ms1Entry (h: SpectrumHeader) =
        match h.ScanTime with
        | Some scanTime when isFinite scanTime && not (String.IsNullOrWhiteSpace h.ID) -> Some (scanTime, h.ID)
        | _ -> None

    /// Sorted by scan time, equal scan times in file order (Seq.sortBy is stable).
    let private sortMs1 (entries: ResizeArray<float*string>) =
        entries |> Seq.sortBy fst |> Array.ofSeq

    let private noMs1 = "The run has no MS1 spectrum with an ID and a finite scan time, which the quantification needs."

    /// Scan time and ID of the MS1 spectra of the run, sorted by scan time (file order on equal scan times), from one
    /// pass over the spectrum descriptions with the streaming parser instead of the object model. A description that
    /// cannot be read and a MS1 without a finite scan time or an ID are not used; one warning gives their number. Fails
    /// if no MS1 is left.
    /// onOther is called with every readable header that is not of a MS1.
    let readMs1Index (reader: IMzIODataReader) (runID: string) (warn: string -> unit) (onOther: SpectrumHeader -> unit) =
        let entries = ResizeArray<float*string>()
        let mutable unreadable = 0
        for header in readHeaders reader runID do
            if not (isNull header.Error) then
                unreadable <- unreadable + 1
            elif header.MsLevel = Some 1 then
                match ms1Entry header with
                | Some entry -> entries.Add entry
                | None -> unreadable <- unreadable + 1
            else
                onOther header
        if unreadable > 0 then
            warn (sprintf "%i spectrum descriptions cannot be read or lack the ID or scan time; these spectra are not used." unreadable)
        if entries.Count = 0 then failwith noMs1
        sortMs1 entries

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

    /// The MS2 header of a readable MS2 description, None (logged with log) without ID, finite scan time, precursor m/z
    /// above the proton mass or with a negative charge.
    let private ms2Header (log: string -> unit) (position: int) (h: SpectrumHeader) =
        if String.IsNullOrWhiteSpace h.ID ||
           not (h.ScanTime |> Option.exists isFinite) ||
           not (h.PrecursorMz |> Option.exists (fun x -> isFinite x && x > BioFSharp.Mass.Table.PMassInU)) ||
           (h.ChargeState |> Option.exists (fun x -> x < 0)) then
            log (sprintf "spectrum at position %i has missing or invalid search metadata; skipped." position)
            None
        else
            Some
                {
                    Id = h.ID
                    Position = position
                    ScanTime = defaultArg h.ScanTime -1.
                    PrecursorMz = defaultArg h.PrecursorMz -1.
                    Charge = defaultArg h.ChargeState 0
                    IonMobility = defaultArg h.IonMobility nan
                }

    let private sortMs2 (headers: ResizeArray<Ms2Header>) =
        let ordered = headers.ToArray()
        Array.sortInPlaceWith (fun a b -> compare (a.ScanTime, a.Position) (b.ScanTime, b.Position)) ordered
        ordered

    /// The MS2 headers of the run in scan time order, equal scan times in file order, from one pass over the spectrum
    /// descriptions. A description that cannot be read and a MS2 without ID, finite scan time, precursor m/z above the
    /// proton mass or with a negative charge are logged with log and skipped.
    let readMs2Headers (reader: IMzIODataReader) (runID: string) (log: string -> unit) =
        let headers = ResizeArray<Ms2Header>()
        let mutable position = 0
        for h in readHeaders reader runID do
            if not (isNull h.Error) then
                log (sprintf "spectrum at position %i cannot be read: %s" position h.Error)
            elif h.MsLevel = Some 2 then
                ms2Header log position h |> Option.iter headers.Add
            position <- position + 1
        sortMs2 headers

    /// The MS1 index of readMs1Index and the MS2 headers of readMs2Headers from one pass over the run. A description
    /// that cannot be read is logged with log; the MS1 spectra without an ID or a finite scan time are counted in one
    /// warning. The MS1 index can be empty.
    let readMs1AndMs2 (reader: IMzIODataReader) (runID: string) (warn: string -> unit) (log: string -> unit) =
        let ms1s = ResizeArray<float*string>()
        let ms2s = ResizeArray<Ms2Header>()
        let mutable skippedMs1 = 0
        let mutable position = 0
        for h in readHeaders reader runID do
            if not (isNull h.Error) then
                log (sprintf "spectrum at position %i cannot be read: %s" position h.Error)
            elif h.MsLevel = Some 1 then
                match ms1Entry h with
                | Some entry -> ms1s.Add entry
                | None -> skippedMs1 <- skippedMs1 + 1
            elif h.MsLevel = Some 2 then
                ms2Header log position h |> Option.iter ms2s.Add
            position <- position + 1
        if skippedMs1 > 0 then
            warn (sprintf "%i MS1 spectra lack the ID or a finite scan time; these spectra are not used." skippedMs1)
        sortMs1 ms1s, sortMs2 ms2s

    /// The item of (scan time, item) pairs sorted by scan time whose scan time is closest to a given one; the first
    /// of equal scan times, the earlier item on equal distance.
    let initClosestMs1 (ms1s: (float*'a) []) =
        let scanTimes = ms1s |> Array.map fst
        fun (scanTime: float) -> snd ms1s.[SortedSearch.closestIndex scanTimes scanTime]
