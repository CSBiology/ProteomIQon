namespace ProteomIQon

open System
open System.Collections.Generic
open System.Data.SQLite
open System.Threading.Tasks
open BioFSharp
open BioFSharp.Mz

/// Peptide table and fragment ion index in the style of MSFragger. The peptides are kept in
/// ascending mass order, so the entries of every fragment bin are in mass order as well and
/// the precursor window of a spectrum becomes a binary search on peptide ranks. Only target
/// fragments are stored. The fragments of a reversed decoy follow from the target fragments
/// by a water shift (decoy b ions are target y ions minus water, decoy y ions are target b
/// ions plus water), so the search serves decoys from the same index.
module PeptideIndex =

    /// Mask that strips the ion type flag from an index entry.
    [<Literal>]
    let RankMask = 0x3FFFFFFF

    /// Bit 30 marks a y ion entry.
    [<Literal>]
    let YFlag = 0x40000000

    /// A residue code stands for one residue together with the modification tokens that
    /// precede it in the ModSequence string, for example "[ac][ox]M".
    type ResidueCode =
        {
            Token  : string
            /// Neutral residue mass per GlobalMod (index 0 unlabeled, 1 with the isotopic label).
            Masses : float[]
        }

    /// Peptides sorted by neutral mass. Sequences are stored as residue codes in one flat array.
    type PeptideTable =
        {
            Count         : int
            ModSequenceId : int[]
            PepSequenceId : int[]
            Mass          : float[]
            GlobalMod     : byte[]
            SeqStart      : int[]
            Residues      : byte[]
            Codes         : ResidueCode[]
            /// Longest residue count of any peptide.
            MaxLength     : int
            /// ModSequence strings of the peptides whose string holds tokens without a residue
            /// (the terminator "*" and the gap "-"), keyed by rank. All other strings can be
            /// rebuilt from the residue codes.
            RawSequence   : IDictionary<int, string>
        }
        member this.Length (rank: int) = this.SeqStart.[rank + 1] - this.SeqStart.[rank]

        /// The ModSequence string of a peptide as stored in the database.
        member this.SequenceString (rank: int) =
            match this.RawSequence.TryGetValue rank with
            | true, raw -> raw
            | _ ->
                Array.sub this.Residues this.SeqStart.[rank] (this.Length rank)
                |> Array.map (fun code -> this.Codes.[int code].Token)
                |> String.concat ""

    /// Fragment index over the target b and y ions of all peptides of a PeptideTable.
    type FragmentIndex =
        {
            BinWidth   : float
            BinCount   : int
            /// Per bin: chunk holding the entries, local start offset, number of entries.
            BinChunk   : int[]
            BinLocal   : int[]
            BinLength  : int[]
            /// Entries are peptide ranks with the YFlag bit set for y ions.
            Chunks     : int[][]
            EntryCount : int64
        }

    let waterMass = Formula.parseFormulaString "H2O" |> Formula.monoisoMass

    /// Splits a ModSequence string such as "[ac][ox]MTAILER" into residue tokens. Every
    /// character outside square brackets is a residue. The terminator "*" and the gap "-"
    /// that BioFSharp writes for some database entries count as residues here as well.
    let tokenize (sequence: string) =
        let rec loop start pos acc =
            if pos >= sequence.Length then
                if pos > start then failwithf "Sequence '%s' ends with a modification token without residue." sequence
                else List.rev acc
            elif sequence.[pos] = '[' then
                match sequence.IndexOf(']', pos) with
                | -1 -> failwithf "Sequence '%s' has an unclosed modification token." sequence
                | close -> loop start (close + 1) acc
            else
                loop (pos + 1) (pos + 1) (sequence.Substring(start, pos + 1 - start) :: acc)
        loop 0 0 []

    /// Computes the residue masses of a token for every GlobalMod with the parser and the mass
    /// function of the search database, so labels and modifications follow the definitions
    /// that produced the RealMass column.
    let createResidueCode (sdbParams: SearchDB.SearchDbParams) (globalMods: int) (token: string) =
        let parse = SearchDB.initOfModAminoAcidString sdbParams.IsotopicMod (sdbParams.FixedMods @ sdbParams.VariableMods)
        let masses =
            Array.init globalMods (fun globalMod ->
                match parse globalMod token with
                | [ residue ] -> sdbParams.MassFunction residue
                | _ -> failwithf "Token '%s' does not parse to a single residue." token)
        { Token = token; Masses = masses }

    [<Struct>]
    type private PeptideRow =
        {
            ModSequenceId : int
            PepSequenceId : int
            RealMass      : float
            GlobalMod     : byte
            Sequence      : byte[]
            /// The database string when it holds tokens that are not residues, otherwise null.
            Raw           : string
        }

    /// Reads every ModSequence row of the peptide database into a mass sorted table. Tokens
    /// that the BioFSharp parser drops (the terminator "*" and the gap "-", both of mass 0)
    /// are left out, so the residues match what the scoring functions see.
    let loadPeptideTable (cn: SQLiteConnection) (sdbParams: SearchDB.SearchDbParams) (log: string -> unit) =
        let parse = SearchDB.initOfModAminoAcidString sdbParams.IsotopicMod (sdbParams.FixedMods @ sdbParams.VariableMods) 0
        let codeIndex = Dictionary<string, int option>()
        let codeOf (token: string) =
            match codeIndex.TryGetValue token with
            | true, code -> code
            | _ ->
                let code =
                    if List.isEmpty (parse token) then None
                    else
                        let assigned = codeIndex.Values |> Seq.choose id |> Seq.length
                        if assigned > 255 then failwith "More than 256 distinct residue tokens are not supported."
                        Some assigned
                codeIndex.[token] <- code
                code
        let encode (sequence: string) =
            let codes = tokenize sequence |> List.map codeOf
            let residues = codes |> List.choose (Option.map byte) |> Array.ofList
            residues, (if codes |> List.exists Option.isNone then sequence else null)
        log "Reading the ModSequence table."
        let rows =
            use cmd = new SQLiteCommand("SELECT ID, PepSequenceID, RealMass, Sequence, GlobalMod FROM ModSequence", cn)
            use reader = cmd.ExecuteReader()
            [|
                while reader.Read() do
                    let residues, raw = encode (reader.GetString 3)
                    yield
                        {
                            ModSequenceId = reader.GetInt32 0
                            PepSequenceId = reader.GetInt32 1
                            RealMass = reader.GetDouble 2
                            GlobalMod = byte (reader.GetInt32 4)
                            Sequence = residues
                            Raw = raw
                        }
            |]
            |> Array.sortBy (fun row -> row.RealMass)
        log (sprintf "%i mod sequences read and sorted by mass." rows.Length)
        let globalMods = 1 + (rows |> Array.fold (fun acc row -> max acc (int row.GlobalMod)) 0)
        let codes =
            codeIndex
            |> Seq.choose (fun kv -> kv.Value |> Option.map (fun code -> code, kv.Key))
            |> Seq.sortBy fst
            |> Seq.map (fun (_, token) -> createResidueCode sdbParams globalMods token)
            |> Array.ofSeq
        let total = rows |> Array.sumBy (fun row -> int64 row.Sequence.Length)
        if total > int64 Int32.MaxValue then failwith "The concatenated peptide sequences exceed the supported size."
        let seqStart = rows |> Array.scan (fun start row -> start + row.Sequence.Length) 0
        let rawSequences =
            rows
            |> Array.indexed
            |> Array.choose (fun (rank, row) -> if isNull row.Raw then None else Some (rank, row.Raw))
            |> dict
        let table =
            {
                Count = rows.Length
                ModSequenceId = rows |> Array.map (fun row -> row.ModSequenceId)
                PepSequenceId = rows |> Array.map (fun row -> row.PepSequenceId)
                Mass = rows |> Array.map (fun row -> row.RealMass)
                GlobalMod = rows |> Array.map (fun row -> row.GlobalMod)
                SeqStart = seqStart
                Residues = rows |> Array.collect (fun row -> row.Sequence)
                Codes = codes
                MaxLength = rows |> Array.fold (fun acc row -> max acc row.Sequence.Length) 0
                RawSequence = rawSequences
            }
        log (sprintf "Peptide table ready: %i peptides, %i residue codes, mass %.3f to %.3f." table.Count codes.Length table.Mass.[0] table.Mass.[table.Count - 1])
        table

    /// Neutral mass of a peptide from its residue codes.
    let peptideMass (table: PeptideTable) (rank: int) =
        let globalMod = int table.GlobalMod.[rank]
        Array.sub table.Residues table.SeqStart.[rank] (table.Length rank)
        |> Array.sumBy (fun code -> table.Codes.[int code].Masses.[globalMod])
        |> (+) waterMass

    /// Checks that the residue masses reproduce the RealMass column on a sample of peptides
    /// and returns the worst deviation.
    let verifyMasses (table: PeptideTable) (maxDeviation: float) =
        let step = max 1 (table.Count / 20000)
        let worst =
            [| 0 .. step .. table.Count - 1 |]
            |> Array.map (fun rank -> abs (peptideMass table rank - table.Mass.[rank]))
            |> Array.max
        if worst > maxDeviation then
            failwithf "Residue masses deviate from the database masses by up to %f Da." worst
        worst

    /// Writes the neutral b ion masses (b1 .. b(L-1)) into buffer positions 0 .. L-2 and the
    /// neutral y ion masses (y1 .. y(L-1)) into positions L-1 .. 2L-3. Returns 2(L-1).
    /// Runs for every peptide and every reported hit, so it fills a reused buffer in place.
    let fragmentMasses (table: PeptideTable) (rank: int) (buffer: float[]) =
        let start = table.SeqStart.[rank]
        let length = table.SeqStart.[rank + 1] - start
        let globalMod = int table.GlobalMod.[rank]
        let residueMass i = table.Codes.[int table.Residues.[start + i]].Masses.[globalMod]
        let mutable acc = 0.
        for i = 0 to length - 2 do
            acc <- acc + residueMass i
            buffer.[i] <- acc
        acc <- waterMass
        for i = 0 to length - 2 do
            acc <- acc + residueMass (length - 1 - i)
            buffer.[length - 1 + i] <- acc
        2 * (length - 1)

    /// Builds the fragment index with the given number of workers. Target b and y ions are
    /// binned by neutral mass. Entries inside a bin end up in peptide rank order because every
    /// worker owns a contiguous rank range and processes it in order, and the workers write to
    /// disjoint offsets inside every bin.
    let buildFragmentIndex (table: PeptideTable) (binWidth: float) (maxFragmentMass: float) (workers: int) (log: string -> unit) =
        let binCount = int (maxFragmentMass / binWidth) + 2
        let workers = max 1 (min workers 64)
        let bufferSize = 2 * max 1 table.MaxLength
        let rankOfWorker = Array.init (workers + 1) (fun w -> int (int64 table.Count * int64 w / int64 workers))
        let binOf (mass: float) =
            let bin = int (mass / binWidth)
            if bin < 0 || bin >= binCount then -1 else bin
        let visit (rank: int) (buffer: float[]) (f: int -> int -> unit) =
            let count = fragmentMasses table rank buffer
            let half = count / 2
            for i = 0 to half - 1 do
                let b = binOf buffer.[i]
                if b >= 0 then f b rank
                let y = binOf buffer.[half + i]
                if y >= 0 then f y (rank ||| YFlag)
        log (sprintf "Counting fragments in %i bins with %i workers." binCount workers)
        let counts = Array.init workers (fun _ -> Array.zeroCreate<int> binCount)
        Parallel.For(0, workers, fun w ->
            let buffer = Array.zeroCreate<float> bufferSize
            let c = counts.[w]
            for rank = rankOfWorker.[w] to rankOfWorker.[w + 1] - 1 do
                visit rank buffer (fun b _ -> c.[b] <- c.[b] + 1)) |> ignore
        let binTotal = Array.init binCount (fun b -> counts |> Array.sumBy (fun c -> int64 c.[b]))
        let binStart = Array.scan (+) 0L binTotal
        let total = binStart.[binCount]
        log (sprintf "%i fragment entries (%.2f GB)." total (float total * 4. / 1e9))
        // Bins are grouped into chunks so that no chunk exceeds the array length limit.
        let maxChunk = int64 (1 <<< 30)
        let binChunk, binLocal =
            binTotal
            |> Array.scan (fun (chunk, local, _) len ->
                if len > maxChunk then failwith "A single fragment bin exceeds the supported size. Use a smaller bin width."
                if local + len > maxChunk then (chunk + 1, len, 0L) else (chunk, local + len, local)) (0, 0L, 0L)
            |> Array.tail
            |> Array.map (fun (chunk, _, local) -> chunk, int local)
            |> Array.unzip
        let chunkSizes =
            Array.zip binChunk binTotal
            |> Array.groupBy fst
            |> Array.sortBy fst
            |> Array.map (fun (_, bins) -> bins |> Array.sumBy snd |> int)
        let chunks = chunkSizes |> Array.map (fun n -> Array.zeroCreate<int> n)
        // Worker w writes the entries of bin b behind the entries of the workers before it,
        // running sums over the workers in one pass per bin.
        let workerOffset = Array.init workers (fun _ -> Array.zeroCreate<int> binCount)
        for b = 0 to binCount - 1 do
            for w = 0 to workers - 1 do
                workerOffset.[w].[b] <- (if w = 0 then binLocal.[b] else workerOffset.[w - 1].[b] + counts.[w - 1].[b])
        log "Filling fragment bins."
        Parallel.For(0, workers, fun w ->
            let buffer = Array.zeroCreate<float> bufferSize
            let cursor = workerOffset.[w]
            for rank = rankOfWorker.[w] to rankOfWorker.[w + 1] - 1 do
                visit rank buffer (fun b entry ->
                    chunks.[binChunk.[b]].[cursor.[b]] <- entry
                    cursor.[b] <- cursor.[b] + 1)) |> ignore
        log (sprintf "Fragment index ready: %i chunks." chunks.Length)
        {
            BinWidth = binWidth
            BinCount = binCount
            BinChunk = binChunk
            BinLocal = binLocal
            BinLength = binTotal |> Array.map int
            Chunks = chunks
            EntryCount = total
        }

    /// First position in [lo, hi) of the chunk whose rank is at least the given rank.
    let lowerBound (chunk: int[]) (lo: int) (hi: int) (rank: int) =
        let rec loop l h =
            if l >= h then l
            else
                let m = l + ((h - l) >>> 1)
                if (chunk.[m] &&& RankMask) < rank then loop (m + 1) h else loop l m
        loop lo hi

    /// First rank whose mass is at least the given mass.
    let rankLowerBound (table: PeptideTable) (mass: float) =
        let rec loop l h =
            if l >= h then l
            else
                let m = l + ((h - l) >>> 1)
                if table.Mass.[m] < mass then loop (m + 1) h else loop l m
        loop 0 table.Count
