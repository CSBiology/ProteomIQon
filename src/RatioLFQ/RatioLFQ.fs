namespace ProteomIQon

open System
open System.IO
open System.Globalization
open System.Collections.Generic
open FSharpAux.IO.SchemaReader
open FSharp.Stats
open ProteomIQon.Domain

/// Light/heavy ratios per protein group and run from a peptide ion table of a labeled experiment.
/// The ratio profile of a protein group comes from pairwise comparisons of runs over shared ions
/// (the MaxLFQ construction applied to ion log2 ratios). The intensity is that ratio times one
/// reference intensity per protein group, taken from the spike-in channel.
module RatioLFQ =

    let private logger = Logging.createLogger "RatioLFQ"

    /// One valid ion measurement. Ion is an index local to its protein group, Run an index into Runs.
    [<Struct>]
    type IonMeasurement =
        {
            Ion   : int
            Run   : int
            Light : float
            Heavy : float
        }

    /// Measurements per protein group, in the order of ProteinGroups.
    type IonTable =
        {
            Runs          : string []
            ProteinGroups : string []
            Measurements  : IonMeasurement [] []
            RowsRead      : int
            RowsIgnored   : int
        }

    /// Ratio and Intensity hold one value per run, nan where the protein group has no valid ion.
    type ProteinResult =
        {
            ProteinGroup : string
            Ratio        : float []
            Intensity    : float []
            RunGroups    : int
        }

    let private isValidIntensity (intensity: float) = Numerics.isFinite intensity && intensity > 0.

    /// The index of the key, a new key gets the next one.
    let private indexOf (indices: Dictionary<'Key, int>) (key: 'Key) =
        match indices.TryGetValue key with
        | true, index -> index
        | _ ->
            let index = indices.Count
            indices.Add(key, index)
            index

    let private keysByIndex (indices: Dictionary<string, int>) =
        indices |> Seq.sortBy (fun entry -> entry.Value) |> Seq.map (fun entry -> entry.Key) |> Array.ofSeq

    /// Reads the ion table by its column names. Runs, protein groups and ions become integer indices while reading,
    /// only the compact measurements are kept.
    let readIonTable (path: string) : IonTable =
        Core.TableHeader.checkColumns<Dto.GlobModAggregationResult> "the GlobModAggregation.txt of LabeledProteinQuantification" path
        let runs = Dictionary<string, int>()
        let proteinGroups = Dictionary<string, int>()
        let ions = ResizeArray<Dictionary<struct (string * int), int>>()
        let measurements = ResizeArray<ResizeArray<IonMeasurement>>()
        let seen = HashSet<struct (int * int * int)>()
        let mutable rowsRead = 0
        let mutable rowsIgnored = 0
        for row in Csv.CsvReader<Dto.GlobModAggregationResult>(SchemaMode=Csv.Fill).ReadFile(path,'\t',true) do
            rowsRead <- rowsRead + 1
            // LabeledProteinQuantification quotes run names with a space and protein groups with a semicolon
            let file, proteinGroup = row.FileName.Trim '"', row.ProteinGroup.Trim '"'
            let run = if file = "" then -1 else indexOf runs file
            let group = if proteinGroup = "" then -1 else indexOf proteinGroups proteinGroup
            if group = measurements.Count then
                ions.Add(Dictionary<struct (string * int), int>())
                measurements.Add(ResizeArray<IonMeasurement>())
            if run >= 0 && group >= 0 && row.StringSequence <> "" && isValidIntensity row.Quant_Light && isValidIntensity row.Quant_Heavy then
                let ion = indexOf ions.[group] (struct (row.StringSequence, row.Charge))
                if not (seen.Add(struct (group, ion, run))) then
                    failwithf "The protein group %s has more than one valid row for the ion %s, charge %i in the run %s. RatioLFQ needs one row per ion and run, the table is not an aggregated ion table."
                        proteinGroup row.StringSequence row.Charge file
                measurements.[group].Add { Ion = ion; Run = run; Light = row.Quant_Light; Heavy = row.Quant_Heavy }
            else
                rowsIgnored <- rowsIgnored + 1
        // protein groups keep the order of their first appearance, also when that row was invalid
        let groupNames = keysByIndex proteinGroups
        let kept = [| for g in 0 .. groupNames.Length - 1 do if measurements.[g].Count > 0 then yield g |]
        {
            Runs          = keysByIndex runs
            ProteinGroups = kept |> Array.map (fun g -> groupNames.[g])
            Measurements  = kept |> Array.map (fun g -> measurements.[g].ToArray())
            RowsRead      = rowsRead
            RowsIgnored   = rowsIgnored
        }

    /// Plain ratio of every run: the median over the ions of the log2 ratio, nan for a run without an ion.
    /// y.[ion].[run] is the log2 ratio or nan.
    let plainRatio (y: float [] []) : float [] =
        let runs = if y.Length = 0 then 0 else y.[0].Length
        Array.init runs (fun run ->
            y
            |> Array.choose (fun ion -> if Double.IsNaN ion.[run] then None else Some ion.[run])
            |> Numerics.median)

    /// Profile per run (pairwise comparison and level) and the number of run groups.
    let private profile (weighting: PairWeighting) (minSharedIons: int) (y: float [] []) : float [] * int =
        let runs = if y.Length = 0 then 0 else y.[0].Length
        let m = plainRatio y
        let present = [| for r in 0 .. runs - 1 do if not (Double.IsNaN m.[r]) then yield r |]
        let k = present.Length
        // pairwise ratios over the present runs only, each valid pair weighted by its number of shared ions or by 1
        let pair = Array2D.create k k nan
        let weight = Array2D.zeroCreate<float> k k
        for a in 0 .. k - 1 do
            for b in a + 1 .. k - 1 do
                let diffs =
                    y
                    |> Array.choose (fun ion ->
                        let ya, yb = ion.[present.[a]], ion.[present.[b]]
                        if Double.IsNaN ya || Double.IsNaN yb then None else Some (ya - yb))
                if diffs.Length >= minSharedIons then
                    let r = Numerics.median diffs
                    pair.[a, b] <- r
                    pair.[b, a] <- -r
                    let w = match weighting with PairWeighting.SharedIonCount -> float diffs.Length | PairWeighting.Unweighted -> 1.
                    weight.[a, b] <- w
                    weight.[b, a] <- w
        let x = Array.create runs nan
        let runGroup = Array.create k -1
        let mutable groups = 0
        for start in 0 .. k - 1 do
            if runGroup.[start] < 0 then
                let members = ResizeArray<int>()
                let queue = Queue<int>()
                runGroup.[start] <- groups
                queue.Enqueue start
                while queue.Count > 0 do
                    let a = queue.Dequeue()
                    members.Add a
                    for b in 0 .. k - 1 do
                        if runGroup.[b] < 0 && not (Double.IsNaN pair.[a, b]) then
                            runGroup.[b] <- groups
                            queue.Enqueue b
                groups <- groups + 1
                let idx = members.ToArray() |> Array.sort
                let size = idx.Length
                if size = 1 then
                    x.[present.[idx.[0]]] <- m.[present.[idx.[0]]]
                else
                    // weighted normal equations L z = b, the added J/k fixes the free additive constant
                    let anchor = 1. / float size
                    let laplacian =
                        Matrix.init size size (fun p q ->
                            if p = q then (idx |> Array.sumBy (fun b -> weight.[idx.[p], b])) + anchor
                            else anchor - weight.[idx.[p], idx.[q]])
                    let b =
                        Vector.init size (fun p ->
                            idx |> Array.sumBy (fun q -> if Double.IsNaN pair.[idx.[p], q] then 0. else weight.[idx.[p], q] * pair.[idx.[p], q]))
                    let z = Algebra.LinearAlgebra.SolveLinearSystem laplacian b |> Vector.toArray
                    let meanZ = Array.average z
                    let z = z |> Array.map (fun v -> v - meanZ)
                    // the level counts every ion measurement of the run group once, so a few poor runs cannot shift it
                    let level =
                        [| for ion in y do
                             for p in 0 .. size - 1 do
                                 let v = ion.[present.[idx.[p]]]
                                 if not (Double.IsNaN v) then yield v - z.[p] |]
                        |> Numerics.median
                    for p in 0 .. size - 1 do
                        x.[present.[idx.[p]]] <- z.[p] + level
        x, groups

    /// Log2 ratio per run from pairwise comparisons over shared ions, weighted by the shared ion count, placed at the median level of the run group, nan for runs without an ion.
    /// y.[ion].[run] is the log2 ratio or nan.
    let pairwiseProfile (minSharedIons: int) (y: float [] []) : float [] =
        profile PairWeighting.SharedIonCount minSharedIons y |> fst

    /// Reference level of a protein group from the per-ion medians of the reference channel: the mean of all of them or the mean of the up to three largest, as the rule says.
    /// reference.[ion].[run] is the linear reference intensity or nan.
    let referenceLevel (rule: ReferenceLevel) (reference: float [] []) : float =
        let h =
            reference
            |> Array.choose (fun ion ->
                let values = ion |> Array.filter (fun v -> not (Double.IsNaN v))
                if values.Length = 0 then None else Some (Numerics.median values))
            |> Array.sortDescending
        if h.Length = 0 then nan
        else
            match rule with
            | ReferenceLevel.MeanOfIonMedians -> Array.average h
            | ReferenceLevel.TopThreeIonMedians -> h |> Array.truncate 3 |> Array.average

    /// Centring, profile, reference level and the output rows for every protein group. The matrices of one protein group cover its present runs
    /// only and are dropped after the group.
    let compute (parameters: RatioLFQParams) (table: IonTable) : ProteinResult [] =
        if parameters.MinSharedIons < 1 then
            failwithf "MinSharedIons has to be at least 1, got %i." parameters.MinSharedIons
        let log2Ratio (m: IonMeasurement) = Math.Log2(m.Light / m.Heavy)
        let shift =
            if parameters.CentreRuns then
                let perRun = Array.init table.Runs.Length (fun _ -> ResizeArray<float>())
                for group in table.Measurements do
                    for m in group do
                        perRun.[m.Run].Add(log2Ratio m)
                perRun |> Array.map (fun values -> if values.Count = 0 then 0. else Numerics.median (values.ToArray()))
            else
                Array.zeroCreate table.Runs.Length
        Array.map2 (fun name (ms: IonMeasurement []) ->
            let presentRuns = ms |> Array.map (fun m -> m.Run) |> Array.distinct |> Array.sort
            let column = Dictionary<int, int>()
            presentRuns |> Array.iteri (fun c run -> column.[run] <- c)
            let ions = (ms |> Array.map (fun m -> m.Ion) |> Array.max) + 1
            let y = Array.init ions (fun _ -> Array.create presentRuns.Length nan)
            let reference = Array.init ions (fun _ -> Array.create presentRuns.Length nan)
            for m in ms do
                y.[m.Ion].[column.[m.Run]] <- log2Ratio m - shift.[m.Run]
                reference.[m.Ion].[column.[m.Run]] <-
                    (match parameters.Reference with
                     | ReferenceChannel.Heavy -> m.Heavy
                     | ReferenceChannel.Light -> m.Light)
            let x, groups = profile parameters.PairWeighting parameters.MinSharedIons y
            let level = referenceLevel parameters.ReferenceLevel reference
            let ratio = Array.create table.Runs.Length nan
            let intensity = Array.create table.Runs.Length nan
            presentRuns |> Array.iteri (fun c run ->
                ratio.[run] <- Math.Pow(2., x.[c])
                intensity.[run] <-
                    (match parameters.Reference with
                     | ReferenceChannel.Heavy -> Math.Pow(2., x.[c]) * level
                     | ReferenceChannel.Light -> Math.Pow(2., -x.[c]) * level))
            { ProteinGroup = name; Ratio = ratio; Intensity = intensity; RunGroups = groups }
        ) table.ProteinGroups table.Measurements

    /// Reads the ion table, computes the ratios and intensities and writes RatioLFQ.txt into the output directory.
    let ratioLFQ (parameters: RatioLFQParams) (outputDirectory: string) (ionTable: string) : unit =
        if parameters.MinSharedIons < 1 then
            failwithf "MinSharedIons has to be at least 1, got %i." parameters.MinSharedIons
        let warn (message: string) =
            logger.Warn message
            logger.Info (sprintf "Warning: %s" message)
        logger.Trace (sprintf "Parameters: %A" parameters)
        let table = readIonTable ionTable
        logger.Trace (sprintf "Rows read: %i, rows ignored: %i, runs: %i" table.RowsRead table.RowsIgnored table.Runs.Length)
        let mutable validCount = 0
        let mutable largestLight = 0.
        let mutable largestHeavy = 0.
        for group in table.Measurements do
            for m in group do
                validCount <- validCount + 1
                largestLight <- max largestLight m.Light
                largestHeavy <- max largestHeavy m.Heavy
        if validCount = 0 then
            warn "The ion table has no valid row. The output file holds the header only."
        elif largestLight < 100. || largestHeavy < 100. then
            failwithf "The largest Quant_Light is %g and the largest Quant_Heavy is %g. The table looks transformed (for example log2). RatioLFQ needs untransformed intensities." largestLight largestHeavy
        let results = compute parameters table
        let cellText (v: float) = if Double.IsNaN v then "" else v.ToString(CultureInfo.InvariantCulture)
        Directory.CreateDirectory outputDirectory |> ignore
        use writer = new StreamWriter(Path.Combine(outputDirectory, "RatioLFQ.txt"))
        let header =
            [ yield "ProteinGroup"
              for run in table.Runs do
                  yield run + ".Ratio_LightByHeavy"
                  yield run + ".LFQ_Intensity" ]
        writer.WriteLine(String.Join("\t", header))
        for r in results do
            let cells =
                [ yield r.ProteinGroup
                  for run in 0 .. table.Runs.Length - 1 do
                      yield cellText r.Ratio.[run]
                      yield cellText r.Intensity.[run] ]
            writer.WriteLine(String.Join("\t", cells))
        let singleRun = results |> Array.filter (fun r -> (r.Ratio |> Array.filter (fun v -> not (Double.IsNaN v))).Length = 1) |> Array.length
        let split = results |> Array.filter (fun r -> r.RunGroups > 1) |> Array.length
        logger.Trace (sprintf "Protein groups written: %i, present in one run only: %i, with more than one run group: %i" results.Length singleRun split)
        logger.Info (sprintf "%i rows read, %i ignored, %i runs, %i protein groups written (%i in one run only, %i with more than one run group)"
                        table.RowsRead table.RowsIgnored table.Runs.Length results.Length singleRun split)
