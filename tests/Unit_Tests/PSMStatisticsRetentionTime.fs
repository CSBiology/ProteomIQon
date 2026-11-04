module ProteomIQon.UnitTestsPSMStatisticsRetentionTime

open Expecto
open ProteomIQon
open ProteomIQon.PSMStatistics
open ProteomIQon.PSMStatistics.RetentionTime

let private residues = "ACDEFGHIKLMNPQRSTVWY"

/// Random tryptic-like sequences with oxidised methionines and carbamidomethylated cysteines, written as the search results write them.
let private sequences (seed: int) n =
    let rnd = System.Random(seed)
    Array.init n (fun _ ->
        let length = rnd.Next(7, 25)
        let tokens =
            Array.init length (fun i ->
                if i = length - 1 then string "KR".[rnd.Next 2]
                else
                    let residue = residues.[rnd.Next residues.Length]
                    match residue with
                    | 'M' when rnd.Next 3 = 0 -> "[ox]M"
                    | 'C' -> "[ca]C"
                    | r -> string r)
        String.concat "" tokens)
    |> Array.distinct

/// A time that depends on the composition and the ends of the peptide, plus noise.
let private timeOf (rnd: System.Random) (sequence: string) =
    let stripped = System.Text.RegularExpressions.Regex.Replace(sequence, @"\[[^\]]*\]", "")
    let weight (c: char) = float (residues.IndexOf c) * 0.4 - 3.
    let composition = stripped |> Seq.sumBy weight
    let ends = weight stripped.[0] + weight stripped.[stripped.Length - 2]
    50. + composition + 2. * ends + (rnd.NextDouble() - 0.5) * 4.

let private strip (sequence: string) = System.Text.RegularExpressions.Regex.Replace(sequence, @"\[[^\]]*\]", "")

let private psm (id: int) (scan: int) (label: int) (score: float) (time: float) (sequence: string) : Dto.PeptideSpectrumMatchingResult =
    {
        PSMId = sprintf "psm_%i" scan; GlobalMod = 0; PepSequenceID = id; ModSequenceID = id; Label = label; ScanNr = scan; ScanTime = time
        Charge = 2; PrecursorMZ = 500.; TheoMass = 1000.; AbsDeltaMass = 0.; PeptideLength = sequence.Length; MissCleavages = 0
        SequestScore = score; SequestNormDeltaBestToRest = 0.; SequestNormDeltaNext = 0.; AndroScore = 0.; AndroNormDeltaBestToRest = 0.
        AndroNormDeltaNext = 0.; XtandemScore = 0.; XtandemNormDeltaBestToRest = 0.; XtandemNormDeltaNext = 0.; StringSequence = sequence
        IonMobility = nan; Hyperscore = nan; Expectscore = nan; MatchedIons = 0; TotalIons = 0
    }

[<Tests>]
let retentionTimeTests =
    testList "RetentionTime" [
        testCase "a target, its decoy and their modified forms share the peptide group" <| fun _ ->
            for sequence in sequences 1 500 do
                let decoy = decoySequence sequence
                Expect.equal (canonicalSequence decoy) (canonicalSequence sequence) (sprintf "%s and its decoy %s" sequence decoy)
                Expect.equal (canonicalSequence (strip sequence)) (canonicalSequence sequence) (sprintf "%s and its unmodified form" sequence)
            // sequences whose reversal without the C-terminal residue is themselves
            for modified, plain in [ "[ox]MSSMR", "MSSMR"; "[ox]MAMK", "MAMK"; "MA[ox]MK", "MAMK"; "[ca]CAAK", "CAAK" ] do
                Expect.equal (canonicalSequence modified) (canonicalSequence plain) (sprintf "%s and %s" modified plain)
                Expect.equal (canonicalSequence (decoySequence modified)) (canonicalSequence plain) (sprintf "the decoy of %s and %s" modified plain)

        testCase "the feature vectors of a target and its decoy are equal" <| fun _ ->
            let all = sequences 2 500
            let modIndex = modificationIndex all
            Expect.equal modIndex.Count 2 "both modifications are numbered"
            for sequence in all do
                let target = featureVector modIndex sequence
                let decoy = featureVector modIndex (decoySequence sequence)
                Expect.equal decoy target (sprintf "%s and its decoy" sequence)
                Expect.equal target.Length (featureWidth modIndex) "the vector has the declared width"

        testCase "the folds are stable and all used" <| fun _ ->
            let groups = sequences 5 2000 |> Array.map canonicalSequence
            let folds = groups |> Array.map (foldOfGroup 1024)
            Expect.equal folds (groups |> Array.map (foldOfGroup 1024)) "the same seed gives the same folds"
            Expect.equal (folds |> Array.distinct |> Array.sort) [| 0 .. nFolds - 1 |] "all folds are used"
            for count in folds |> Array.countBy id |> Array.map snd do
                Expect.isGreaterThan count 200 "the folds are about equal"

        testCase "a group does not train the predictor that predicts it: without signal the out-of-fold R2 is about zero" <| fun _ ->
            // Each peptide group has a random time and four modified forms, and the time has nothing to do with the sequence. A predictor that had
            // seen another form of a group would recognise the group and give a high R2.
            let rnd = System.Random(7)
            let rows =
                Array.init 400 (fun _ ->
                    let residue () = residues.[rnd.Next residues.Length]
                    let head = System.String(Array.init 4 (fun _ -> residue ()))
                    let tail = System.String(Array.init 4 (fun _ -> residue ()))
                    let time = rnd.NextDouble() * 100.
                    [| for first in [ "M"; "[ox]M" ] do
                        for second in [ "M"; "[ox]M" ] do
                            yield head + first + "A" + second + tail + "K", time |])
                |> Array.concat
            let predictor = fit 1024 rows
            Expect.isLessThan predictor.R2 0.15 "the group is not known to its predictor"

        testCase "the predictor explains a time that depends on the sequence, and sigma comes from the out-of-fold residuals" <| fun _ ->
            let rnd = System.Random(3)
            let targets = sequences 3 1500
            let rows = targets |> Array.map (fun sequence -> sequence, timeOf rnd sequence)
            let predictor = fit 1024 rows
            let predictions = predictor.Predict (rows |> Array.map fst)
            let residuals = Array.map2 (fun (_, time) (prediction: float32) -> time - float prediction) rows predictions
            let median (values: float []) = FSharp.Stats.Array.median values
            let centre = median residuals
            let sigma = 1.4826 * median (residuals |> Array.map (fun r -> abs (r - centre)))
            Expect.floatClose Accuracy.high predictor.Sigma sigma "sigma is 1.4826 times the MAD of the out-of-fold residuals"
            Expect.isGreaterThan predictor.R2 minR2 "the composition explains the time"
            Expect.equal (predictor.Predict (targets |> Array.map decoySequence)) (predictor.Predict targets) "pred(decoy) = pred(target)"

        testCase "tryInitDeviation gives each decoy the prediction of its own sequence" <| fun _ ->
            // targets have even ModSequenceIDs; the decoy of the ID n is n + 1 and holds the sequence of another target, as an ID
            // scheme in which target and decoy keys could collide
            let rnd = System.Random(11)
            let targets = sequences 9 800
            let times = targets |> Array.map (timeOf rnd)
            let n = targets.Length
            let build decoyOffset =
                [|
                    for i in 0 .. n - 1 do
                        yield psm (2 * i) i 1 (100. + rnd.NextDouble()) times.[i] targets.[i]
                    for i in 0 .. n - 1 do
                        let other = (i + 1) % n
                        yield psm (2 * i + 1 + decoyOffset) (10000 + i) -1 (rnd.NextDouble()) times.[other] targets.[other]
                |]
            let score (x: Dto.PeptideSpectrumMatchingResult) = float32 x.SequestScore
            let logger = NLog.LogManager.GetLogger "RetentionTimeTests"
            let colliding = build 0
            let separate = build 1000000
            let deviation psms = RetentionTime.tryInitDeviation logger 1024 score psms |> Option.get
            let collidingDeviation = deviation colliding
            let separateDeviation = deviation separate
            for i in 0 .. 2 * n - 1 do
                Expect.equal (collidingDeviation colliding.[i]) (separateDeviation separate.[i]) (sprintf "PSM %i" i)
            // a decoy of the sequence of a target with the time of that target has the deviation of that target
            for i in 0 .. n - 1 do
                let other = (i + 1) % n
                Expect.floatClose Accuracy.high (float (collidingDeviation colliding.[n + i])) (float (collidingDeviation colliding.[other])) "decoy and target"
    ]
