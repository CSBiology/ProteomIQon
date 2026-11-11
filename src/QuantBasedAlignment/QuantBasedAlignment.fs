namespace ProteomIQon

open System
open System.IO
open System.Collections.Generic
open FSharpAux.IO.SchemaReader
open Plotly.NET
open BioFSharp
open Dto
open Dto.QuantificationResult


/// Robust model of the shift between two runs (target minus source scan time, or ion mobility) as a function of the source value:
/// a running median over the k nearest neighbours, fitted after trimming gross outliers, with k chosen by grouped cross-validation.
module RobustShift =

    let isFinite (v: float) = not (Double.IsNaN v || Double.IsInfinity v)

    /// Median of the values, nan for an empty array.
    let median (values: float []) =
        if values.Length = 0 then
            nan
        else
            let sorted = Array.sort values
            let n = sorted.Length
            if n % 2 = 1 then sorted.[n / 2] else 0.5 * (sorted.[n / 2 - 1] + sorted.[n / 2])

    /// Percentile p (0 to 1) with linear interpolation between the order statistics, nan for an empty array.
    let percentile (p: float) (values: float []) =
        if values.Length = 0 then
            nan
        else
            let sorted = Array.sort values
            let pos = p * float (sorted.Length - 1)
            let lo = int (floor pos)
            let hi = min (lo + 1) (sorted.Length - 1)
            sorted.[lo] + (pos - float lo) * (sorted.[hi] - sorted.[lo])

    type Settings =
        {
            /// Candidate numbers of neighbours of the running median.
            KGrid            : int []
            NFolds           : int
            TrimRounds       : int
            /// Points with |residual| > max(TrimMadFactor * 1.4826 * MAD, MinTrimThreshold) are left out of the fit.
            TrimMadFactor    : float
            MinTrimThreshold : float
        }

    /// The running median at the distinct training values, linearly interpolated in between and constant beyond the ends.
    type Curve =
        {
            X     : float []
            Shift : float []
        }

    let predict (curve: Curve) (x: float) =
        let xs = curve.X
        let n = xs.Length
        if n = 0 || Double.IsNaN x then nan
        elif x <= xs.[0] then curve.Shift.[0]
        elif x >= xs.[n - 1] then curve.Shift.[n - 1]
        else
            // xs.[lo] <= x < xs.[hi]
            let hi = SortedSearch.upperBound xs x
            let lo = hi - 1
            let t = (x - xs.[lo]) / (xs.[hi] - xs.[lo])
            curve.Shift.[lo] + t * (curve.Shift.[hi] - curve.Shift.[lo])

    /// Fits the running median of k neighbours to points sorted by x. Of two equally distant neighbours the left one is taken.
    let private fitSorted (k: int) (x: float []) (shift: float []) =
        let n = x.Length
        let k = min k n
        let window = Array.zeroCreate k
        let fitted =
            Array.init n (fun i ->
                let mutable lo = i
                let mutable hi = i
                while hi - lo + 1 < k do
                    if lo = 0 then hi <- hi + 1
                    elif hi = n - 1 then lo <- lo - 1
                    elif x.[i] - x.[lo - 1] <= x.[hi + 1] - x.[i] then lo <- lo - 1
                    else hi <- hi + 1
                Array.blit shift lo window 0 k
                Array.Sort window
                if k % 2 = 1 then window.[k / 2] else 0.5 * (window.[k / 2 - 1] + window.[k / 2])
                )
        // one knot per distinct x
        let knots =
            Array.zip x fitted
            |> Array.groupBy fst
            |> Array.map (fun (xv, group) -> xv, group |> Array.averageBy snd)
        {
            X     = Array.map fst knots
            Shift = Array.map snd knots
        }

    /// Fits the curve to the points with the given indices.
    let fitCurve (k: int) (x: float []) (shift: float []) (indices: int []) =
        let sorted = indices |> Array.sortWith (fun a b -> compare (x.[a], shift.[a]) (x.[b], shift.[b]))
        fitSorted k (sorted |> Array.map (fun i -> x.[i])) (sorted |> Array.map (fun i -> shift.[i]))

    /// The shift of every point predicted by the curve fitted to the selected points of the other folds. If all selected
    /// points are in the fold of a point, by the curve fitted to the selected points without the point itself.
    let outOfFold (settings: Settings) (k: int) (x: float []) (shift: float []) (fold: int []) (selected: bool []) =
        let predicted = Array.create x.Length nan
        for f in 0 .. settings.NFolds - 1 do
            let train = [| for i in 0 .. x.Length - 1 do if selected.[i] && fold.[i] <> f then yield i |]
            if train.Length > 0 then
                let curve = fitCurve k x shift train
                for i in 0 .. x.Length - 1 do
                    if fold.[i] = f then predicted.[i] <- predict curve x.[i]
            else
                for i in 0 .. x.Length - 1 do
                    if fold.[i] = f then
                        let others = [| for j in 0 .. x.Length - 1 do if selected.[j] && j <> i then yield j |]
                        if others.Length > 0 then predicted.[i] <- predict (fitCurve k x shift others) x.[i]
        predicted

    /// Median absolute out-of-fold error of the selected points.
    let private cvError (settings: Settings) (k: int) (x: float []) (shift: float []) (fold: int []) (selected: bool []) =
        let predicted = outOfFold settings k x shift fold selected
        let errors = [| for i in 0 .. x.Length - 1 do if selected.[i] && isFinite predicted.[i] then yield abs (shift.[i] - predicted.[i]) |]
        let m = median errors
        if isFinite m then m else infinity

    /// The k of the grid with the smallest median absolute out-of-fold error, the smallest k on ties.
    let chooseK (settings: Settings) (x: float []) (shift: float []) (fold: int []) (selected: bool []) =
        settings.KGrid
        |> Array.map (fun k -> k, cvError settings k x shift fold selected)
        |> Array.minBy snd
        |> fst

    type Fit =
        {
            Curve     : Curve
            K         : int
            /// The points that survived the trimming; the curve is fitted to them.
            Kept      : bool []
            /// Out-of-fold shift of every point, trimmed points included.
            OutOfFold : float []
            /// Absolute out-of-fold error of every point.
            AbsError  : float []
        }

    /// Robust fit: k by cross-validation on all points, up to TrimRounds rounds of trimming with that k,
    /// then k by cross-validation on the kept points and the final curve on the kept points.
    let fit (settings: Settings) (x: float []) (shift: float []) (fold: int []) =
        let n = x.Length
        let indicesOf (mask: bool []) = Array.init n id |> Array.filter (fun i -> mask.[i])
        let all = Array.create n true
        let k0 = chooseK settings x shift fold all
        let minKept = max settings.KGrid.[0] settings.NFolds
        let rec trim round (kept: bool []) =
            if round = settings.TrimRounds then
                kept
            else
                let curve = fitCurve k0 x shift (indicesOf kept)
                let residual = Array.init n (fun i -> shift.[i] - predict curve x.[i])
                let keptResidual = indicesOf kept |> Array.map (fun i -> residual.[i])
                let center = median keptResidual
                let mad = 1.4826 * median (keptResidual |> Array.map (fun r -> abs (r - center)))
                let threshold = max (settings.TrimMadFactor * mad) settings.MinTrimThreshold
                let next = residual |> Array.map (fun r -> abs r <= threshold)
                if next = kept || (next |> Array.filter id |> Array.length) < minKept then kept else trim (round + 1) next
        let kept = trim 0 all
        let k = chooseK settings x shift fold kept
        let outOfFoldShift = outOfFold settings k x shift fold kept
        {
            Curve     = fitCurve k x shift (indicesOf kept)
            K         = k
            Kept      = kept
            OutOfFold = outOfFoldShift
            AbsError  = Array.init n (fun i -> abs (shift.[i] - outOfFoldShift.[i]))
        }

    /// Robust local error scale: 1.4826 * median absolute out-of-fold error of the points within +-HalfWindow of x.
    type ErrorScale =
        {
            X          : float []
            AbsError   : float []
            /// Used where fewer than MinPoints points are in the window.
            Global     : float
            HalfWindow : float
            MinPoints  : int
            Floor      : float
        }

    let initErrorScale (halfWindow: float) (minPoints: int) (floor: float) (x: float []) (absError: float []) =
        let points =
            Array.zip x absError
            |> Array.filter (fun (a, b) -> isFinite a && isFinite b)
            |> Array.sortBy fst
        {
            X          = Array.map fst points
            AbsError   = Array.map snd points
            Global     = 1.4826 * median (Array.map snd points)
            HalfWindow = halfWindow
            MinPoints  = minPoints
            Floor      = floor
        }

    /// The local scale at x, with the absolute error ownError of a point of the scale left out of the window (nan: none).
    let localScaleWithout (scale: ErrorScale) (ownError: float) (x: float) =
        let local =
            if not (isFinite x) then
                scale.Global
            else
                let lo = SortedSearch.lowerBound scale.X (x - scale.HalfWindow)
                let hi = SortedSearch.upperBound scale.X (x + scale.HalfWindow)
                let window = scale.AbsError.[lo .. hi - 1]
                let window =
                    match Array.tryFindIndex (fun e -> e = ownError) window with
                    | Some j -> Array.append window.[.. j - 1] window.[j + 1 ..]
                    | None -> window
                if window.Length >= scale.MinPoints then 1.4826 * median window else scale.Global
        if isFinite local then max scale.Floor local else scale.Floor


module QuantBasedAlignment =

    open RobustShift

    ///
    type PeptideIon =
        {
            Sequence             : string
            GlobalMod            : int
            Charge               : int
        }

    ///
    let toPeptideIon (qp:QuantificationResult) =
        {
            Sequence             = qp.StringSequence
            GlobalMod            = qp.GlobalMod
            Charge               = qp.Charge
        }

    /// Scan time model (minutes).
    let scanTimeSettings =
        {
            KGrid            = [|5; 7; 11; 15; 21; 31; 41; 51; 71; 101|]
            NFolds           = 5
            TrimRounds       = 3
            TrimMadFactor    = 3.
            MinTrimThreshold = 0.05
        }

    /// Ion mobility model (1/K0).
    let ionMobilitySettings = { scanTimeSettings with MinTrimThreshold = 0.005 }

    /// Source and target pairs with fewer shared ions are not modelled.
    let minSharedIons = 20

    /// The local error scale of a scan time is taken from the shared ions within +-5 minutes, and is at least 0.01 minutes.
    let initScanTimeScale = initErrorScale 5. 10 0.01

    /// The local error scale of an ion mobility is taken from the shared ions within +-0.05 1/K0, and is at least 0.001.
    let initIonMobilityScale = initErrorScale 0.05 10 0.001

    /// A row of a .quant file that passes the filter with the scan time of the peak apex of its target label.
    type QuantIon =
        {
            Row       : QuantificationResult
            Apex      : float
        }

    ///
    type QuantFile =
        {
            FileName  : string
            /// The ions that pass the filter; of several rows of an ion the last one is kept.
            Ions      : Dictionary<PeptideIon,QuantIon>
            /// The keys of Ions in file order.
            IonOrder  : PeptideIon []
        }

    /// The scan time of the peak apex of the target label of an ion that is used, else the reason why it is not used.
    let apexOrDropReason (qp:QuantificationResult) =
        if qp.GlobalMod <> 0 && qp.GlobalMod <> 1 then
            Error "GlobalMod is neither 0 nor 1"
        else
            let peakParams = if qp.GlobalMod = 0 then qp.Params_Light else qp.Params_Heavy
            if Array.isEmpty peakParams then
                Error "no fitted peak"
            elif peakParams.Length < 3 then
                Error "fit without scan time and sd"
            else
                let scanTime = PeakShape.apex peakParams
                let sd = PeakShape.width peakParams
                let difference = getTargetApexDifferenceAt scanTime qp
                if not (isFinite scanTime && isFinite sd && sd > 0. && isFinite difference) then
                    Error "non-finite fitted scan time, sd or search to fitted scan time difference"
                // the fit is a good estimate of the scan time
                elif abs difference / sd >= 2. then
                    Error "search scan time 2 or more sd away from the fitted scan time"
                else
                    Ok scanTime

    /// Reads a .quant file and keeps the ions whose fitted scan time can be trusted, as QuantBasedAlignment does.
    let readQuantFile (quantFilePath:string) =
        let fileName = Path.GetFileNameWithoutExtension quantFilePath
        let logger = Logging.createLogger fileName
        let rows =
            Csv.CsvReader<QuantificationResult>().ReadFile(quantFilePath,'\t',false,1)
            |> Array.ofSeq
        let ions = Dictionary<PeptideIon,QuantIon>(HashIdentity.Structural)
        let ionOrder = ResizeArray<PeptideIon>()
        let dropCounts = SortedDictionary<string,int>()
        let mutable duplicates = 0
        for qp in rows do
            match apexOrDropReason qp with
            | Error reason ->
                logger.Debug (sprintf "Ion not used, %s: %s, GlobalMod %i, charge %i" reason qp.StringSequence qp.GlobalMod qp.Charge)
                dropCounts.[reason] <- (match dropCounts.TryGetValue reason with | true, c -> c + 1 | _ -> 1)
            | Ok apex ->
                let ion = toPeptideIon qp
                if ions.ContainsKey ion then
                    logger.Debug (sprintf "Ion has several rows, the last one is used: %s, GlobalMod %i, charge %i" qp.StringSequence qp.GlobalMod qp.Charge)
                    duplicates <- duplicates + 1
                else
                    ionOrder.Add ion
                ions.[ion] <- { Row = qp; Apex = apex }
        let drops =
            dropCounts
            |> Seq.map (fun kv -> sprintf "%s: %i" kv.Key kv.Value)
            |> String.concat ", "
        logger.Info (sprintf "%s: %i rows, %i ions used; rows not used: %s; repeated rows of used ions (last one kept): %i" fileName rows.Length ions.Count (if drops = "" then "none" else drops) duplicates)
        {
            FileName  = fileName
            Ions      = ions
            IonOrder  = ionOrder.ToArray()
        }

    /// Assigns every peptide (PepSequenceID, so all modified forms of a sequence) to one of nFolds cross-validation folds by a
    /// seeded permutation of the sorted ids.
    let assignFolds (seed:int) (nFolds:int) (peptides: seq<int>) =
        let sorted =
            peptides
            |> Seq.distinct
            |> Array.ofSeq
            |> Array.sort
        let rnd = Random(seed)
        for i in sorted.Length - 1 .. -1 .. 1 do
            let j = rnd.Next(i + 1)
            let tmp = sorted.[i]
            sorted.[i] <- sorted.[j]
            sorted.[j] <- tmp
        let foldOf = Dictionary<int,int>()
        sorted |> Array.iteri (fun i peptide -> foldOf.[peptide] <- i % nFolds)
        foldOf

    ///
    type IonMobilityModel =
        {
            Fit       : Fit
            /// Source ion mobility of the ions the model was trained on.
            X         : float []
            /// Shift of the ions the model was trained on.
            Shift     : float []
            /// Out-of-fold predicted target ion mobility of the ions the model was trained on.
            OutOfFold : Dictionary<PeptideIon,float>
            /// Absolute out-of-fold error of the ions the model was trained on.
            AbsError  : Dictionary<PeptideIon,float>
            Scale     : ErrorScale
        }

    ///
    type PairModel =
        {
            Source     : QuantFile
            SharedIons : PeptideIon []
            Fold       : int []
            /// Source scan time of the shared ions.
            X          : float []
            /// Target minus source scan time of the shared ions.
            Shift      : float []
            ScanTime   : Fit
            /// Out-of-fold predicted target scan time of the shared ions.
            OutOfFold  : Dictionary<PeptideIon,float>
            /// Absolute out-of-fold scan time error of the shared ions.
            AbsError   : Dictionary<PeptideIon,float>
            Scale      : ErrorScale
            /// Median absolute out-of-fold scan time error over all shared ions (minutes); ranks the sources.
            ModelError : float
            IonMobility: IonMobilityModel option
        }

    /// Models the target scan time (and ion mobility) as a function of the source scan time (and ion mobility) on the ions both runs share.
    let alignPair (logger:NLog.Logger) (foldOf:Dictionary<int,int>) (target:QuantFile) (source:QuantFile) =
        let pairName = sprintf "%s -> %s" source.FileName target.FileName
        let shared = target.IonOrder |> Array.filter source.Ions.ContainsKey
        if shared.Length < minSharedIons then
            logger.Warn (sprintf "%s: %i shared ions, at least %i are needed; the source is not used." pairName shared.Length minSharedIons)
            None
        else
        let fold  = shared |> Array.map (fun ion -> foldOf.[target.Ions.[ion].Row.PepSequenceID])
        let x     = shared |> Array.map (fun ion -> source.Ions.[ion].Apex)
        let shift = shared |> Array.mapi (fun i ion -> target.Ions.[ion].Apex - x.[i])
        let scanTime = fit scanTimeSettings x shift fold
        let outOfFold = Dictionary<PeptideIon,float>(HashIdentity.Structural)
        shared |> Array.iteri (fun i ion -> if isFinite scanTime.OutOfFold.[i] then outOfFold.[ion] <- x.[i] + scanTime.OutOfFold.[i])
        let absError = Dictionary<PeptideIon,float>(HashIdentity.Structural)
        shared |> Array.iteri (fun i ion -> absError.[ion] <- scanTime.AbsError.[i])
        let errors = scanTime.AbsError |> Array.filter isFinite
        let modelError = median errors
        let trimmed = scanTime.Kept |> Array.filter not |> Array.length
        logger.Info (sprintf "%s: %i shared ions, %i trimmed, k = %i, out-of-fold absolute scan time error median %.4f / p90 %.4f min" pairName shared.Length trimmed scanTime.K modelError (percentile 0.9 errors))
        let ionMobility =
            let withIonMobility =
                [|
                    for i in 0 .. shared.Length - 1 do
                        if isFinite source.Ions.[shared.[i]].Row.IonMobility && isFinite target.Ions.[shared.[i]].Row.IonMobility then yield i
                |]
            if withIonMobility.Length = 0 then
                logger.Trace (sprintf "%s: no ion mobility." pairName)
                None
            elif withIonMobility.Length < minSharedIons then
                logger.Warn (sprintf "%s: %i shared ions with ion mobility in both runs, at least %i are needed; no ion mobility model." pairName withIonMobility.Length minSharedIons)
                None
            else
                let ionsIM  = withIonMobility |> Array.map (fun i -> shared.[i])
                let xIM     = ionsIM |> Array.map (fun ion -> source.Ions.[ion].Row.IonMobility)
                let shiftIM = ionsIM |> Array.mapi (fun i ion -> target.Ions.[ion].Row.IonMobility - xIM.[i])
                let foldIM  = withIonMobility |> Array.map (fun i -> fold.[i])
                let fitIM = fit ionMobilitySettings xIM shiftIM foldIM
                let outOfFoldIM = Dictionary<PeptideIon,float>(HashIdentity.Structural)
                ionsIM |> Array.iteri (fun i ion -> if isFinite fitIM.OutOfFold.[i] then outOfFoldIM.[ion] <- xIM.[i] + fitIM.OutOfFold.[i])
                let absErrorIM = Dictionary<PeptideIon,float>(HashIdentity.Structural)
                ionsIM |> Array.iteri (fun i ion -> absErrorIM.[ion] <- fitIM.AbsError.[i])
                let errorsIM = fitIM.AbsError |> Array.filter isFinite
                logger.Info (
                    sprintf "%s: ion mobility model on %i shared ions, %i trimmed, k = %i, out-of-fold absolute ion mobility error median %.5f / p90 %.5f"
                        pairName ionsIM.Length (fitIM.Kept |> Array.filter not |> Array.length) fitIM.K (median errorsIM) (percentile 0.9 errorsIM)
                    )
                Some {
                    Fit       = fitIM
                    X         = xIM
                    Shift     = shiftIM
                    OutOfFold = outOfFoldIM
                    AbsError  = absErrorIM
                    Scale     = initIonMobilityScale xIM fitIM.AbsError
                }
        Some {
            Source      = source
            SharedIons  = shared
            Fold        = fold
            X           = x
            Shift       = shift
            ScanTime    = scanTime
            OutOfFold   = outOfFold
            AbsError    = absError
            Scale       = initScanTimeScale x scanTime.AbsError
            ModelError  = modelError
            IonMobility = ionMobility
        }

    /// Transfers an ion of the source of the pair. Ions the model was trained on get their out-of-fold prediction and
    /// a local error scale without their own error, so that they show the same error as ions the target lacks.
    let createAlignmentResult (pair:PairModel) (ion:PeptideIon) : AlignmentResult =
        let quantifiedPeptide = pair.Source.Ions.[ion].Row
        let ownError (absError: Dictionary<PeptideIon,float>) = match absError.TryGetValue ion with | true, e -> e | _ -> nan
        let sourceScanTime = pair.Source.Ions.[ion].Apex
        let predictedScanTime =
            match pair.OutOfFold.TryGetValue ion with
            | true, v -> v
            | _ -> sourceScanTime + predict pair.ScanTime.Curve sourceScanTime
        let predictedIonMobility, predictedIonMobilitySd =
            match pair.IonMobility with
            | Some model when isFinite quantifiedPeptide.IonMobility ->
                let sourceIonMobility = quantifiedPeptide.IonMobility
                let predicted =
                    match model.OutOfFold.TryGetValue ion with
                    | true, v -> v
                    | _ -> sourceIonMobility + predict model.Fit.Curve sourceIonMobility
                predicted, localScaleWithout model.Scale (ownError model.AbsError) sourceIonMobility
            | _ -> nan, nan
        {
            StringSequence                                  = quantifiedPeptide.StringSequence
            GlobalMod                                       = quantifiedPeptide.GlobalMod
            Charge                                          = quantifiedPeptide.Charge
            PepSequenceID                                   = quantifiedPeptide.PepSequenceID
            ModSequenceID                                   = quantifiedPeptide.ModSequenceID
            Mz                                              = Mass.toMZ (quantifiedPeptide.TheoMass) (float quantifiedPeptide.Charge)
            ProteinNames                                    = quantifiedPeptide.ProteinNames
            PredictedScanTime                               = predictedScanTime
            ScanTime_SourceFile                             = sourceScanTime
            ApexIntensity_SourceFile                        =
                if quantifiedPeptide.GlobalMod = 1 then
                    quantifiedPeptide.MeasuredApex_Heavy
                else
                    quantifiedPeptide.MeasuredApex_Light
            Quant_SourceFile                                = getTargetIntensity quantifiedPeptide
            RtTrace_SourceFile                              = getTargetRtTrace quantifiedPeptide
            IntensityTrace_SourceFile                       = getTargetIntensityTrace quantifiedPeptide
            IsotopicPatternMz_SourceFile                    = getIsotopicPatternMz quantifiedPeptide
            IsotopicPatternIntensity_Observed_SourceFile    = getIsotopicPatternIntensity_Observed quantifiedPeptide
            IonMobility                                     = quantifiedPeptide.IonMobility
            PredictedScanTimeSd                             = localScaleWithout pair.Scale (ownError pair.AbsError) sourceScanTime
            SourceFileName                                  = pair.Source.FileName
            SourceModelError                                = pair.ModelError
            PredictedIonMobility                            = predictedIonMobility
            PredictedIonMobilitySd                          = predictedIonMobilitySd
        }

    /// The out-of-fold predictions of the shared ions in the first cross-validation fold, in the .alignmetric layout.
    let createMetrics (target:QuantFile) (pair:PairModel) : AlignmentMetricsDTO [] =
        pair.SharedIons
        |> Array.mapi (fun i ion -> i, ion)
        |> Array.filter (fun (i, ion) -> pair.Fold.[i] = 0 && pair.OutOfFold.ContainsKey ion)
        |> Array.map (fun (i, ion) ->
            let sourcePep = pair.Source.Ions.[ion].Row
            let targetPep = target.Ions.[ion].Row
            let predicted = pair.OutOfFold.[ion]
            {
                Sequence                             = ion.Sequence
                GlobalMod                            = ion.GlobalMod
                Charge                               = ion.Charge
                PepSequenceID                        = sourcePep.PepSequenceID
                ModSequenceID                        = sourcePep.ModSequenceID
                X_FileName                           = pair.Source.FileName
                X_Intensities                        = getTargetIntensity sourcePep
                X_Stabw                              = getTargetPeakWidth sourcePep
                X_Test                               = pair.X.[i]
                X_IsotopicPatternMz                  = getIsotopicPatternMz sourcePep
                X_IsotopicPatternIntensity_Observed  = getIsotopicPatternIntensity_Observed sourcePep
                X_RtTrace                            = getTargetRtTrace sourcePep
                X_IntensityTrace                     = getTargetIntensityTrace sourcePep
                Y_Test                               = target.Ions.[ion].Apex
                YHat_Test                            = predicted
                YHat_Refined_Test                    = predicted
                Y_Intensity                          = getTargetIntensity targetPep
                Y_IsotopicPatternMz                  = getIsotopicPatternMz targetPep
                Y_IsotopicPatternIntensity_Observed  = getIsotopicPatternIntensity_Observed targetPep
                Y_RtTrace                            = getTargetRtTrace targetPep
                Y_IntensityTrace                     = getTargetIntensityTrace targetPep
                DtwDistanceBefore                    = nan
                IonMobility                          = sourcePep.IonMobility
            }
            )

    /// Per pair: the shifts of the shared ions with the fitted curve, and the out-of-fold errors.
    let createModelChart (target:QuantFile) (pairs:PairModel []) =
        pairs
        |> Array.collect (fun pair ->
            let kept, trimmed = Array.init pair.X.Length id |> Array.partition (fun i -> pair.ScanTime.Kept.[i])
            let points (indices:int []) (values:float []) = indices |> Array.map (fun i -> pair.X.[i], values.[i])
            let shiftChart =
                [
                    Chart.Point(points kept pair.Shift)    |> Chart.withTraceInfo (sprintf "%s used" pair.Source.FileName)
                    Chart.Point(points trimmed pair.Shift) |> Chart.withTraceInfo (sprintf "%s trimmed" pair.Source.FileName)
                    Chart.Line(Array.zip pair.ScanTime.Curve.X pair.ScanTime.Curve.Shift) |> Chart.withTraceInfo (sprintf "%s running median, k = %i" pair.Source.FileName pair.ScanTime.K)
                ]
                |> Chart.combine
                |> Chart.withXAxisStyle "source scan time"
                |> Chart.withYAxisStyle "target - source scan time"
            let errorChart =
                Chart.Point(points (Array.append kept trimmed) (Array.map2 (-) pair.Shift pair.ScanTime.OutOfFold))
                |> Chart.withTraceInfo (sprintf "%s out-of-fold error" pair.Source.FileName)
                |> Chart.withXAxisStyle "source scan time"
                |> Chart.withYAxisStyle "target scan time - out-of-fold prediction"
            [|shiftChart; errorChart|]
            )
        |> Chart.Grid(pairs.Length, 2)
        |> Chart.withSize(2000, 600 * pairs.Length)
        |> Chart.withTitle target.FileName

    /// Aligns all sources to the target and writes the .align and .alignmetric files of the target.
    let alignFiles diagCharts (outputDir:string) (foldOf:Dictionary<int,int>) (sources:QuantFile []) (target:QuantFile) =
        let logger = Logging.createLogger target.FileName
        logger.Info (sprintf "%s: aligning %i source runs" target.FileName sources.Length)
        let pairs =
            sources
            |> Array.choose (alignPair logger foldOf target)
            |> Array.sortBy (fun pair -> pair.ModelError, pair.Source.FileName)
        if pairs.Length = 0 then
            logger.Warn (sprintf "%s: no source run shares at least %i ions with it; no ion is transferred." target.FileName minSharedIons)
        logger.Info (
            pairs
            |> Array.map (fun pair -> sprintf "%s (%.4f min)" pair.Source.FileName pair.ModelError)
            |> String.concat ", "
            |> sprintf "%s: sources ranked by median absolute out-of-fold scan time error: %s" target.FileName
            )
        if diagCharts && pairs.Length > 0 then
            createModelChart target pairs
            |> Chart.saveHtml(Path.Combine [|outputDir; target.FileName + "_ScanTimeModels"|])
        // every ion of the target and the sources, from the best ranked source that has it
        let ionsToAlign =
            Seq.append target.IonOrder (sources |> Seq.collect (fun source -> source.IonOrder))
            |> Seq.distinct
            |> Array.ofSeq
        let transferred =
            ionsToAlign
            |> Array.choose (fun ion ->
                pairs
                |> Array.tryFind (fun pair -> pair.Source.Ions.ContainsKey ion)
                |> Option.map (fun pair -> ion, createAlignmentResult pair ion)
                )
        let targetOnly = ionsToAlign |> Array.filter (fun ion -> sources |> Array.forall (fun source -> source.Ions.ContainsKey ion |> not)) |> Array.length
        let unusedSourcesOnly = ionsToAlign.Length - transferred.Length - targetOnly
        let overlap = transferred |> Array.filter (fun (ion, _) -> target.Ions.ContainsKey ion) |> Array.length
        let results = transferred |> Array.map snd
        logger.Info (
            sprintf "%s: %i ions transferred (%i of them used in the target, %i new), per source: %s; median PredictedScanTimeSd %.4f min; %i with predicted ion mobility; not transferred: %i quantified in the target only, %i only in sources that are not used"
                target.FileName results.Length overlap (results.Length - overlap)
                (pairs |> Array.map (fun pair -> sprintf "%s %i" pair.Source.FileName (results |> Array.filter (fun r -> r.SourceFileName = pair.Source.FileName) |> Array.length)) |> String.concat ", ")
                (results |> Array.map (fun r -> r.PredictedScanTimeSd) |> median)
                (results |> Array.filter (fun r -> isFinite r.PredictedIonMobility) |> Array.length)
                targetOnly unusedSourcesOnly
            )
        let alignFilePath = Path.Combine [|outputDir; target.FileName + ".align"|]
        File.WriteAllLines(alignFilePath, results |> SeqIO'.csv "\t" true false)
        let metricFilePath = Path.Combine [|outputDir; target.FileName + ".alignmetric"|]
        File.WriteAllLines(metricFilePath, pairs |> Array.collect (createMetrics target) |> SeqIO'.csv "\t" true false)
        logger.Trace (sprintf "Results written to %s and %s" alignFilePath metricFilePath)

    /// Reads the target and source .quant files, assigns the cross-validation folds by peptide and aligns each target
    /// file to the source files other than itself, writing its .align and .alignmetric files to outputDir. c is the
    /// number of files processed in parallel.
    let alignRuns (diagCharts: bool) (outputDir: string) (seed: int) (c: int) (targetFiles: string []) (sourceFiles: string []) =
        let quantFiles =
            Array.append targetFiles sourceFiles
            |> Array.distinct
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> FSharpAux.PSeq.map (fun path -> path, readQuantFile path)
            |> dict
        let foldOf =
            quantFiles.Values
            |> Seq.collect (fun quantFile -> quantFile.IonOrder |> Seq.map (fun ion -> quantFile.Ions.[ion].Row.PepSequenceID))
            |> assignFolds seed scanTimeSettings.NFolds
        targetFiles
        |> FSharpAux.PSeq.withDegreeOfParallelism c
        |> FSharpAux.PSeq.iter (fun targetFile ->
            let sources =
                sourceFiles
                |> Array.filter (fun sourceFile -> sourceFile <> targetFile)
                |> Array.map (fun sourceFile -> quantFiles.[sourceFile])
            alignFiles diagCharts outputDir foldOf sources quantFiles.[targetFile]
            )
