namespace ProteomIQon

open System
open System.IO
open System.Data.SQLite
open ProteomIQon.Core
open Core.MzIO
open Dto
open FSharp.Stats
open BioFSharp.Mz.Quantification
open BioFSharp.Mz
open FSharpAux.IO.SchemaReader
open Plotly.NET
open BioFSharp
open MzIO.Processing
open BioFSharp.Mz.SearchDB
open Gabor3D.algorithm
open PeakQuantification
open PSMQuantification

module PSMBasedQuantificationTIMs =
    module Query = 
        open System
        open System.Linq
        open MzIO.Commons.Arrays
        open MzIO.Processing.MzIOLinq
        open MzIO.Binary
         
        /// Reads the peaks of a spectrum sorted by m/z.
        let readSortedPeaks (reader: MzIO.MzSQL.MzSQL) (spectrumID: string) =
            let peaks = Xic.readSortedPeaks reader spectrumID
            // a peak without ion mobility gets NaN and falls outside every ion mobility window
            if peaks.IonMobility.Length = 0 then { peaks with IonMobility = Array.create peaks.Mz.Length nan } else peaks

        // extracts retention time profile for given mz and rt range and ion mobility range. 
        let initRTProfileGabor(readSpecPeaks: string -> Peak1DArray)(rtIndex: IMzIOArray<RtIndexEntry>)(rtRange: RangeQuery)(mzRange: RangeQuery)(ionMobilityRange: RangeQuery) =
            let entries =
                RtIndexEntry.Search(rtIndex, rtRange).ToArray()
            [|
                for entry in entries do
                    let peaks =
                        (readSpecPeaks entry.SpectrumID).Peaks
                    yield
                        RtIndexEntry.MzSearch(peaks, mzRange)
                        |> Seq.choose (fun peak ->
                            match peak.IonMobility with
                            | Some mobility
                                when
                                    mobility > ionMobilityRange.LowValue
                                    && mobility < ionMobilityRange.HighValue
                                ->

                                Some (
                                    RtIndexEntry.AsPeak2D(
                                        peak,
                                        entry.Rt
                                    )
                                )

                            | _ ->
                                None
                        )
                        |> Seq.toArray
                    
            |]

    /// Reading of the .qpsm of PSMStatistics and of the ion mobility of the PSM spectra.
    module QpsmInput =

        /// Column names of the .qpsm in the order the positional reader expects them.
        let expectedColumns =
            Reflection.FSharpType.GetRecordFields(typeof<PSMStatisticsResult>)
            |> Array.map (fun field -> field.Name)

        /// Fails when the header of the file is not the column layout of PSMStatisticsResult,
        /// for example a .qpsm in the former FragPipe layout of MsFraggerToPSM.
        let checkHeader (path: string) =
            let header =
                use reader = new StreamReader(path)
                match reader.ReadLine() with
                | null -> [||]
                | line -> line.TrimEnd('\r').Split('\t')
            if header <> expectedColumns then
                let missing = expectedColumns |> Array.filter (fun c -> not (Array.contains c header))
                let unexpected = header |> Array.filter (fun c -> not (Array.contains c expectedColumns))
                failwithf
                    "%s does not have the column layout of the PSMStatistics .qpsm. Missing columns: [%s]. Unexpected columns: [%s]. Expected in this order: [%s]."
                    path
                    (String.Join(", ", missing))
                    (String.Join(", ", unexpected))
                    (String.Join(", ", expectedColumns))

        /// Reads the .qpsm after checking its header. An empty file has no PSMs.
        let read (path: string) =
            if FileInfo(path).Length = 0L then
                [||]
            else
                checkHeader path
                Csv.CsvReader<PSMStatisticsResult>(SchemaMode=Csv.Fill).ReadFile(path,'\t',false,1)
                |> Array.ofSeq

        /// The retention time index and the scan times of the MS1 spectra, both sorted by scan time (SpectrumIndex.readMs1Index),
        /// and the inverse reduced ion mobility (MS:1002815, NaN when absent) of the spectra in psmSpectrumIds, from one
        /// pass over the spectrum descriptions with the streaming parser.
        let readSpectrumIndex (reader: MzIO.MzSQL.MzSQL) (runID: string) (warn: string -> unit) (psmSpectrumIds: Collections.Generic.HashSet<string>) =
            let mobilities = Collections.Generic.Dictionary<string,float>()
            let ms1SortedByScanTime =
                SpectrumIndex.readMs1Index reader runID warn (fun header ->
                    if psmSpectrumIds.Contains header.ID then
                        mobilities.[header.ID] <- header.IonMobility |> Option.defaultValue nan)
            let rtIndex =
                ms1SortedByScanTime
                |> Array.map (fun (scanTime, id) -> MzIO.Processing.MzIOLinq.RtIndexEntry(scanTime, id))
                |> MzIO.Commons.Arrays.MzIOArray.ToMzIOArray
            rtIndex, ms1SortedByScanTime, mobilities

    ///
    let getSpec (reader:MzIO.IO.IMzIODataReader) (ms1ID: string)  =
        Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(ms1ID).Peaks)
        |> fun (mzData,intensityData) -> PeakArray.zip mzData intensityData

    ///
    let initGetProcessedXIC (baseLineCorrection:Domain.BaseLineCorrection option) readSortedPeaks idx scanTimeWindow mzWindow_Da imWindow meanScanTime meanPrecMz meanIM =
        let rtQuery = Query.createRangeQuery meanScanTime scanTimeWindow
        let rt, intensity = Xic.rtProfile readSortedPeaks idx rtQuery meanPrecMz mzWindow_Da meanIM imWindow
        Xic.processXic baseLineCorrection rt intensity

    /// The XIC of an ion located at its PSMs: scan time, m/z and expected ion mobility (the one of the PSMs) to the ion mobility
    /// centre of the ion (None if there is no signal to take it from) and the XIC summed within +-halfWidth around it.
    /// The raw XIC is extracted once over +-halfWidth around the expected ion mobility and the pairs of the mobility profile are
    /// recorded. The centre is Xic.tryMobilityCentre of the raw XIC; with a centre the quantified XIC is the profile summed around
    /// it (Xic.summedTrace), without one it is the raw XIC at the expected ion mobility; both go through Xic.processXic. The profile
    /// belongs to the returned function and is emptied for every ion, so it serves the ions of a file, which are quantified one after
    /// the other. An expected ion mobility that is not finite gives the XIC at it (initGetProcessedXIC).
    let initGetTargetXic (baseLineCorrection: Domain.BaseLineCorrection option) readSortedPeaks idx scanTimeWindow mzWindow_Da (halfWidth: float) =
        let profile = Xic.MobilityProfile()
        fun (scanTime: float) (mz: float) (expectedMobility: float) ->
            if not (Numerics.isFinite expectedMobility) then
                None, initGetProcessedXIC baseLineCorrection readSortedPeaks idx scanTimeWindow mzWindow_Da halfWidth scanTime mz expectedMobility
            else
                let rtQuery = Query.createRangeQuery scanTime scanTimeWindow
                let rt, intensity = Xic.rtProfileWithMobilityProfile readSortedPeaks idx rtQuery mz mzWindow_Da expectedMobility halfWidth profile
                match Xic.tryMobilityCentre profile rt intensity scanTime expectedMobility halfWidth with
                | Some centre -> Some centre, Xic.processXic baseLineCorrection rt (fst (Xic.summedTrace profile centre halfWidth))
                | None -> None, Xic.processXic baseLineCorrection rt intensity


    [<Literal>]
    let rtBinWidth = 0.01

    [<Literal>]
    let mobilityBinWidth = 0.0005

    [<Literal>]
    let gaborRtSearchRadiusBins = 50

    [<Literal>]
    let gaborMobilitySearchRadiusBins = 20

    [<Literal>]
    let gaussianFitRadiusInSigma = 3.0

    [<Literal>]
    let gaussianFitMaximumIterations = 250

    /// Scale factors used to build a compact Gabor filter bank around the configured sigma values.
    let private gaborSigmaScaleFactors =
        [| 0.75; 1.00; 1.25 |]

    type GaborObservation =
        {
            Rt: float
            Mobility: float
            Intensity: float
        }

    /// A regular RT-mobility grid. Array indices are stored separately from physical coordinates.
    type GaborGrid =
        {
            RtOrigin: float
            MobilityOrigin: float
            RetentionTimes: float[]
            Mobilities: float[]
            RtBinCount: int
            MobilityBinCount: int
            Intensities: float[][]
        }

    type GaborGridPosition =
        {
            RtBin: int
            MobilityBin: int
        }

    /// Stores the kernel, score, and physical position belonging to one filter-bank match.
    type GaborKernelMatch =
        {
            Parameters: Domain.Gabor3DParams
            Position: GaborGridPosition
            RetentionTime: float
            IonMobility: float
            CorrectedIntensityAtPosition: float
            Magnitude: float
        }

    type GaborFilterBankResult =
        {
            BestMatch: GaborKernelMatch
            SecondBestMatch: GaborKernelMatch option
            ScoreRatio: float option
        }

    /// Contains the baseline surface assembled from independent RT-wise AsLS corrections.
    type BaselineCorrection2DResult =
        {
            Baseline: float[][]
            CorrectedIntensities: float[][]
        }

    type Gaussian2DObservation =
        {
            RetentionTime: float
            IonMobility: float
            Intensity: float
        }

    /// Parameters of an axis-aligned two-dimensional Gaussian peak model.
    type Gaussian2DParameters =
        {
            Height: float
            MeanRt: float
            MeanMobility: float
            SigmaRt: float
            SigmaMobility: float
        }

    type GaborFitWindow =
        {
            MinRtBin: int
            MaxRtBin: int
            MinMobilityBin: int
            MaxMobilityBin: int
        }

    /// Final Gabor quantification result. The old connected-region representation is intentionally removed.
    type Gabor2DQuantification =
        {
            Volume: float
            DirectWindowVolume: float
            ApexIntensity: float
            ApexRetentionTime: float
            ApexIonMobility: float
            WeightedRetentionTime: float
            WeightedIonMobility: float
            ProjectedRetentionTimes: float[]
            ProjectedIntensities: float[]
            BestGaborMatch: GaborKernelMatch
            SecondBestGaborMatch: GaborKernelMatch option
            GaborScoreRatio: float option
            InitialGaussianParameters: Gaussian2DParameters
            FittedGaussianParameters: Gaussian2DParameters
            FitRmse: float
        }

    type GaborTargetEstimate =
        {
            ExpectedRt: float
            ExpectedMz: float
            ExpectedMobility: float
        }

    let private isUsableIntensity value =
        Double.IsFinite value && value > 0.0

    let private tryWeightedMeanOrMean (weights: float[]) (values: float[]) =
        if Array.isEmpty values || weights.Length <> values.Length then
            None
        else
            let usablePairs =
                Array.zip weights values
                |> Array.filter (fun (weight, value) ->
                    Double.IsFinite weight
                    && Double.IsFinite value
                    && weight >= 0.0
                )

            if Array.isEmpty usablePairs then
                None
            else
                let usableWeights, usableValues =
                    Array.unzip usablePairs

                let weightSum =
                    Array.sum usableWeights

                let result =
                    if Double.IsFinite weightSum && weightSum > 0.0 then
                        Array.map2
                            (fun weight value -> weight * value)
                            usableWeights
                            usableValues
                        |> Array.sum
                        |> fun weightedSum -> weightedSum / weightSum
                    else
                        Array.average usableValues

                if Double.IsFinite result then
                    Some result
                else
                    None

    /// ionMobilityOf gives the ion mobility of the MS2 spectrum of a PSM by its PSMId.
    let estimateGaborTargetFromPsms
        (centrePSMs: int)
        (ionMobilityOf: string -> float)
        (correctionAt: float -> float -> float)
        theoMz
        (psms: (PSMStatisticsResult * float)[])
        =

        let usablePsms =
            psms
            |> Array.filter (fun (psm, weight) ->
                Double.IsFinite psm.ScanTime
                && Double.IsFinite (ionMobilityOf psm.PSMId)
                && Double.IsFinite weight
                && weight >= 0.0
            )

        if Array.isEmpty usablePsms || not (Double.IsFinite theoMz) then
            None
        else
            let bestPsms =
                usablePsms
                |> Array.sortByDescending snd
                |> Array.truncate centrePSMs

            let weights =
                bestPsms |> Array.map snd

            let scanTimes =
                bestPsms |> Array.map (fun (psm, _) -> psm.ScanTime)

            let mobilities =
                bestPsms |> Array.map (fun (psm, _) -> ionMobilityOf psm.PSMId)

            match
                tryWeightedMeanOrMean weights scanTimes,
                tryWeightedMeanOrMean weights mobilities
            with
            | Some expectedRt, Some expectedMobility ->
                let expectedMz =
                    correctionAt expectedRt theoMz + theoMz

                if Double.IsFinite expectedMz then
                    Some
                        {
                            ExpectedRt = expectedRt
                            ExpectedMz = expectedMz
                            ExpectedMobility = expectedMobility
                        }
                else
                    None

            | _ ->
                None

    let private toRtBin rtOrigin rt =
        int (floor ((rt - rtOrigin) / rtBinWidth))

    let private toMobilityBin mobilityOrigin mobility =
        int (floor ((mobility - mobilityOrigin) / mobilityBinWidth))

    let private rtBinCenter rtOrigin rtBin =
        rtOrigin + ((float rtBin + 0.5) * rtBinWidth)

    let private mobilityBinCenter mobilityOrigin mobilityBin =
        mobilityOrigin + ((float mobilityBin + 0.5) * mobilityBinWidth)

    let private clamp minValue maxValue value =
        if value < minValue then minValue
        elif value > maxValue then maxValue
        else value

    let private clampBin binCount bin =
        clamp 0 (binCount - 1) bin

    let private copyJaggedArray (data: float[][]) =
        data |> Array.map Array.copy

    let private validateRectangularGrid argumentName (intensities: float[][]) =
        if isNull intensities || intensities.Length = 0 then
            invalidArg argumentName "The intensity grid must not be empty."

        if isNull intensities.[0] || intensities.[0].Length = 0 then
            invalidArg argumentName "The intensity grid must contain at least one mobility bin."

        let mobilityCount =
            intensities.[0].Length

        let isRectangular =
            intensities
            |> Array.forall (fun row ->
                not (isNull row) && row.Length = mobilityCount
            )

        if not isRectangular then
            invalidArg argumentName "All rows of the intensity grid must have equal length."

    let extractGaborObservations (peakGroups: MzIO.Binary.Peak2D[][]) =
        peakGroups
        |> Array.collect id
        |> Array.choose (fun peak ->
            match peak.IonMobility with
            | Some mobility
                when Double.IsFinite peak.Rt
                  && Double.IsFinite mobility
                  && isUsableIntensity peak.Intensity ->
                Some
                    {
                        Rt = peak.Rt
                        Mobility = mobility
                        Intensity = peak.Intensity
                    }
            | _ ->
                None
        )

    /// Bins irregular MS1 observations into a regular RT-mobility matrix required by the Gabor convolution.
    let tryCreateGaborGrid (peakGroups: MzIO.Binary.Peak2D[][]) =
        let observations =
            extractGaborObservations peakGroups

        if Array.isEmpty observations then
            None
        else
            let rtOrigin =
                observations
                |> Array.minBy (fun observation -> observation.Rt)
                |> fun observation -> observation.Rt

            let mobilityOrigin =
                observations
                |> Array.minBy (fun observation -> observation.Mobility)
                |> fun observation -> observation.Mobility

            let binnedIntensities =
                System.Collections.Generic.Dictionary<int * int, float>()

            observations
            |> Array.iter (fun observation ->
                let rtBin =
                    toRtBin rtOrigin observation.Rt

                let mobilityBin =
                    toMobilityBin mobilityOrigin observation.Mobility

                let key =
                    rtBin, mobilityBin

                let oldIntensity =
                    match binnedIntensities.TryGetValue key with
                    | true, value -> value
                    | false, _ -> 0.0

                binnedIntensities.[key] <- oldIntensity + observation.Intensity
            )

            // The dictionary key is (rtBin, mobilityBin); therefore fst belongs to RT and snd to mobility.
            let rtBinCount =
                binnedIntensities.Keys
                |> Seq.maxBy fst
                |> fst
                |> fun maxBin -> maxBin + 1

            let mobilityBinCount =
                binnedIntensities.Keys
                |> Seq.maxBy snd
                |> snd
                |> fun maxBin -> maxBin + 1

            let intensities =
                Array.init rtBinCount (fun _ ->
                    Array.zeroCreate<float> mobilityBinCount
                )

            binnedIntensities
            |> Seq.iter (fun entry ->
                let rtBin, mobilityBin =
                    entry.Key

                intensities.[rtBin].[mobilityBin] <- entry.Value
            )

            let retentionTimes =
                Array.init rtBinCount (rtBinCenter rtOrigin)

            let mobilities =
                Array.init mobilityBinCount (mobilityBinCenter mobilityOrigin)

            Some
                {
                    RtOrigin = rtOrigin
                    MobilityOrigin = mobilityOrigin
                    RetentionTimes = retentionTimes
                    Mobilities = mobilities
                    RtBinCount = rtBinCount
                    MobilityBinCount = mobilityBinCount
                    Intensities = intensities
                }

    /// Applies one-dimensional AsLS along RT independently for every mobility bin.
    let correctBaselineAlongRt
        (baseLineCorrection: Domain.BaseLineCorrection option)
        (intensities: float[][])
        =

        validateRectangularGrid "intensities" intensities

        let rtCount =
            intensities.Length

        let mobilityCount =
            intensities.[0].Length

        let baselineGrid =
            Array.init rtCount (fun _ ->
                Array.zeroCreate<float> mobilityCount
            )

        let correctedGrid =
            Array.init rtCount (fun _ ->
                Array.zeroCreate<float> mobilityCount
            )

        match baseLineCorrection with
        | None ->
            {
                Baseline = baselineGrid
                CorrectedIntensities = copyJaggedArray intensities
            }

        | Some parameters ->
            for mobilityBin in 0 .. mobilityCount - 1 do
                let rawRtTrace =
                    Array.init rtCount (fun rtBin ->
                        let value = intensities.[rtBin].[mobilityBin]
                        if Double.IsFinite value && value > 0.0 then value else 0.0
                    )

                let containsSignal =
                    rawRtTrace |> Array.exists (fun value -> value > 0.0)

                // AsLS needs several samples to estimate curvature. Short or empty traces remain unchanged.
                if containsSignal && rawRtTrace.Length >= 3 then
                    let baselineTrace =
                        try
                            FSharp.Stats.Signal.Baseline.baselineAls'
                                parameters.MaxIterations
                                parameters.Lambda
                                parameters.P
                                rawRtTrace
                            |> Array.ofSeq
                        with
                        | _ ->
                            Array.zeroCreate<float> rawRtTrace.Length

                    let usableBaseline =
                        if baselineTrace.Length = rawRtTrace.Length then
                            baselineTrace
                        else
                            Array.zeroCreate<float> rawRtTrace.Length

                    for rtBin in 0 .. rtCount - 1 do
                        let baselineValue =
                            let value = usableBaseline.[rtBin]
                            if Double.IsFinite value then value else 0.0

                        baselineGrid.[rtBin].[mobilityBin] <- baselineValue
                        correctedGrid.[rtBin].[mobilityBin] <-
                            max 0.0 (rawRtTrace.[rtBin] - baselineValue)
                else
                    for rtBin in 0 .. rtCount - 1 do
                        correctedGrid.[rtBin].[mobilityBin] <- rawRtTrace.[rtBin]

            {
                Baseline = baselineGrid
                CorrectedIntensities = correctedGrid
            }

    let private validateGaborParameters (parameters: Domain.Gabor3DParams) =
        if parameters.sizeX <= 0 then
            invalidArg "sizeX" "Gabor sizeX must be greater than zero."

        if parameters.sizeY <= 0 then
            invalidArg "sizeY" "Gabor sizeY must be greater than zero."

        if not (Double.IsFinite parameters.sigmaX) || parameters.sigmaX <= 0.0 then
            invalidArg "sigmaX" "Gabor sigmaX must be positive and finite."

        if not (Double.IsFinite parameters.sigmaY) || parameters.sigmaY <= 0.0 then
            invalidArg "sigmaY" "Gabor sigmaY must be positive and finite."

        if not (Double.IsFinite parameters.frequency) then
            invalidArg "frequency" "Gabor frequency must be finite."

        if not (Double.IsFinite parameters.theta) then
            invalidArg "theta" "Gabor theta must be finite."

    let private ensureOdd value =
        let positiveValue = max 3 value
        if positiveValue % 2 = 0 then positiveValue + 1 else positiveValue

    let private kernelSizeFromSigma configuredSize sigma =
        let sigmaBasedSize =
            2 * int (ceil (3.0 * sigma)) + 1
            |> ensureOdd

        // Keep the filter bank computationally bounded by the explicitly configured kernel support.
        min (ensureOdd configuredSize) sigmaBasedSize
        |> ensureOdd

    /// Creates a sigma filter bank while retaining frequency and orientation from the configured kernel.
    let createGaborParameterBank (baseParameters: Domain.Gabor3DParams) =
        validateGaborParameters baseParameters

        [|
            for sigmaXScale in gaborSigmaScaleFactors do
                for sigmaYScale in gaborSigmaScaleFactors do
                    let sigmaX =
                        baseParameters.sigmaX * sigmaXScale

                    let sigmaY =
                        baseParameters.sigmaY * sigmaYScale

                    yield
                        {
                            baseParameters with
                                sizeX = kernelSizeFromSigma baseParameters.sizeX sigmaX
                                sizeY = kernelSizeFromSigma baseParameters.sizeY sigmaY
                                sigmaX = sigmaX
                                sigmaY = sigmaY
                        }
        |]
        |> Array.distinctBy (fun parameters ->
            parameters.sizeX,
            parameters.sizeY,
            parameters.sigmaX,
            parameters.sigmaY,
            parameters.frequency,
            parameters.theta
        )

    /// Reapplies complex L2 normalization so filter-bank scores are comparable even with older package builds.
    let private normalizeGaborKernelInPlace (kernel: Create.GaborWavelet3D) =
        let real =
            kernel.RealAnteil

        let imaginary =
            kernel.Imaginary

        validateRectangularGrid "kernel.RealAnteil" real
        validateRectangularGrid "kernel.Imaginary" imaginary

        if real.Length <> imaginary.Length || real.[0].Length <> imaginary.[0].Length then
            invalidArg "kernel" "Real and imaginary Gabor kernel parts must have identical dimensions."

        let mutable squaredNorm =
            0.0

        for x in 0 .. real.Length - 1 do
            for y in 0 .. real.[x].Length - 1 do
                let realValue = real.[x].[y]
                let imaginaryValue = imaginary.[x].[y]
                squaredNorm <-
                    squaredNorm
                    + realValue * realValue
                    + imaginaryValue * imaginaryValue

        let norm =
            sqrt squaredNorm

        if not (Double.IsFinite norm) || norm <= 0.0 then
            invalidArg "kernel" "Gabor kernel L2 energy must be positive and finite."

        for x in 0 .. real.Length - 1 do
            for y in 0 .. real.[x].Length - 1 do
                real.[x].[y] <- real.[x].[y] / norm
                imaginary.[x].[y] <- imaginary.[x].[y] / norm

        kernel

    let createGaborKernel (parameters: Domain.Gabor3DParams) =
        validateGaborParameters parameters

        Create.createGaborWavelet2D
            parameters.sizeX
            parameters.sizeY
            parameters.sigmaX
            parameters.sigmaY
            parameters.frequency
            parameters.theta
        |> normalizeGaborKernelInPlace

    let applyGaborToGrid (parameters: Domain.Gabor3DParams) (grid: GaborGrid) =
        let kernel =
            createGaborKernel parameters

        Create.applyGaborWavelet2D grid.Intensities kernel

    let private isInsideGrid (grid: GaborGrid) position =
        position.RtBin >= 0
        && position.RtBin < grid.RtBinCount
        && position.MobilityBin >= 0
        && position.MobilityBin < grid.MobilityBinCount

    let private isKernelFullyInsideGrid
        (grid: GaborGrid)
        (parameters: Domain.Gabor3DParams)
        position
        =

        let rtRadius =
            parameters.sizeX / 2

        let mobilityRadius =
            parameters.sizeY / 2

        position.RtBin - rtRadius >= 0
        && position.RtBin + rtRadius < grid.RtBinCount
        && position.MobilityBin - mobilityRadius >= 0
        && position.MobilityBin + mobilityRadius < grid.MobilityBinCount

    let private gaborMagnitudeAt (response: Create.gaborResponse) position =
        response.Magnitude.[position.RtBin].[position.MobilityBin]

    let private correctedIntensityAt (grid: GaborGrid) position =
        grid.Intensities.[position.RtBin].[position.MobilityBin]

    let private tryExpectedGridPosition
        (grid: GaborGrid)
        expectedRt
        expectedMobility
        =

        if
            not (Double.IsFinite expectedRt)
            || not (Double.IsFinite expectedMobility)
            || grid.RtBinCount <= 0
            || grid.MobilityBinCount <= 0
        then
            None
        else
            let expectedPosition =
                {
                    RtBin = toRtBin grid.RtOrigin expectedRt
                    MobilityBin = toMobilityBin grid.MobilityOrigin expectedMobility
                }

            if isInsideGrid grid expectedPosition then
                Some expectedPosition
            else
                None

    let private isValidGaborResponse
        (grid: GaborGrid)
        (response: Create.gaborResponse)
        =

        let magnitude =
            response.Magnitude

        not (isNull magnitude)
        && magnitude.Length = grid.RtBinCount
        && magnitude.Length > 0
        && (
            magnitude
            |> Array.forall (fun row ->
                not (isNull row)
                && row.Length = grid.MobilityBinCount
                && row.Length > 0
                && (row |> Array.forall Double.IsFinite)
            )
        )

    /// Finds the strongest response for one kernel inside the expected PSM search window.
    let private tryFindHighestGaborPositionNearExpectedPosition
        (grid: GaborGrid)
        (parameters: Domain.Gabor3DParams)
        (response: Create.gaborResponse)
        expectedRt
        expectedMobility
        =

        tryExpectedGridPosition grid expectedRt expectedMobility
        |> Option.bind (fun expectedPosition ->
            let minRtBin =
                expectedPosition.RtBin - gaborRtSearchRadiusBins
                |> clampBin grid.RtBinCount

            let maxRtBin =
                expectedPosition.RtBin + gaborRtSearchRadiusBins
                |> clampBin grid.RtBinCount

            let minMobilityBin =
                expectedPosition.MobilityBin - gaborMobilitySearchRadiusBins
                |> clampBin grid.MobilityBinCount

            let maxMobilityBin =
                expectedPosition.MobilityBin + gaborMobilitySearchRadiusBins
                |> clampBin grid.MobilityBinCount

            let createCandidates requireFullKernelSupport =
                [|
                    for rtBin in minRtBin .. maxRtBin do
                        for mobilityBin in minMobilityBin .. maxMobilityBin do
                            let position =
                                {
                                    RtBin = rtBin
                                    MobilityBin = mobilityBin
                                }

                            let supportIsUsable =
                                not requireFullKernelSupport
                                || isKernelFullyInsideGrid grid parameters position

                            if supportIsUsable then
                                let score =
                                    gaborMagnitudeAt response position

                                if Double.IsFinite score && score > 0.0 then
                                    yield position, score
                |]

            let fullSupportCandidates =
                createCandidates true

            // Fall back to the package's zero-padded edge behavior only when no fully supported position exists.
            let candidates =
                if Array.isEmpty fullSupportCandidates then
                    createCandidates false
                else
                    fullSupportCandidates

            if Array.isEmpty candidates then
                None
            else
                candidates
                |> Array.maxBy snd
                |> fst
                |> Some
        )

    let private tryEvaluateGaborKernel
        (grid: GaborGrid)
        expectedRt
        expectedMobility
        (parameters: Domain.Gabor3DParams)
        =

        if parameters.sizeX > grid.RtBinCount || parameters.sizeY > grid.MobilityBinCount then
            None
        else
            let response =
                applyGaborToGrid parameters grid

            if not (isValidGaborResponse grid response) then
                None
            else
                tryFindHighestGaborPositionNearExpectedPosition
                    grid
                    parameters
                    response
                    expectedRt
                    expectedMobility
                |> Option.bind (fun position ->
                    let magnitude =
                        gaborMagnitudeAt response position

                    if not (Double.IsFinite magnitude) || magnitude <= 0.0 then
                        None
                    else
                        Some
                            {
                                Parameters = parameters
                                Position = position
                                RetentionTime = grid.RetentionTimes.[position.RtBin]
                                IonMobility = grid.Mobilities.[position.MobilityBin]
                                CorrectedIntensityAtPosition = correctedIntensityAt grid position
                                Magnitude = magnitude
                            }
                )

    /// Executes the complete filter bank and retains both the best and second-best kernel match.
    let tryFindBestGaborMatch
        (parameterBank: Domain.Gabor3DParams[])
        (grid: GaborGrid)
        expectedRt
        expectedMobility
        =

        let sortedMatches =
            parameterBank
            |> Array.choose (
                tryEvaluateGaborKernel
                    grid
                    expectedRt
                    expectedMobility
            )
            |> Array.sortByDescending (fun matchResult -> matchResult.Magnitude)

        if Array.isEmpty sortedMatches then
            None
        else
            let bestMatch =
                sortedMatches.[0]

            let secondBestMatch =
                if sortedMatches.Length > 1 then
                    Some sortedMatches.[1]
                else
                    None

            let scoreRatio =
                secondBestMatch
                |> Option.bind (fun second ->
                    if Double.IsFinite second.Magnitude && second.Magnitude > 0.0 then
                        Some (bestMatch.Magnitude / second.Magnitude)
                    else
                        None
                )

            Some
                {
                    BestMatch = bestMatch
                    SecondBestMatch = secondBestMatch
                    ScoreRatio = scoreRatio
                }

    /// Crops the baseline-corrected grid to the PSM search area plus enough kernel padding.
    let tryCreateLocalGaborSearchGrid
        (parameterBank: Domain.Gabor3DParams[])
        (grid: GaborGrid)
        expectedRt
        expectedMobility
        =

        tryExpectedGridPosition grid expectedRt expectedMobility
        |> Option.map (fun expectedPosition ->
            let maximumRtKernelRadius =
                parameterBank
                |> Array.map (fun parameters -> parameters.sizeX / 2)
                |> fun values -> if Array.isEmpty values then 0 else Array.max values

            let maximumMobilityKernelRadius =
                parameterBank
                |> Array.map (fun parameters -> parameters.sizeY / 2)
                |> fun values -> if Array.isEmpty values then 0 else Array.max values

            let minRtBin =
                expectedPosition.RtBin
                - gaborRtSearchRadiusBins
                - maximumRtKernelRadius
                |> clampBin grid.RtBinCount

            let maxRtBin =
                expectedPosition.RtBin
                + gaborRtSearchRadiusBins
                + maximumRtKernelRadius
                |> clampBin grid.RtBinCount

            let minMobilityBin =
                expectedPosition.MobilityBin
                - gaborMobilitySearchRadiusBins
                - maximumMobilityKernelRadius
                |> clampBin grid.MobilityBinCount

            let maxMobilityBin =
                expectedPosition.MobilityBin
                + gaborMobilitySearchRadiusBins
                + maximumMobilityKernelRadius
                |> clampBin grid.MobilityBinCount

            let localRtCount =
                maxRtBin - minRtBin + 1

            let localMobilityCount =
                maxMobilityBin - minMobilityBin + 1

            let localIntensities =
                Array.init localRtCount (fun localRtBin ->
                    Array.init localMobilityCount (fun localMobilityBin ->
                        grid.Intensities.[minRtBin + localRtBin].[minMobilityBin + localMobilityBin]
                    )
                )

            let localRtOrigin =
                grid.RtOrigin + float minRtBin * rtBinWidth

            let localMobilityOrigin =
                grid.MobilityOrigin + float minMobilityBin * mobilityBinWidth

            {
                RtOrigin = localRtOrigin
                MobilityOrigin = localMobilityOrigin
                RetentionTimes =
                    grid.RetentionTimes.[minRtBin .. maxRtBin]
                Mobilities =
                    grid.Mobilities.[minMobilityBin .. maxMobilityBin]
                RtBinCount = localRtCount
                MobilityBinCount = localMobilityCount
                Intensities = localIntensities
            }
        )

    /// Converts the rotated Gabor envelope widths into axis-aligned standard deviations in grid bins.
    let private effectiveGaussianSigmasInBins (parameters: Domain.Gabor3DParams) =
        let cosTheta =
            cos parameters.theta

        let sinTheta =
            sin parameters.theta

        let rtEnvelopeVariance =
            parameters.sigmaX * parameters.sigmaX * cosTheta * cosTheta
            + parameters.sigmaY * parameters.sigmaY * sinTheta * sinTheta

        let mobilityEnvelopeVariance =
            parameters.sigmaX * parameters.sigmaX * sinTheta * sinTheta
            + parameters.sigmaY * parameters.sigmaY * cosTheta * cosTheta

        // The Gabor implementation uses exp(-x²/s²), whereas a standard Gaussian uses exp(-x²/(2σ²)).
        sqrt (rtEnvelopeVariance / 2.0),
        sqrt (mobilityEnvelopeVariance / 2.0)

    let createGaussianFitWindow
        (grid: GaborGrid)
        (bestMatch: GaborKernelMatch)
        =

        let sigmaRtBins, sigmaMobilityBins =
            effectiveGaussianSigmasInBins bestMatch.Parameters

        let rtRadius =
            max 2 (int (ceil (gaussianFitRadiusInSigma * sigmaRtBins)))

        let mobilityRadius =
            max 2 (int (ceil (gaussianFitRadiusInSigma * sigmaMobilityBins)))

        {
            MinRtBin = clampBin grid.RtBinCount (bestMatch.Position.RtBin - rtRadius)
            MaxRtBin = clampBin grid.RtBinCount (bestMatch.Position.RtBin + rtRadius)
            MinMobilityBin = clampBin grid.MobilityBinCount (bestMatch.Position.MobilityBin - mobilityRadius)
            MaxMobilityBin = clampBin grid.MobilityBinCount (bestMatch.Position.MobilityBin + mobilityRadius)
        }

    let extractGaussianObservations
        (grid: GaborGrid)
        (window: GaborFitWindow)
        =

        [|
            for rtBin in window.MinRtBin .. window.MaxRtBin do
                for mobilityBin in window.MinMobilityBin .. window.MaxMobilityBin do
                    let intensity =
                        grid.Intensities.[rtBin].[mobilityBin]

                    if Double.IsFinite intensity && intensity >= 0.0 then
                        yield
                            {
                                RetentionTime = grid.RetentionTimes.[rtBin]
                                IonMobility = grid.Mobilities.[mobilityBin]
                                Intensity = intensity
                            }
        |]

    let private gaussian2DValue
        (parameters: Gaussian2DParameters)
        retentionTime
        ionMobility
        =

        if
            parameters.Height < 0.0
            || parameters.SigmaRt <= 0.0
            || parameters.SigmaMobility <= 0.0
        then
            nan
        else
            let rtDistance =
                (retentionTime - parameters.MeanRt) / parameters.SigmaRt

            let mobilityDistance =
                (ionMobility - parameters.MeanMobility) / parameters.SigmaMobility

            parameters.Height
            * exp (-0.5 * (rtDistance * rtDistance + mobilityDistance * mobilityDistance))

    let private gaussianSse
        (parameters: Gaussian2DParameters)
        (observations: Gaussian2DObservation[])
        =

        if Array.isEmpty observations then
            infinity
        else
            observations
            |> Array.sumBy (fun observation ->
                let predicted =
                    gaussian2DValue
                        parameters
                        observation.RetentionTime
                        observation.IonMobility

                if not (Double.IsFinite predicted) then
                    infinity
                else
                    let residual =
                        observation.Intensity - predicted

                    residual * residual
            )

    let private gaussianParameterAt index parameters =
        match index with
        | 0 -> parameters.Height
        | 1 -> parameters.MeanRt
        | 2 -> parameters.MeanMobility
        | 3 -> parameters.SigmaRt
        | 4 -> parameters.SigmaMobility
        | _ -> invalidArg "index" "Gaussian parameter index must be between zero and four."

    let private setGaussianParameter index value parameters =
        match index with
        | 0 -> { parameters with Height = value }
        | 1 -> { parameters with MeanRt = value }
        | 2 -> { parameters with MeanMobility = value }
        | 3 -> { parameters with SigmaRt = value }
        | 4 -> { parameters with SigmaMobility = value }
        | _ -> invalidArg "index" "Gaussian parameter index must be between zero and four."

    let private clampGaussianParameters
        (lowerBounds: float[])
        (upperBounds: float[])
        (parameters: Gaussian2DParameters)
        : Gaussian2DParameters =

        [| 0 .. 4 |]
        |> Array.fold
            (fun (current: Gaussian2DParameters) (index: int) ->
                let boundedValue =
                    gaussianParameterAt index current
                    |> clamp lowerBounds.[index] upperBounds.[index]

                setGaussianParameter index boundedValue current
            )
            parameters

    let createInitialGaussianParameters
        (bestMatch: GaborKernelMatch)
        (observations: Gaussian2DObservation[])
        =

        let sigmaRtBins, sigmaMobilityBins =
            effectiveGaussianSigmasInBins bestMatch.Parameters

        let maximumIntensity =
            observations
            |> Array.map (fun observation -> observation.Intensity)
            |> Array.filter Double.IsFinite
            |> fun values ->
                if Array.isEmpty values then 0.0 else Array.max values

        {
            Height = maximumIntensity
            MeanRt = bestMatch.RetentionTime
            MeanMobility = bestMatch.IonMobility
            SigmaRt = max (rtBinWidth * 0.5) (sigmaRtBins * rtBinWidth)
            SigmaMobility = max (mobilityBinWidth * 0.5) (sigmaMobilityBins * mobilityBinWidth)
        }

    /// Fits the five Gaussian parameters with bounded coordinate pattern search.
    /// This keeps the implementation deterministic and avoids adding another optimizer dependency.
    let tryFitGaussian2D
        (window: GaborFitWindow)
        (grid: GaborGrid)
        (initialParameters: Gaussian2DParameters)
        (observations: Gaussian2DObservation[])
        =

        if observations.Length < 9 then
            None
        else
            let maximumIntensity =
                observations
                |> Array.maxBy (fun observation -> observation.Intensity)
                |> fun observation -> observation.Intensity

            if not (Double.IsFinite maximumIntensity) || maximumIntensity <= 0.0 then
                None
            else
                let minRt =
                    grid.RetentionTimes.[window.MinRtBin]

                let maxRt =
                    grid.RetentionTimes.[window.MaxRtBin]

                let minMobility =
                    grid.Mobilities.[window.MinMobilityBin]

                let maxMobility =
                    grid.Mobilities.[window.MaxMobilityBin]

                let rtSpan =
                    max rtBinWidth (maxRt - minRt + rtBinWidth)

                let mobilitySpan =
                    max mobilityBinWidth (maxMobility - minMobility + mobilityBinWidth)

                let lowerBounds =
                    [|
                        0.0
                        minRt
                        minMobility
                        rtBinWidth * 0.5
                        mobilityBinWidth * 0.5
                    |]

                let upperBounds =
                    [|
                        max (maximumIntensity * 5.0) (initialParameters.Height * 5.0)
                        maxRt
                        maxMobility
                        rtSpan
                        mobilitySpan
                    |]

                let mutable current =
                    clampGaussianParameters lowerBounds upperBounds initialParameters

                let mutable currentSse =
                    gaussianSse current observations

                if not (Double.IsFinite currentSse) then
                    None
                else
                    let mutable steps =
                        [|
                            max (maximumIntensity * 0.10) 1e-9
                            max (rtSpan * 0.10) rtBinWidth
                            max (mobilitySpan * 0.10) mobilityBinWidth
                            max (current.SigmaRt * 0.25) (rtBinWidth * 0.5)
                            max (current.SigmaMobility * 0.25) (mobilityBinWidth * 0.5)
                        |]

                    let minimumSteps =
                        [|
                            max (maximumIntensity * 1e-6) 1e-9
                            rtBinWidth * 0.01
                            mobilityBinWidth * 0.01
                            rtBinWidth * 0.01
                            mobilityBinWidth * 0.01
                        |]

                    let mutable iteration =
                        0

                    let mutable converged =
                        false

                    while iteration < gaussianFitMaximumIterations && not converged do
                        let mutable improvedInIteration =
                            false

                        for parameterIndex in 0 .. 4 do
                            let currentValue =
                                gaussianParameterAt parameterIndex current

                            let candidates =
                                [|
                                    current
                                    setGaussianParameter parameterIndex (currentValue + steps.[parameterIndex]) current
                                    setGaussianParameter parameterIndex (currentValue - steps.[parameterIndex]) current
                                |]
                                |> Array.map (clampGaussianParameters lowerBounds upperBounds)
                                |> Array.distinct
                                |> Array.map (fun candidate ->
                                    candidate, gaussianSse candidate observations
                                )
                                |> Array.filter (fun (_, sse) -> Double.IsFinite sse)

                            if not (Array.isEmpty candidates) then
                                let bestCandidate, bestSse =
                                    candidates |> Array.minBy snd

                                if bestSse < currentSse then
                                    current <- bestCandidate
                                    currentSse <- bestSse
                                    improvedInIteration <- true

                        if not improvedInIteration then
                            steps <- steps |> Array.map (fun step -> step * 0.5)

                        converged <-
                            Array.forall2 (fun step minimumStep -> step <= minimumStep) steps minimumSteps

                        iteration <- iteration + 1

                    let fittedValues =
                        [|
                            current.Height
                            current.MeanRt
                            current.MeanMobility
                            current.SigmaRt
                            current.SigmaMobility
                        |]

                    if
                        (fittedValues |> Array.forall Double.IsFinite)
                        && current.Height > 0.0
                        && current.SigmaRt > 0.0
                        && current.SigmaMobility > 0.0
                    then
                        Some current
                    else
                        None

    let private calculateGaussianVolume parameters =
        2.0
        * Math.PI
        * parameters.Height
        * parameters.SigmaRt
        * parameters.SigmaMobility

    let private calculateDirectWindowVolume (observations: Gaussian2DObservation[]) =
        observations
        |> Array.sumBy (fun observation -> observation.Intensity)
        |> fun intensitySum -> intensitySum * rtBinWidth * mobilityBinWidth

    let private calculateGaussianRmse parameters (observations: Gaussian2DObservation[]) =
        if Array.isEmpty observations then
            nan
        else
            gaussianSse parameters observations
            |> fun sse -> sqrt (sse / float observations.Length)

    let private projectFitWindowToRtTrace
        (grid: GaborGrid)
        (window: GaborFitWindow)
        =

        [| window.MinRtBin .. window.MaxRtBin |]
        |> Array.map (fun rtBin ->
            let intensity =
                [| window.MinMobilityBin .. window.MaxMobilityBin |]
                |> Array.sumBy (fun mobilityBin ->
                    grid.Intensities.[rtBin].[mobilityBin]
                )

            grid.RetentionTimes.[rtBin], intensity
        )
        |> Array.unzip

    /// Quantifies one PSM using RT-wise AsLS, a normalized Gabor filter bank, and a fitted 2D Gaussian.
    let tryQuantifyGabor2D
        (baseLineCorrection: Domain.BaseLineCorrection option)
        (baseParameters: Domain.Gabor3DParams)
        expectedRt
        expectedMobility
        (peakGroups: MzIO.Binary.Peak2D[][])
        =

        peakGroups
        |> tryCreateGaborGrid
        |> Option.bind (fun rawGrid ->
            let baselineResult =
                correctBaselineAlongRt
                    baseLineCorrection
                    rawGrid.Intensities

            let correctedGrid =
                {
                    rawGrid with
                        Intensities = baselineResult.CorrectedIntensities
                }

            let parameterBank =
                createGaborParameterBank baseParameters

            tryCreateLocalGaborSearchGrid
                parameterBank
                correctedGrid
                expectedRt
                expectedMobility
            |> Option.bind (fun searchGrid ->
                tryFindBestGaborMatch
                    parameterBank
                    searchGrid
                    expectedRt
                    expectedMobility
                |> Option.bind (fun filterBankResult ->
                    let fitWindow =
                        createGaussianFitWindow
                            searchGrid
                            filterBankResult.BestMatch

                    let observations =
                        extractGaussianObservations
                            searchGrid
                            fitWindow

                    let initialParameters =
                        createInitialGaussianParameters
                            filterBankResult.BestMatch
                            observations

                    tryFitGaussian2D
                        fitWindow
                        searchGrid
                        initialParameters
                        observations
                    |> Option.bind (fun fittedParameters ->
                        let volume =
                            calculateGaussianVolume fittedParameters

                        let directWindowVolume =
                            calculateDirectWindowVolume observations

                        let fitRmse =
                            calculateGaussianRmse fittedParameters observations

                        if
                            not (Double.IsFinite volume)
                            || volume <= 0.0
                            || not (Double.IsFinite directWindowVolume)
                            || directWindowVolume <= 0.0
                            || not (Double.IsFinite fitRmse)
                        then
                            None
                        else
                            let measuredApex =
                                observations
                                |> Array.maxBy (fun observation -> observation.Intensity)

                            let projectedRetentionTimes, projectedIntensities =
                                projectFitWindowToRtTrace searchGrid fitWindow

                            Some
                                {
                                    Volume = volume
                                    DirectWindowVolume = directWindowVolume
                                    ApexIntensity = measuredApex.Intensity
                                    ApexRetentionTime = measuredApex.RetentionTime
                                    ApexIonMobility = measuredApex.IonMobility
                                    WeightedRetentionTime = fittedParameters.MeanRt
                                    WeightedIonMobility = fittedParameters.MeanMobility
                                    ProjectedRetentionTimes = projectedRetentionTimes
                                    ProjectedIntensities = projectedIntensities
                                    BestGaborMatch = filterBankResult.BestMatch
                                    SecondBestGaborMatch = filterBankResult.SecondBestMatch
                                    GaborScoreRatio = filterBankResult.ScoreRatio
                                    InitialGaussianParameters = initialParameters
                                    FittedGaussianParameters = fittedParameters
                                    FitRmse = fitRmse
                                }
                    )
                )
            )
        )

    let initGetGabor2DQuantification
        getGaborPeaks
        idx
        scanTimeWindow
        mzWindow
        mobilityWindow
        baseLineCorrection
        =

        fun
            (parameters: Domain.Gabor3DParams)
            expectedRt
            expectedMz
            expectedMobility ->

            let rtQuery =
                Query.createRangeQuery expectedRt scanTimeWindow

            let mzQuery =
                Query.createRangeQuery expectedMz mzWindow

            let ionMobilityQuery =
                Query.createRangeQuery expectedMobility mobilityWindow

            getGaborPeaks idx rtQuery mzQuery ionMobilityQuery
            |> tryQuantifyGabor2D
                baseLineCorrection
                parameters
                expectedRt
                expectedMobility

        
    ///
    let quantifyPeptides diagCharts zipCharts (processParams:Domain.QuantificationParams) (outputDir:string) (cn:SQLiteConnection) (instrumentOutput:string) (scoredPSMs:string)  =
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension scoredPSMs)
        logger.Trace (sprintf "Input file: %s" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        logger.Trace (sprintf "Now performing Quantification using: %s and %s, Results will be written to: %s" instrumentOutput scoredPSMs outputDir)
        logger.Trace "Read scored PSMs"
        ///
        let qpsms = QpsmInput.read scoredPSMs
        logger.Trace "Read scored PSMs:finished"
        // initialize Reader and Transaction
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".quant"
            Path.Combine [|outputDir;fileName|]
        logger.Trace (sprintf "outFilePath:%s" outFilePath)
        if Array.isEmpty qpsms then
            logger.Warn (sprintf "%s contains no PSMs, %s gets no quantifications." scoredPSMs outFilePath)
            Seq.empty<QuantificationResult>
            |> SeqIO'.csv "\t" true false
            |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        else
        //
        let plotDirectory =
            let fileName = sprintf "%s_plots" (Path.GetFileNameWithoutExtension instrumentOutput)
            let path = Path.Combine [|outputDir;fileName|]
            if System.IO.Directory.Exists path then
                path
            else
                System.IO.Directory.CreateDirectory path |> ignore
                path
        logger.Trace (sprintf "plotDirectory:%s" plotDirectory)
        // the peptides are looked up in the data base file, a copy of the data base in memory would take its whole size (GBs) for two lookups per ion;
        // every file has its own connection and the lookups run without a transaction, so that files in parallel do not lock each other
        use peptideDB = new SQLiteConnection(cn.ConnectionString)
        peptideDB.Open()
        logger.Trace "Get peptide lookUp function"
        let dBParams     = getSDBParamsByCn peptideDB
        //let massLookUp = prepareSelectMassByModSequenceAndGlobalMod peptideDB
        let peptideLookUp = getThreadSafePeptideLookUpFromFileBySequenceAndGMod peptideDB dBParams
        use isotopePatterns = new IsotopePatterns.IsotopePatternLookUp(logger, peptideDB)
        let calcIonSeries aal = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal
        logger.Trace "Get peptide lookUp function: finished"
        // initialize Reader and Transaction
        logger.Trace "Init connection to mass spectrum data."

        if
            not (
                instrumentOutput.EndsWith(
                    ".mzlite",
                    StringComparison.OrdinalIgnoreCase
                )
            )
        then
            invalidArg
                "instrumentOutput"
                "PSMBasedQuantificationTIMs currently requires an .mzlite input file."

        let inReader =
            Core.MzIO.Reader.getReader instrumentOutput
            :?> MzIO.MzSQL.MzSQL
        inReader.Connection.Open()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader       
        let inTr = inReader.BeginTransaction()
        let countMatchedMasses (peptide: LookUpResult<AminoAcids.AminoAcid>)(psms: PSMStatisticsResult []) =
            let frag = 
                let ionSeries = (calcIonSeries peptide.BioSequence).TargetMasses
                [1. .. 2.]
                |> List.collect (fun ch -> 
                    ionSeries 
                    |> List.map (fun x -> Mass.toMZ x.MainPeak.Mass ch)
                )                  
            psms
            |> Array.map (fun psm -> 
                let spec = inReader.ReadSpectrumPeaks psm.PSMId
                let sum = 
                    spec.Peaks 
                    |> Seq.filter (fun peak -> 
                        frag
                        |> List.exists (fun ion -> abs (ion - peak.Mz) <= (Mass.deltaMassByPpm 100. peak.Mz))
                    )
                    |> Seq.sumBy (fun x -> x.Intensity)
                psm,sum
            )

        logger.Trace "Create RetentionTime index, read the ms1s and the ion mobility of the PSM spectra"
        /// retTimeIdxed and ms1SortedByScanTime: the MS1s sorted by their scanTime. ionMobilityBySpectrum: the ion
        /// mobility of every spectrum a PSM refers to, which the .qpsm does not carry.
        let retTimeIdxed, ms1SortedByScanTime, ionMobilityBySpectrum =
            let psmSpectrumIds = Collections.Generic.HashSet<string>(qpsms |> Array.map (fun x -> x.PSMId))
            QpsmInput.readSpectrumIndex inReader inRunID (fun message -> logger.Warn message) psmSpectrumIds
        let ionMobilityOf (psmId: string) =
            match ionMobilityBySpectrum.TryGetValue psmId with
            | true, mobility -> mobility
            | _ -> nan
        let psmsNotInMzlite = qpsms |> Array.filter (fun x -> not (ionMobilityBySpectrum.ContainsKey x.PSMId)) |> Array.length
        let psmsWithIonMobility = qpsms |> Array.filter (fun x -> Double.IsFinite (ionMobilityOf x.PSMId)) |> Array.length
        let psmsWithoutIonMobility = qpsms.Length - psmsNotInMzlite - psmsWithIonMobility
        logger.Trace (sprintf "Create RetentionTime index, read the ms1s and the ion mobility of the PSM spectra:finished, %i MS1s, %i of %i PSMs have an ion mobility in their spectrum" ms1SortedByScanTime.Length psmsWithIonMobility qpsms.Length)
        if qpsms.Length > 0 && psmsWithIonMobility = 0 then
            failwithf "No PSM of %s refers to a spectrum of %s with an ion mobility (MS:1002815): %i PSMIds not found or not readable in the mzlite (wrong mzlite?), %i spectra without ion mobility. The PSMId has to be the ID of an MS2 spectrum of this mzlite." scoredPSMs instrumentOutput psmsNotInMzlite psmsWithoutIonMobility
        if psmsWithIonMobility < qpsms.Length then
            logger.Warn (sprintf "%i PSMIds not found or not readable in the mzlite (wrong mzlite?), %i spectra without ion mobility of %i PSMs in %s, they do not count for the ion mobility of their ion." psmsNotInMzlite psmsWithoutIonMobility qpsms.Length instrumentOutput)
        
        let calibration = calibratePsms logger diagCharts plotDirectory qpsms

        logger.Trace "init lookup functions"        
        ///
        
        let mzWindow = MzCalibration.mzWindow processParams.XicExtraction.MzWindow_Da calibration
        logger.Trace (sprintf "XIC %s" (MzCalibration.describeWindow processParams.XicExtraction.MzWindow_Da mzWindow (calibration.HalfWindow = MzCalibration.fallbackMzWindow)))

        let readSpecPeaksWithMem = FSharpAux.Memoization.memoize inReader.ReadSpectrumPeaks
        ///
        let readSortedPeaks = Query.readSortedPeaks inReader
        // The ions are quantified in the order of their scan time and their XICs read the MS1 spectra of a window of
        // 2 x ScanTimeWindow, so the spectra of twice the widest window are kept instead of all MS1 spectra of the file.
        let peakCache =
            let cacheSize, widestWindow = Xic.peakCacheSize (ms1SortedByScanTime |> Array.map fst) processParams.XicExtraction.ScanTimeWindow
            logger.Trace (sprintf "MS1 peaks cache: %i spectra (%i MS1s in the widest XIC window)" cacheSize widestWindow)
            Xic.PeakCache(readSortedPeaks, cacheSize)
        let readSortedPeaksWithMem = peakCache.Get
        let closestMs1 = SpectrumIndex.initClosestMs1 ms1SortedByScanTime
        let compareGaborPredictedAndMeasuredIsotopicCluster =
            initComparePredictedAndMeasuredIsotopicCluster
                isotopePatterns
                (fun id lowMz highMz -> getSpec inReader id |> Array.filter (fun x -> x.Mz < highMz && x.Mz > lowMz))
                closestMs1
                (calibration.HalfWindow / MzCalibration.mzWindowSdFactor)
        let ionMobilityHalfWindow =
            match processParams.XicExtraction.IonMobilityWindow with
            | Domain.Window.Fixed v -> v
            | Domain.Window.Estimate -> failwith "IonMobilityWindow Estimate is only available in AlignmentBasedQuantification; use Fixed."
        let getXIC = initGetProcessedXIC processParams.BaseLineCorrection readSortedPeaksWithMem retTimeIdxed processParams.XicExtraction.ScanTimeWindow mzWindow ionMobilityHalfWindow
        // the XIC is summed over +-IonMobilityWindow around the ion mobility centre of the ion, taken from the signal around the ion mobility of the PSMs
        logger.Trace (sprintf "Ion mobility: XIC summed over +-%.4f around the ion mobility centre of the ion" ionMobilityHalfWindow)
        // for the log, the ions with a finite expected ion mobility are counted
        let centreTally = Xic.MobilityCentreTally()
        let getTargetXic =
            let tryGetTargetXic = initGetTargetXic processParams.BaseLineCorrection readSortedPeaksWithMem retTimeIdxed processParams.XicExtraction.ScanTimeWindow mzWindow ionMobilityHalfWindow
            fun scanTime mz expectedMobility ->
                let centre, xic = tryGetTargetXic scanTime mz expectedMobility
                if Numerics.isFinite expectedMobility then centreTally.Add(expectedMobility, centre)
                defaultArg centre expectedMobility, xic
        // the isotope pattern is read from the XIC points of its peaks inside the ion mobility window of the XIC
        let comparePredictedAndMeasuredIsotopicCluster =
            initComparePredictedAndMeasuredIsotopicClusterIonMobility
                isotopePatterns
                (fun id lockMz halfWindow imLow imHigh -> Xic.extractPoint (readSortedPeaksWithMem id) lockMz halfWindow imLow imHigh)
                closestMs1
                mzWindow
                ionMobilityHalfWindow
        ///
        let identifyPeaks = initIdentifyPeaks processParams.XicExtraction.XicProcessing
        logger.Trace "init lookup functions:finished"
        
        let getGaborPeaks =
            Query.initRTProfileGabor readSpecPeaksWithMem

        let getGabor2DQuantification =
            initGetGabor2DQuantification
                getGaborPeaks
                retTimeIdxed
                processParams.XicExtraction.ScanTimeWindow
                mzWindow
                ionMobilityHalfWindow
                processParams.BaseLineCorrection
        logger.Trace "init quantification functions"

        let emptyClusterComparison =
            {
                PeakComparisons = [||]
                KLDiv_UnCorrected = nan
                KLDiv_Corrected = nan
            }

        let compareGaborIsotopicCluster
            (quant: Gabor2DQuantification)
            charge
            peptide
            targetMz
            =

            try
                if Array.isEmpty quant.ProjectedRetentionTimes || Array.isEmpty quant.ProjectedIntensities then
                    emptyClusterComparison
                else
                    compareGaborPredictedAndMeasuredIsotopicCluster
                        quant.ProjectedRetentionTimes
                        quant.ProjectedIntensities
                        quant.ProjectedIntensities
                        charge
                        peptide
                        quant.WeightedRetentionTime
                        targetMz
            with
            | ex ->
                logger.Trace (sprintf "Gabor isotope cluster comparison failed: %A" ex)
                emptyClusterComparison

        let quantMzAtScanTime (peptide: LookUpResult<AminoAcids.AminoAcid>) charge scanTime =
            let theoMz = Mass.toMZ peptide.Mass (float charge)
            calibration.CorrectionAt scanTime theoMz + theoMz

        let createGaborQuantSide
            charge
            (peptide: LookUpResult<AminoAcids.AminoAcid>)
            expectedRt
            quantMz
            (quant: Gabor2DQuantification option)
            =

            match quant with
            | Some q ->
                let cluster =
                    compareGaborIsotopicCluster q charge peptide quantMz

                q.Volume,
                q.ApexIntensity,
                q.FitRmse,
                // Parameter order: kernel sizeX, sizeY, sigmaX, sigmaY, frequency, theta,
                // Gabor magnitude, best/second-best score ratio, Gaussian height, mean RT,
                // mean mobility, sigma RT, sigma mobility.
                [|
                    float q.BestGaborMatch.Parameters.sizeX
                    float q.BestGaborMatch.Parameters.sizeY
                    q.BestGaborMatch.Parameters.sigmaX
                    q.BestGaborMatch.Parameters.sigmaY
                    q.BestGaborMatch.Parameters.frequency
                    q.BestGaborMatch.Parameters.theta
                    q.BestGaborMatch.Magnitude
                    (q.GaborScoreRatio |> Option.defaultValue nan)
                    q.FittedGaussianParameters.Height
                    q.FittedGaussianParameters.MeanRt
                    q.FittedGaussianParameters.MeanMobility
                    q.FittedGaussianParameters.SigmaRt
                    q.FittedGaussianParameters.SigmaMobility
                |],
                expectedRt - q.WeightedRetentionTime,
                cluster,
                q.ProjectedRetentionTimes,
                q.ProjectedIntensities
            | None ->
                nan,
                nan,
                nan,
                ([||] : float[]),
                nan,
                emptyClusterComparison,
                ([||] : float[]),
                ([||] : float[])

        let gaborLabeledQuantification
            (gaborParams: Domain.Gabor3DParams)
            (pepIon: PeptideIon)
            (psms: PSMStatisticsResult[])
            =

            try
                let bestQValue, bestPepValue, prots =
                    psms
                    |> bestPsmByQValue
                    |> fun x -> x.QValue, x.PEPValue, x.ProteinNames

                let unlabledPeptide =
                    peptideLookUp pepIon.Sequence 0

                let labeledPeptide =
                    peptideLookUp pepIon.Sequence 1

                let targetPeptide =
                    if pepIon.GlobalMod = 0 then unlabledPeptide else labeledPeptide

                let psmsWithMatchedSums =
                    countMatchedMasses targetPeptide psms

                let theoMz =
                    Mass.toMZ targetPeptide.Mass (float pepIon.Charge)

                match estimateGaborTargetFromPsms processParams.CentrePSMs ionMobilityOf calibration.CorrectionAt theoMz psmsWithMatchedSums with
                | None ->
                    logger.Trace (sprintf "Gabor quant failed: No PSM target estimate, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                    None
                | Some target ->
                    match getGabor2DQuantification gaborParams target.ExpectedRt target.ExpectedMz target.ExpectedMobility with
                    | None ->
                        logger.Trace (sprintf "Gabor quant failed: No valid Gabor-Gaussian fit, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    | Some targetQuant ->
                        let meanScore =
                            psmsWithMatchedSums
                            |> Array.averageBy (fun (psm, _) -> psm.ModelScore)

                        if pepIon.GlobalMod = 0 then
                            let quantMzLight =
                                quantMzAtScanTime unlabledPeptide pepIon.Charge targetQuant.WeightedRetentionTime

                            let mzHeavy =
                                quantMzAtScanTime labeledPeptide pepIon.Charge targetQuant.WeightedRetentionTime

                            let inferredHeavy =
                                getGabor2DQuantification
                                    gaborParams
                                    targetQuant.WeightedRetentionTime
                                    mzHeavy
                                    targetQuant.WeightedIonMobility

                            let qLight, apexLight, seoLight, paramsLight, diffLight, clusterLight, rtTraceLight, intensityTraceLight =
                                createGaborQuantSide pepIon.Charge unlabledPeptide target.ExpectedRt quantMzLight (Some targetQuant)

                            let qHeavy, apexHeavy, seoHeavy, paramsHeavy, diffHeavy, clusterHeavy, rtTraceHeavy, intensityTraceHeavy =
                                createGaborQuantSide pepIon.Charge labeledPeptide targetQuant.WeightedRetentionTime mzHeavy inferredHeavy

                            let avgMass =
                                Mass.ofMZ quantMzLight (float pepIon.Charge)

                            {
                                StringSequence                              = pepIon.Sequence
                                GlobalMod                                   = pepIon.GlobalMod
                                Charge                                      = pepIon.Charge
                                PepSequenceID                               = pepIon.PepSequenceID
                                ModSequenceID                               = pepIon.ModSequenceID
                                PrecursorMZ                                 = quantMzLight
                                MeasuredMass                                = avgMass
                                TheoMass                                    = unlabledPeptide.Mass
                                AbsDeltaMass                                = abs(avgMass - unlabledPeptide.Mass)
                                MeanPercolatorScore                         = meanScore
                                QValue                                      = bestQValue
                                PEPValue                                    = bestPepValue
                                ProteinNames                                = prots
                                QuantMz_Light                               = quantMzLight
                                Quant_Light                                 = qLight
                                MeasuredApex_Light                          = apexLight
                                Seo_Light                                   = seoLight
                                Params_Light                                = paramsLight
                                Difference_SearchRT_FittedRT_Light          = diffLight
                                KLDiv_Observed_Theoretical_Light            = clusterLight.KLDiv_UnCorrected
                                KLDiv_CorrectedObserved_Theoretical_Light   = clusterLight.KLDiv_Corrected
                                QuantMz_Heavy                               = mzHeavy
                                Quant_Heavy                                 = qHeavy
                                MeasuredApex_Heavy                          = apexHeavy
                                Seo_Heavy                                   = seoHeavy
                                Params_Heavy                                = paramsHeavy
                                Difference_SearchRT_FittedRT_Heavy          = diffHeavy
                                KLDiv_Observed_Theoretical_Heavy            = clusterHeavy.KLDiv_UnCorrected
                                KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterHeavy.KLDiv_Corrected
                                Correlation_Light_Heavy                     = nan
                                QuantificationSource                        = QuantificationSource.PSM
                                IsotopicPatternMz_Light                     = clusterLight.PeakComparisons |> Array.map (fun x -> x.Mz)
                                IsotopicPatternIntensity_Observed_Light     = clusterLight.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                                IsotopicPatternIntensity_Corrected_Light    = clusterLight.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                                RtTrace_Light                               = rtTraceLight
                                IntensityTrace_Observed_Light               = intensityTraceLight
                                IntensityTrace_Corrected_Light              = intensityTraceLight
                                IsotopicPatternMz_Heavy                     = clusterHeavy.PeakComparisons |> Array.map (fun x -> x.Mz)
                                IsotopicPatternIntensity_Observed_Heavy     = clusterHeavy.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                                IsotopicPatternIntensity_Corrected_Heavy    = clusterHeavy.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                                RtTrace_Heavy                               = rtTraceHeavy
                                IntensityTrace_Observed_Heavy               = intensityTraceHeavy
                                IntensityTrace_Corrected_Heavy              = intensityTraceHeavy
                                AlignmentScore                              = nan
                                AlignmentQValue                             = 0.
                                IonMobility                                 = targetQuant.WeightedIonMobility
                                IonMobilityCentre                           = nan
                            }
                            |> Some
                        else
                            let quantMzHeavy =
                                quantMzAtScanTime labeledPeptide pepIon.Charge targetQuant.WeightedRetentionTime

                            let mzLight =
                                quantMzAtScanTime unlabledPeptide pepIon.Charge targetQuant.WeightedRetentionTime

                            let inferredLight =
                                getGabor2DQuantification
                                    gaborParams
                                    targetQuant.WeightedRetentionTime
                                    mzLight
                                    targetQuant.WeightedIonMobility

                            let qLight, apexLight, seoLight, paramsLight, diffLight, clusterLight, rtTraceLight, intensityTraceLight =
                                createGaborQuantSide pepIon.Charge unlabledPeptide targetQuant.WeightedRetentionTime mzLight inferredLight

                            let qHeavy, apexHeavy, seoHeavy, paramsHeavy, diffHeavy, clusterHeavy, rtTraceHeavy, intensityTraceHeavy =
                                createGaborQuantSide pepIon.Charge labeledPeptide target.ExpectedRt quantMzHeavy (Some targetQuant)

                            let avgMass =
                                Mass.ofMZ quantMzHeavy (float pepIon.Charge)

                            {
                                StringSequence                              = pepIon.Sequence
                                GlobalMod                                   = pepIon.GlobalMod
                                Charge                                      = pepIon.Charge
                                PepSequenceID                               = pepIon.PepSequenceID
                                ModSequenceID                               = pepIon.ModSequenceID
                                PrecursorMZ                                 = quantMzHeavy
                                MeasuredMass                                = avgMass
                                TheoMass                                    = labeledPeptide.Mass
                                AbsDeltaMass                                = abs(avgMass - labeledPeptide.Mass)
                                MeanPercolatorScore                         = meanScore
                                QValue                                      = bestQValue
                                PEPValue                                    = bestPepValue
                                ProteinNames                                = prots
                                QuantMz_Light                               = mzLight
                                Quant_Light                                 = qLight
                                MeasuredApex_Light                          = apexLight
                                Seo_Light                                   = seoLight
                                Params_Light                                = paramsLight
                                Difference_SearchRT_FittedRT_Light          = diffLight
                                KLDiv_Observed_Theoretical_Light            = clusterLight.KLDiv_UnCorrected
                                KLDiv_CorrectedObserved_Theoretical_Light   = clusterLight.KLDiv_Corrected
                                QuantMz_Heavy                               = quantMzHeavy
                                Quant_Heavy                                 = qHeavy
                                MeasuredApex_Heavy                          = apexHeavy
                                Seo_Heavy                                   = seoHeavy
                                Params_Heavy                                = paramsHeavy
                                Difference_SearchRT_FittedRT_Heavy          = diffHeavy
                                KLDiv_Observed_Theoretical_Heavy            = clusterHeavy.KLDiv_UnCorrected
                                KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterHeavy.KLDiv_Corrected
                                Correlation_Light_Heavy                     = nan
                                QuantificationSource                        = QuantificationSource.PSM
                                IsotopicPatternMz_Light                     = clusterLight.PeakComparisons |> Array.map (fun x -> x.Mz)
                                IsotopicPatternIntensity_Observed_Light     = clusterLight.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                                IsotopicPatternIntensity_Corrected_Light    = clusterLight.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                                RtTrace_Light                               = rtTraceLight
                                IntensityTrace_Observed_Light               = intensityTraceLight
                                IntensityTrace_Corrected_Light              = intensityTraceLight
                                IsotopicPatternMz_Heavy                     = clusterHeavy.PeakComparisons |> Array.map (fun x -> x.Mz)
                                IsotopicPatternIntensity_Observed_Heavy     = clusterHeavy.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                                IsotopicPatternIntensity_Corrected_Heavy    = clusterHeavy.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                                RtTrace_Heavy                               = rtTraceHeavy
                                IntensityTrace_Observed_Heavy               = intensityTraceHeavy
                                IntensityTrace_Corrected_Heavy              = intensityTraceHeavy
                                AlignmentScore                              = nan
                                AlignmentQValue                             = 0.
                                IonMobility                                 = targetQuant.WeightedIonMobility
                                IonMobilityCentre                           = nan
                            }
                            |> Some
            with
            | ex ->
                logger.Trace (sprintf "Gabor quantfailed: %A" ex)
                None

        let gaborLabelFreeQuantification
            (gaborParams: Domain.Gabor3DParams)
            (pepIon: PeptideIon)
            (psms: PSMStatisticsResult[])
            =

            try
                if pepIon.GlobalMod <> 0 then
                    None
                else
                    let bestQValue, bestPepValue, prots =
                        psms
                        |> bestPsmByQValue
                        |> fun x -> x.QValue, x.PEPValue, x.ProteinNames

                    let unlabledPeptide =
                        peptideLookUp pepIon.Sequence 0

                    let psmsWithMatchedSums =
                        countMatchedMasses unlabledPeptide psms

                    let theoMz =
                        Mass.toMZ unlabledPeptide.Mass (float pepIon.Charge)

                    match estimateGaborTargetFromPsms processParams.CentrePSMs ionMobilityOf calibration.CorrectionAt theoMz psmsWithMatchedSums with
                    | None ->
                        logger.Trace (sprintf "Gabor quant failed: No PSM target estimate, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    | Some target ->
                        match getGabor2DQuantification gaborParams target.ExpectedRt target.ExpectedMz target.ExpectedMobility with
                        | None ->
                            logger.Trace (sprintf "Gabor quant failed: No valid Gabor-Gaussian fit, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                            None
                        | Some targetQuant ->
                            let meanScore =
                                psmsWithMatchedSums
                                |> Array.averageBy (fun (psm, _) -> psm.ModelScore)

                            let quantMzLight =
                                quantMzAtScanTime unlabledPeptide pepIon.Charge targetQuant.WeightedRetentionTime

                            let qLight, apexLight, seoLight, paramsLight, diffLight, clusterLight, rtTraceLight, intensityTraceLight =
                                createGaborQuantSide pepIon.Charge unlabledPeptide target.ExpectedRt quantMzLight (Some targetQuant)

                            let avgMass =
                                Mass.ofMZ quantMzLight (float pepIon.Charge)

                            {
                                StringSequence                              = pepIon.Sequence
                                GlobalMod                                   = pepIon.GlobalMod
                                Charge                                      = pepIon.Charge
                                PepSequenceID                               = pepIon.PepSequenceID
                                ModSequenceID                               = pepIon.ModSequenceID
                                PrecursorMZ                                 = quantMzLight
                                MeasuredMass                                = avgMass
                                TheoMass                                    = unlabledPeptide.Mass
                                AbsDeltaMass                                = abs(avgMass - unlabledPeptide.Mass)
                                MeanPercolatorScore                         = meanScore
                                QValue                                      = bestQValue
                                PEPValue                                    = bestPepValue
                                ProteinNames                                = prots
                                QuantMz_Light                               = quantMzLight
                                Quant_Light                                 = qLight
                                MeasuredApex_Light                          = apexLight
                                Seo_Light                                   = seoLight
                                Params_Light                                = paramsLight
                                Difference_SearchRT_FittedRT_Light          = diffLight
                                KLDiv_Observed_Theoretical_Light            = clusterLight.KLDiv_UnCorrected
                                KLDiv_CorrectedObserved_Theoretical_Light   = clusterLight.KLDiv_Corrected
                                QuantMz_Heavy                               = nan
                                Quant_Heavy                                 = nan
                                MeasuredApex_Heavy                          = nan
                                Seo_Heavy                                   = nan
                                Params_Heavy                                = [||]
                                Difference_SearchRT_FittedRT_Heavy          = nan
                                KLDiv_Observed_Theoretical_Heavy            = nan
                                KLDiv_CorrectedObserved_Theoretical_Heavy   = nan
                                Correlation_Light_Heavy                     = nan
                                QuantificationSource                        = QuantificationSource.PSM
                                IsotopicPatternMz_Light                     = clusterLight.PeakComparisons |> Array.map (fun x -> x.Mz)
                                IsotopicPatternIntensity_Observed_Light     = clusterLight.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                                IsotopicPatternIntensity_Corrected_Light    = clusterLight.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                                RtTrace_Light                               = rtTraceLight
                                IntensityTrace_Observed_Light               = intensityTraceLight
                                IntensityTrace_Corrected_Light              = intensityTraceLight
                                IsotopicPatternMz_Heavy                     = [||]
                                IsotopicPatternIntensity_Observed_Heavy     = [||]
                                IsotopicPatternIntensity_Corrected_Heavy    = [||]
                                RtTrace_Heavy                               = [||]
                                IntensityTrace_Observed_Heavy               = [||]
                                IntensityTrace_Corrected_Heavy              = [||]
                                AlignmentScore                              = nan
                                AlignmentQValue                             = 0.
                                IonMobility                                 = targetQuant.WeightedIonMobility
                                IonMobilityCentre                           = nan
                            }
                            |> Some
            with
            | ex ->
                logger.Trace (sprintf "Gabor quantfailed: %A" ex)
                None
        let quantificationFunctions : QuantificationFunctions =
            {
                Logger = logger
                DiagCharts = diagCharts
                PlotDirectory = plotDirectory
                ScanTimeWindow = processParams.XicExtraction.ScanTimeWindow
                PeptideLookUp = peptideLookUp
                CountMatchedMasses = countMatchedMasses
                CentrePSMs = processParams.CentrePSMs
                RequirePeakAtScanTime = processParams.RequirePeakAtScanTime
                CorrectionAt = calibration.CorrectionAt
                IonMobilityOf = ionMobilityOf
                GetTargetXic = getTargetXic
                GetXicAtMobility = getXIC
                IdentifyPeaks = identifyPeaks
                ComparePredictedAndMeasuredIsotopicCluster = comparePredictedAndMeasuredIsotopicCluster
            }
        logger.Trace "init quantification functions:finished"
        
        logger.Trace "executing quantification"
        //group after peptide ion
        // -> all PSMs that have the same modifiction peptide with the same charge are getting grouped
        let quantResults = 
            qpsms 
            |> Array.groupBy (fun x -> 
                {
                    Sequence             = x.StringSequence     
                    GlobalMod            = x.GlobalMod    
                    Charge               = x.Charge       
                    ModSequenceID        = x.ModSequenceID
                    PepSequenceID        = x.PepSequenceID
                }
                )
            |> Array.map (fun (pepIon,psms) -> 
                    match processParams.XicExtraction.TopKPSMs with 
                    | Some x when psms.Length > x -> 
                        pepIon,
                        psms
                        // NaN sorts last
                        |> Array.sortByDescending (fun x -> x.ModelScore)
                        |> Array.take x
                    | _ -> 
                        pepIon,
                        psms                  
                )
            |> Xic.mapiOrderedBy (fun (_, psms: PSMStatisticsResult []) -> if psms.Length = 0 then infinity else psms |> Array.averageBy (fun x -> x.ScanTime)) (fun i (pepIon,psms) -> 
                if i % 100 = 0 then logger.Trace (sprintf "%i peptides quantified" i)
                // check if unlabeled or labeled quantification should be performed and if gabor processing is requested
                match processParams.XicExtraction.XicProcessing, processParams.PerformLabeledQuantification with 
                | Domain.XicProcessing.Gabor3D gaborParams, (Domain.Labeling.N15Labeling | Domain.Labeling.N15LabelingOnly | Domain.Labeling.Labelshift) ->
                    gaborLabeledQuantification gaborParams pepIon psms
                | Domain.XicProcessing.Gabor3D gaborParams, Domain.Labeling.Unlabeled ->
                    gaborLabelFreeQuantification gaborParams pepIon psms
                | _, (Domain.Labeling.N15Labeling | Domain.Labeling.N15LabelingOnly | Domain.Labeling.Labelshift) ->
                    quantifyLabeled quantificationFunctions pepIon psms
                | _, Domain.Labeling.Unlabeled ->
                    quantifyLabelFree quantificationFunctions pepIon psms
                )
            |> Array.choose id
                    
        let filteredResults = labelingQualityFilter processParams.PerformLabeledQuantification processParams.QualityFilterLog2Limit quantResults
        filteredResults 
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        if zipCharts then
            plotDirectory
            |> Zipping.zipDirectory "*.html" logger
            |> fun zipped ->
                match zipped with
                | Ok byteArr ->
                    byteArr
                    |> Zipping.saveZippedDirectory outputDir logger ((Path.GetFileNameWithoutExtension instrumentOutput) + "_plots")
                    |> fun saved ->
                        match saved with
                        | Ok save -> 
                            save
                            Directory.Delete (plotDirectory, true)
                        | Error ex -> logger.Trace (sprintf "Error saving zipped directory: %A" ex)
                | Error ex -> logger.Trace (sprintf "Error zipping directory: %A" ex)
        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()
        if peakCache.Reads > 0 then logger.Trace (sprintf "%i MS1 spectra were read for %i MS1s of the file (cache of %i spectra)" peakCache.Reads ms1SortedByScanTime.Length peakCache.Capacity)
        logger.Trace (sprintf "Isotope patterns: %i read from the peptide data base, %i computed" isotopePatterns.Read isotopePatterns.Computed)
        logger.Trace (sprintf "Ion mobility: %s" (centreTally.Summary ionMobilityHalfWindow))
        logger.Trace "executing quantification:finished"
        
