(*
    Temporary ProteomIQon / Gabor3D integration sketch.

    This file is intentionally NOT referenced by PSMBasedQuantificationTIMs.fsproj.
    It is a blueprint for the future implementation in PSMBasedQuantificationTIMs.fs.

    The important design choice:

        Wavelet / SecondDerivative:
            1D RT trace -> 1D peak detection -> HULQ area

        Gabor3D:
            2D RT x ion-mobility image -> Gabor response -> 2D region -> 2D volume

    This means the final output table can stay the same, but internally Gabor does
    not use HULQ.quantifyPeak. Gabor gets its own 2D quantification path.
*)

module PSMBasedQuantificationTIMsGabor3DIntegrationSketch =

    open System
    open System.Collections.Generic
    open System.Linq
    open Gabor3D.algorithm
    open MzIO.Binary
    open MzIO.Commons.Arrays
    open MzIO.Processing.MzIOLinq
    open ProteomIQon

    [<Literal>]
    let rtBinWidth = 0.01

    [<Literal>]
    let mobilityBinWidth = 0.0005

    [<Literal>]
    let gaborRelativeRegionThreshold = 0.25

    [<Literal>]
    let gaborSeedRtSearchRadiusBins = 50

    [<Literal>]
    let gaborSeedMobilitySearchRadiusBins = 20

    [<Literal>]
    let gaborMaxRtRegionRadiusBins = 120

    [<Literal>]
    let gaborMaxMobilityRegionRadiusBins = 60

    [<Literal>]
    let gaborMinRegionPixels = 3

    module Query =

        (*
            Gabor extraction does not need an ion-mobility range.

            This mirrors the logic in /home/paulinehans/Dokumente/Repos/Gabor3D/src/2D.fsx:
            - select by RT window
            - select by mz window
            - keep all peaks that actually have IonMobility
            - build the ion-mobility axis from the data itself

            An IM range can still be added later as an optional noise/performance
            filter, but it should not be required for the core Gabor path.
        *)
        let initRTProfileGabor
            (readSpecPeaks: string -> Peak1DArray)
            (rtIndex: IMzIOArray<RtIndexEntry>)
            (rtRange: RangeQuery)
            (mzRange: RangeQuery)
            =

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
                            | Some _ ->
                                Some (RtIndexEntry.AsPeak2D(peak, entry.Rt))
                            | None ->
                                None
                        )
                        |> Array.ofSeq
            |]

    type GaborObservation =
        {
            Rt: float
            Mobility: float
            Intensity: float
        }

    type GaborGrid =
        {
            RtOrigin: float
            MobilityOrigin: float
            RtBinCount: int
            MobilityBinCount: int

            // rows = ion mobility, cols = retention time
            Intensities: float[][]
        }

    type GaborPixel =
        {
            MobilityBin: int
            RtBin: int
        }

    type Gabor2DQuantification =
        {
            RawIntensitySum: float
            Volume: float
            ApexIntensity: float
            ApexRetentionTime: float
            ApexIonMobility: float
            WeightedRetentionTime: float
            WeightedIonMobility: float
            Region: GaborPixel[]
            ProjectedRetentionTimes: float[]
            ProjectedIntensities: float[]
            Seed: GaborPixel
            GaborSeedScore: float
            GaborThreshold: float
            Response: Create.gaborResponse
        }

    type GaborTargetEstimate =
        {
            ExpectedRt: float
            ExpectedMz: float
            ExpectedMobility: float
        }

    let private isUsableIntensity value =
        Double.IsFinite value && value > 0.0

    let private weightedMeanOrMean (weights: float[]) (values: float[]) =
        let weightSum =
            Array.sum weights

        if weightSum > 0.0 then
            Array.map2 (*) weights values
            |> Array.sum
            |> fun weightedSum -> weightedSum / weightSum
        else
            Array.average values

    let estimateGaborTargetFromFragPipePsms
        scanTimeToMzCorrection
        theoMz
        (psms: (Dto.PSMStatisticsResultFragpipe * float)[])
        =

        if Array.isEmpty psms then
            None
        else
            let bestPsms =
                psms
                |> Array.sortByDescending snd
                |> fun sorted ->
                    if sorted.Length > 3 then sorted.[..2] else sorted

            let weights =
                bestPsms |> Array.map snd

            let scanTimes =
                bestPsms |> Array.map (fun (psm, _) -> psm.ScanTime)

            let mobilities =
                bestPsms |> Array.map (fun (psm, _) -> psm.IonMobility)

            let expectedRt =
                weightedMeanOrMean weights scanTimes

            let expectedMobility =
                weightedMeanOrMean weights mobilities

            Some
                {
                    ExpectedRt = expectedRt
                    ExpectedMz = scanTimeToMzCorrection expectedRt + theoMz
                    ExpectedMobility = expectedMobility
                }

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

    let extractGaborObservations (peakGroups: Peak2D[][]) =
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

    let tryCreateGaborGrid (peakGroups: Peak2D[][]) =
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
                Dictionary<int * int, float>()

            observations
            |> Array.iter (fun observation ->
                let rtBin =
                    toRtBin rtOrigin observation.Rt

                let mobilityBin =
                    toMobilityBin mobilityOrigin observation.Mobility

                // key = mobility row, rt column
                let key =
                    mobilityBin, rtBin

                let oldIntensity =
                    match binnedIntensities.TryGetValue key with
                    | true, value -> value
                    | false, _ -> 0.0

                binnedIntensities.[key] <- oldIntensity + observation.Intensity
            )

            let mobilityBinCount =
                binnedIntensities.Keys
                |> Seq.maxBy fst
                |> fst
                |> fun maxBin -> maxBin + 1

            let rtBinCount =
                binnedIntensities.Keys
                |> Seq.maxBy snd
                |> snd
                |> fun maxBin -> maxBin + 1

            let intensities =
                Array.init mobilityBinCount (fun _ ->
                    Array.zeroCreate<float> rtBinCount
                )

            binnedIntensities
            |> Seq.iter (fun entry ->
                let mobilityBin, rtBin =
                    entry.Key

                intensities.[mobilityBin].[rtBin] <- entry.Value
            )

            Some
                {
                    RtOrigin = rtOrigin
                    MobilityOrigin = mobilityOrigin
                    RtBinCount = rtBinCount
                    MobilityBinCount = mobilityBinCount
                    Intensities = intensities
                }

    let createGaborKernel (parameters: Domain.Gabor3DParams) =
        Create.createGaborWavelet2D
            parameters.sizeX
            parameters.sizeY
            parameters.sigmaX
            parameters.sigmaY
            parameters.frequency
            parameters.theta

    let applyGaborToGrid (parameters: Domain.Gabor3DParams) (grid: GaborGrid) =
        let kernel =
            createGaborKernel parameters

        Create.applyGaborWavelet2D grid.Intensities kernel

    let private isInsideGrid (grid: GaborGrid) pixel =
        pixel.MobilityBin >= 0
        && pixel.MobilityBin < grid.MobilityBinCount
        && pixel.RtBin >= 0
        && pixel.RtBin < grid.RtBinCount

    let private gaborMagnitudeAt (response: Create.gaborResponse) pixel =
        response.Magnitude.[pixel.MobilityBin].[pixel.RtBin]

    let private originalIntensityAt (grid: GaborGrid) pixel =
        grid.Intensities.[pixel.MobilityBin].[pixel.RtBin]

    let private expectedPositionSeed
        (grid: GaborGrid)
        expectedRt
        expectedMobility
        =

        let expectedRtBin =
            expectedRt
            |> toRtBin grid.RtOrigin
            |> clampBin grid.RtBinCount

        let expectedMobilityBin =
            expectedMobility
            |> toMobilityBin grid.MobilityOrigin
            |> clampBin grid.MobilityBinCount

        {
            MobilityBin = expectedMobilityBin
            RtBin = expectedRtBin
        }

    let private findHighestGaborSeedNearExpectedPosition
        (grid: GaborGrid)
        (response: Create.gaborResponse)
        expectedRt
        expectedMobility
        =

        let expectedSeed =
            expectedPositionSeed grid expectedRt expectedMobility

        let minRtBin =
            expectedSeed.RtBin - gaborSeedRtSearchRadiusBins
            |> clampBin grid.RtBinCount

        let maxRtBin =
            expectedSeed.RtBin + gaborSeedRtSearchRadiusBins
            |> clampBin grid.RtBinCount

        let minMobilityBin =
            expectedSeed.MobilityBin - gaborSeedMobilitySearchRadiusBins
            |> clampBin grid.MobilityBinCount

        let maxMobilityBin =
            expectedSeed.MobilityBin + gaborSeedMobilitySearchRadiusBins
            |> clampBin grid.MobilityBinCount

        seq {
            for mobilityBin in minMobilityBin .. maxMobilityBin do
                for rtBin in minRtBin .. maxRtBin do
                    yield
                        {
                            MobilityBin = mobilityBin
                            RtBin = rtBin
                        }
        }
        |> Seq.maxBy (gaborMagnitudeAt response)

    let private isInsideRegionSearchBox seed candidate =
        abs (candidate.RtBin - seed.RtBin) <= gaborMaxRtRegionRadiusBins
        && abs (candidate.MobilityBin - seed.MobilityBin) <= gaborMaxMobilityRegionRadiusBins

    let private neighbors8 pixel =
        [|
            for mobilityDelta in -1 .. 1 do
                for rtDelta in -1 .. 1 do
                    if mobilityDelta <> 0 || rtDelta <> 0 then
                        yield
                            {
                                MobilityBin = pixel.MobilityBin + mobilityDelta
                                RtBin = pixel.RtBin + rtDelta
                            }
        |]

    let findConnectedGaborRegion
        (grid: GaborGrid)
        (response: Create.gaborResponse)
        (seed: GaborPixel)
        threshold
        =

        let visited =
            HashSet<int * int>()

        let queue =
            Queue<GaborPixel>()

        let region =
            ResizeArray<GaborPixel>()

        let addIfNew pixel =
            if isInsideGrid grid pixel && isInsideRegionSearchBox seed pixel then
                let key =
                    pixel.MobilityBin, pixel.RtBin

                if visited.Add key then
                    queue.Enqueue pixel

        addIfNew seed

        while queue.Count > 0 do
            let pixel =
                queue.Dequeue()

            if gaborMagnitudeAt response pixel >= threshold then
                region.Add pixel

                pixel
                |> neighbors8
                |> Array.iter addIfNew

        region.ToArray()

    let private projectRegionToRtTrace (grid: GaborGrid) (region: GaborPixel[]) =
        if Array.isEmpty region then
            [||], [||]
        else
            let minRtBin =
                region |> Array.minBy (fun pixel -> pixel.RtBin) |> fun pixel -> pixel.RtBin

            let maxRtBin =
                region |> Array.maxBy (fun pixel -> pixel.RtBin) |> fun pixel -> pixel.RtBin

            let regionByRt =
                region
                |> Array.groupBy (fun pixel -> pixel.RtBin)
                |> dict

            [| minRtBin .. maxRtBin |]
            |> Array.map (fun rtBin ->
                let intensity =
                    match regionByRt.TryGetValue rtBin with
                    | true, pixels ->
                        pixels
                        |> Array.sumBy (originalIntensityAt grid)
                    | false, _ ->
                        0.0

                rtBinCenter grid.RtOrigin rtBin, intensity
            )
            |> Array.unzip

    let quantifyRegionFromOriginalSignal
        (grid: GaborGrid)
        (response: Create.gaborResponse)
        seed
        threshold
        (region: GaborPixel[])
        =

        if region.Length < gaborMinRegionPixels then
            None
        else
            let rawIntensitySum =
                region
                |> Array.sumBy (originalIntensityAt grid)

            if rawIntensitySum <= 0.0 then
                None
            else
                let volume =
                    rawIntensitySum * rtBinWidth * mobilityBinWidth

                let apex =
                    region
                    |> Array.maxBy (originalIntensityAt grid)

                let apexIntensity =
                    originalIntensityAt grid apex

                let weightedRetentionTime =
                    region
                    |> Array.sumBy (fun pixel ->
                        originalIntensityAt grid pixel
                        * rtBinCenter grid.RtOrigin pixel.RtBin
                    )
                    |> fun weightedSum -> weightedSum / rawIntensitySum

                let weightedIonMobility =
                    region
                    |> Array.sumBy (fun pixel ->
                        originalIntensityAt grid pixel
                        * mobilityBinCenter grid.MobilityOrigin pixel.MobilityBin
                    )
                    |> fun weightedSum -> weightedSum / rawIntensitySum

                let projectedRetentionTimes, projectedIntensities =
                    projectRegionToRtTrace grid region

                Some
                    {
                        RawIntensitySum = rawIntensitySum
                        Volume = volume
                        ApexIntensity = apexIntensity
                        ApexRetentionTime = rtBinCenter grid.RtOrigin apex.RtBin
                        ApexIonMobility = mobilityBinCenter grid.MobilityOrigin apex.MobilityBin
                        WeightedRetentionTime = weightedRetentionTime
                        WeightedIonMobility = weightedIonMobility
                        Region = region
                        ProjectedRetentionTimes = projectedRetentionTimes
                        ProjectedIntensities = projectedIntensities
                        Seed = seed
                        GaborSeedScore = gaborMagnitudeAt response seed
                        GaborThreshold = threshold
                        Response = response
                    }

    let tryQuantifyGabor2D
        (parameters: Domain.Gabor3DParams)
        expectedRt
        expectedMobility
        (peakGroups: Peak2D[][])
        =

        peakGroups
        |> tryCreateGaborGrid
        |> Option.bind (fun grid ->
            let response =
                applyGaborToGrid parameters grid

            let seed =
                findHighestGaborSeedNearExpectedPosition
                    grid
                    response
                    expectedRt
                    expectedMobility

            let threshold =
                (gaborMagnitudeAt response seed) * gaborRelativeRegionThreshold

            let region =
                findConnectedGaborRegion grid response seed threshold

            quantifyRegionFromOriginalSignal
                grid
                response
                seed
                threshold
                region
        )

    let initGetGabor2DQuantification
        getGaborPeaks
        idx
        scanTimeWindow
        mzWindow
        =

        fun
            (parameters: Domain.Gabor3DParams)
            expectedRt
            expectedMz
            expectedMobility ->

            let rtQuery =
                MzIO.Processing.Query.createRangeQuery expectedRt scanTimeWindow

            let mzQuery =
                MzIO.Processing.Query.createRangeQuery expectedMz mzWindow

            getGaborPeaks idx rtQuery mzQuery
            |> tryQuantifyGabor2D
                parameters
                expectedRt
                expectedMobility

    (*
        How this would plug into quantifyPeptides:

        let readSpecPeaksWithMem =
            FSharpAux.Memoization.memoize inReader.ReadSpectrumPeaks

        let getPeaks =
            Query.initRTProfile readSpecPeaksWithMem

        let getXIC =
            initGetProcessedXIC
                logger
                processParams.BaseLineCorrection
                getPeaks
                retTimeIdxed
                processParams.XicExtraction.ScanTimeWindow
                mzWindow
                0.05

        let identifyPeaks =
            initIdentifyPeaks processParams.XicExtraction.XicProcessing

        let getGaborPeaks =
            Query.initRTProfileGabor readSpecPeaksWithMem

        let getGabor2DQuantification =
            initGetGabor2DQuantification
                getGaborPeaks
                retTimeIdxed
                processParams.XicExtraction.ScanTimeWindow
                mzWindow

        Then inside labledQuantification / labelFreeQuantification:

        match processParams.XicExtraction.XicProcessing with
        | Domain.XicProcessing.Gabor3D gaborParams ->
            match
                estimateGaborTargetFromFragPipePsms
                    scanTimeToMzCorrection
                    theoMz
                    psmsWithMatchedSums
            with
            | Some target ->
                match
                    getGabor2DQuantification
                        gaborParams
                        target.ExpectedRt
                        target.ExpectedMz
                        target.ExpectedMobility
                with
                | Some gaborQuant ->
                    // Use gaborQuant.Volume as Quant_Light / Quant_Heavy.
                    // Use gaborQuant.ApexIntensity as MeasuredApex_Light / Heavy.
                    // Use gaborQuant.WeightedRetentionTime for RT diagnostics.
                    // Use gaborQuant.WeightedIonMobility for IonMobility output.
                    // Use gaborQuant.ProjectedRetentionTimes and ProjectedIntensities
                    // for RtTrace_* and IntensityTrace_* diagnostics.
                    //
                    // Important: this path should NOT call HULQ.quantifyPeak, because
                    // the area is no longer a 1D fitted area. It is a 2D volume.
                    ()

                | None ->
                    // Same failure behaviour as "noPeaks" in the current 1D path.
                    ()

            | None ->
                // No PSM seed available.
                ()

        | Domain.XicProcessing.SecondDerivative _
        | Domain.XicProcessing.Wavelet _ ->
            let averagePSM =
                average getXIC scanTimeToMzCorrection theoMz psmsWithMatchedSums

            let peaks =
                identifyPeaks averagePSM.X_Xic averagePSM.Y_Xic

            // Existing 1D path stays as it is:
            // let peakToQuantify = HULQ.getPeakBy peaks averagePSM.WeightedAvgScanTime
            // let quantP = HULQ.quantifyPeak peakToQuantify
            ()

        The reason the Gabor branch can avoid ionMobilityRange:
            estimateGaborTargetFromFragPipePsms gives one expected IM seed from
            the FragPipe PSM statistics values.
            Query.initRTProfileGabor extracts all IM values in the mz/RT window.
            The Gabor grid then discovers the actual 2D region by choosing the
            highest Gabor response near that seed.

        Alternative if too much noise enters:
            Add an optional IM range later, but make it a filter around the seed,
            not a required part of the basic Gabor extraction.

        For inferred heavy/light:
            - Light PSM target -> quantify light with target mz.
            - Infer heavy by calling getGabor2DQuantification at mzHeavy.
            - Use target.ExpectedMobility, or better the target
              gaborQuant.WeightedIonMobility, as expectedMobility for inferred.

        Output compatibility:
            The output DTO can stay the same, but document that Quant_* is a 2D
            RT/IM volume for Gabor3D and a 1D HULQ area for Wavelet/SecondDerivative.

        Parameter selection:
            Do not ask interactively during quantificationTIMS. The choice already
            belongs in the param JSON via:

                "XicProcessing": {
                  "Case": "Gabor3D",
                  "Fields": [{
                    "sizeX": 41,
                    "sizeY": 41,
                    "sigmaX": 8.0,
                    "sigmaY": 8.0,
                    "frequency": 2.0,
                    "theta": 0.06
                  }]
                }

            If theta is meant as degrees, convert in createGaborKernel with
            Create.degToRad parameters.theta.
    *)
