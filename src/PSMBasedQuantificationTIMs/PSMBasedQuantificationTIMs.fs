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
open SearchDB'
open Gabor3D.algorithm

module PSMBasedQuantificationTIMs =
    module Query = 
        open System
        open System.Collections.Generic
        open System.Linq
        open MzIO.Commons.Arrays
        open MzIO.Processing.MzIOLinq
        open MzIO.Binary
         
        /// Extract a rt profile for specified target mass and rt range.
        /// Mz range peak aggregation is closest lock mz.
        /// Profile array with index corresponding to continous mass spectra over rt range and mz range given.
        let initRTProfile (readspecPeaks:string -> Peak1DArray)  (rtIndex: IMzIOArray<RtIndexEntry>) (rtRange: RangeQuery) (mzRange: RangeQuery) (ionMobilityRange: RangeQuery)=
            let entries = RtIndexEntry.Search(rtIndex, rtRange)
            entries
            |>Seq.map (fun entry ->
                let peaks = (readspecPeaks entry.SpectrumID).Peaks
                let p = 
                    peaks
                    |> Seq.filter (fun x ->
                        x.IonMobility.Value < ionMobilityRange.HighValue && x.IonMobility.Value > ionMobilityRange.LowValue &&
                        x.Mz < mzRange.HighValue && x.Mz > mzRange.LowValue
                    )
                    |> fun s ->
                        if (s |> Seq.tryHead).IsNone then
                            RtIndexEntry.AsPeak2D ((new Peak1D(0., mzRange.LockValue, ionMobilityRange.LockValue)),entry.Rt)
                        else
                            s
                            |> Seq.groupBy (fun x -> x.Mz)
                            |> Seq.map (fun (mz,peaks) -> 
                                let i = peaks |> Seq.sumBy (fun x -> x.Intensity)
                                new Peak1D(i, mz)
                            )
                            |> fun x -> 
                                RtIndexEntry.ClosestMz (x, mzRange.LockValue)
                            |> fun x -> RtIndexEntry.AsPeak2D (x, entry.Rt)
                p
            )
            |> Array.ofSeq

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
                        |> Seq.toArray
            |]

    type PeptideIon = 
        {
            Sequence             : string
            GlobalMod            : int
            Charge               : int
            ModSequenceID        : int
            PepSequenceID        : int
        }
            
    type AveragePSM = {
        MeanPrecMz   : float
        MeanScanTime : float
        WeightedAvgScanTime:float
        WeightedAvgIM:float
        MeanScore   : float
        X_Xic         : float []
        Y_Xic         : float []
        Y_Xic_uncorrected: float []
        }

    ///
    let createAveragePSM meanPrecMz meanScanTime weightedAvgScanTime weightedAvgIM meanScore xXic yXic yXic_uncorrected = {
        MeanPrecMz    = meanPrecMz
        MeanScanTime  = meanScanTime
        WeightedAvgScanTime= weightedAvgScanTime
        WeightedAvgIM = weightedAvgIM
        MeanScore = meanScore
        X_Xic         = xXic
        Y_Xic         = yXic
        Y_Xic_uncorrected= yXic_uncorrected
        }


    type PeakComparison = {
        Mz: float
        MeasuredIntensity: float
        MeasuredIntensityCorrected:float
        PredictedRelFrequency:float
        }


    type ClusterComparison = {
        PeakComparisons     : PeakComparison []
        KLDiv_UnCorrected   : float
        KLDiv_Corrected     : float
        }

    ///
    let initSpline binningWindowWidth (scanTimeVsDiff: (float*float) []) = 
        let getBinIdx width scantime = int ((scantime / width))    
        let knotX,train,test = 
            let knotX,train,test = 
                scanTimeVsDiff
                |> Array.groupBy (fun (s,d) -> getBinIdx binningWindowWidth (s))
                |> Array.map (fun (binIdx,ions) -> 
                    let dataS = ions |> Array.shuffleFisherYates
                    let knotX = 
                        dataS 
                        |> Array.map fst
                        |> Array.max
                    let train,test =
                        let nTest = 
                            (float dataS.Length) * 0.9
                            |> int
                        dataS.[.. nTest], dataS.[nTest+1 ..] 
                    float knotX, train, test
                    )
                |> Array.unzip3
            knotX |> Array.sort, train |> Array.concat |> Array.sortBy fst, test |> Array.concat |> Array.sortBy fst
        let trainer lambda = 
            let train' = train 
            let test' = test 
            let fit = FSharp.Stats.Fitting.Spline.smoothingSpline train' (knotX) lambda 
            let rSquared = 
                let x,y,yHat = 
                    test'
                    |> Array.map (fun (x,y) -> 
                        x, y, fit x
                        )
                    |> Array.unzip3
                let rs = FSharp.Stats.Fitting.GoodnessOfFit.calculateDeterminationFromValue yHat y
                rs
            rSquared, fit     
        let rSquared,model = 
            [|0.01.. 0.05 .. 0.5|]
            |> Array.map trainer
            |> Array.maxBy fst          
        rSquared, model

    ///
    let getBaseLineCorrectionOffsetAt tarRT x_Xic y_Xic y_Xic_uncorrected =
        let (rt,y,yUncorr) = 
            Array.zip3 x_Xic y_Xic y_Xic_uncorrected 
            |> Array.minBy (fun (rt,y,yUncorr) -> abs (rt - tarRT))
        yUncorr - y
       
    ///
    let getClosestMs1 (ms1s: (float*MzIO.Model.MassSpectrum) []) scanTime = 
         ms1s
         |> Array.minBy (fun ms -> abs (fst ms - scanTime))
         |> snd

    ///
    let getSpec (reader:MzIO.IO.IMzIODataReader) (ms1: MzIO.Model.MassSpectrum)  =
        Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(ms1.ID).Peaks)
        |> fun (mzData,intensityData) -> PeakArray.zip mzData intensityData
    
    ///
    let lightQualityFilter lowerBorder upperBorder (quantResults:QuantificationResult[]) =
        let medianApexIntensities = 
            quantResults
            |> Array.map (fun x -> x.MeasuredApex_Light)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        let medianQuantIntensities = 
            quantResults
            |> Array.map (fun x -> x.Quant_Light)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        quantResults
        |> Array.filter (fun x -> 
            let qualR = 
                let apexNorm = x.MeasuredApex_Light / medianApexIntensities
                let quantNorm =  x.Quant_Light / medianQuantIntensities
                quantNorm / apexNorm
                |> log2
            (qualR > lowerBorder && qualR < upperBorder) || (nan.Equals(qualR) )
            ) 

    ///
    let heavyQualityFilter lowerBorder upperBorder (quantResults:QuantificationResult[]) =
        let medianApexIntensities = 
            quantResults
            |> Array.map (fun x -> x.MeasuredApex_Heavy)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        let medianQuantIntensities = 
            quantResults
            |> Array.map (fun x -> x.Quant_Heavy)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        quantResults
        |> Array.filter (fun x -> 
            let qualR = 
                let apexNorm = x.MeasuredApex_Heavy / medianApexIntensities
                let quantNorm =  x.Quant_Heavy / medianQuantIntensities
                quantNorm / apexNorm
                |> log2
            (qualR > lowerBorder && qualR < upperBorder) || (nan.Equals(qualR) )
            )
    
    /// Calculates the Kullback-Leibler divergence Dkl(p||q) from q (theory, model, description, or approximation of p) 
    /// to p (the "true" distribution of data, observations, or a precisely calculated theoretical distribution).
    let klDiv (p:float []) (q:float []) = 
        Array.fold2 (fun acc p q -> (System.Math.Log(p/q)*p) + acc ) 0. p q
     
    ///
    let substractBaseLine (logger: NLog.Logger) (baseLineParams:Domain.BaseLineCorrection) (yData:float []) =
        if yData.Length > 500 then
            yData
        else
            let baseLine = FSharp.Stats.Signal.Baseline.baselineAls' baseLineParams.MaxIterations baseLineParams.Lambda baseLineParams.P yData |> Array.ofSeq
            Array.map2 (fun y b ->
                           let c = y - b
                           if c < 0. then 0. else c
                       ) yData baseLine

    ///
    let initGetProcessedXIC logger (baseLineCorrection:Domain.BaseLineCorrection option) getPeaks idx scanTimeWindow mzWindow_Da imWindow meanScanTime meanPrecMz meanIM =
        let rtQuery = Query.createRangeQuery meanScanTime scanTimeWindow
        let mzQuery = Query.createRangeQuery meanPrecMz mzWindow_Da
        let imQuery = Query.createRangeQuery meanIM imWindow
        let retData',itzData' =
            let tmp =
                getPeaks idx rtQuery mzQuery imQuery
                |> Array.map (fun (p:MzIO.Binary.Peak2D) -> p.Rt , p.Intensity)
            tmp
            |> Array.mapi (fun i (rt,intensity) ->
                            if i = 0 || i = tmp.Length-1 || intensity > 0. then
                                Some (rt,intensity)
                            else
                                let rt',intensity' = tmp.[i-1]
                                if intensity' = 0. then
                                    Some (rt,intensity)
                                elif intensity' > (100. * (intensity+1.)) then
                                    None
                                else
                                    Some (rt,intensity)
                          )
            |> Array.choose id
            |> Array.unzip
        match baseLineCorrection with
        | Some baseLineParams ->
            retData', substractBaseLine logger baseLineParams itzData', itzData'
        | None ->
            retData',itzData',itzData'


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
            Array.map2 (fun weight value -> weight * value) weights values
            |> Array.sum
            |> fun weightedSum -> weightedSum / weightSum
        else
            Array.average values

    let estimateGaborTargetFromFragPipePsms
        scanTimeToMzCorrection
        theoMz
        (psms: (PSMStatisticsResultFragpipe * float)[])
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
            System.Collections.Generic.HashSet<int * int>()

        let queue =
            System.Collections.Generic.Queue<GaborPixel>()

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
                System.Collections.Generic.Dictionary<int, GaborPixel[]>()

            region
            |> Array.groupBy (fun pixel -> pixel.RtBin)
            |> Array.iter (fun (rtBin, pixels) ->
                regionByRt.[rtBin] <- pixels
            )

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
                    }

    let tryQuantifyGabor2D
        (parameters: Domain.Gabor3DParams)
        expectedRt
        expectedMobility
        (peakGroups: MzIO.Binary.Peak2D[][])
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
                Query.createRangeQuery expectedRt scanTimeWindow

            let mzQuery =
                Query.createRangeQuery expectedMz mzWindow

            getGaborPeaks idx rtQuery mzQuery
            |> tryQuantifyGabor2D
                parameters
                expectedRt
                expectedMobility
        
    ///
    let initGetIsotopicEnvelope reader idx scanTimeWindow mzWindow_Da ch meanScanTime meanPrecMz =
        let rtQuery   = Query.createRangeQuery meanScanTime scanTimeWindow
        let mzQueries = 
            [|
                for i = 0 to 2 do 
                    let mz = meanPrecMz + (float i) * (Mass.Table.PMassInU / (ch|> float)) 
                    Query.createRangeQuery mz mzWindow_Da 
            |]
            |> Array.filter (fun x ->  abs (meanPrecMz-x.LockValue) < 1.)                 
        let retData',itzData' =
            let query = Query.getXICs reader idx rtQuery mzQueries 
            [|
                for i = 0 to query.[*,0].Length-1 do 
                let tmp = query.[i,*]
                yield (tmp.[0].Rt,tmp |> Array.sumBy (fun p -> p.Intensity))

            |]
            |> Array.unzip 
        retData',itzData'

    ///
    let weightedMean (weights:seq<float>) (items:seq<float>) =
        let sum,n = Seq.fold2 (fun (sum,n) w i -> w*i+sum,n + w ) (0.,0.) weights items
        sum / n
        
    ///
    let average getXic scanTimeToMzCorrection theoMz (psms:(PSMStatisticsResultFragpipe*float) []) =
            //let meanPrecMz   = psms |> Seq.meanBy (fun (psm,m) -> psm.PrecursorMZ)
            
            let meanScanTime = psms |> Seq.meanBy (fun (psm,m) -> psm.ScanTime)
            let meanScore = psms |> Seq.averageBy (fun (psm,m) -> psm.Hyperscore)
            let psms' = 
                let tmp = Array.sortByDescending (fun (psm,m) -> m) psms
                if tmp.Length > 3 then tmp.[..2] else tmp 
            let weightedAvgScanTime =
                let scanTimes = psms' |> Array.map (fun (psm,m) -> psm.ScanTime)
                let weights = psms' |> Array.map snd
                weightedMean weights scanTimes
            let weightedAvgIM =
                let ims = psms' |> Array.map (fun (psm,m) -> psm.IonMobility)
                let weights = psms' |> Array.map snd
                weightedMean weights ims
            let correctedMz = scanTimeToMzCorrection weightedAvgScanTime + theoMz
            let (retData,itzDataCorrected,ItzDataUncorrected) = getXic weightedAvgScanTime correctedMz weightedAvgIM
            createAveragePSM correctedMz meanScanTime weightedAvgScanTime weightedAvgIM meanScore retData itzDataCorrected ItzDataUncorrected



    type InferredXic = {
        X_Xic                       :float[]
        Y_Xic                       :float[]
        Y_Xic_uncorrected           :float[]
        }
    
    let getInferredXic getXic targetScanTime targetMz targetIM =
        let (retData,itzData,uncorrectedItzData)   =
                getXic targetScanTime targetMz targetIM
        {
        X_Xic               = retData
        Y_Xic               = itzData
        Y_Xic_uncorrected   = uncorrectedItzData
        }

    type InferredQuantification = {
        Model                       :HULQ.PeakModel option
        Area                        :float
        StandardErrorOfPrediction   :float
        MeasuredApexIntensity       :float
        Correlation_Light_Heavy     :float
        SearchRTMinusFittedRT       :float
        ClusterComparison           :ClusterComparison               
        EstimatedParams             :float[]
        X_Xic                       :float[]
        Y_Xic                       :float[]
        Y_Xic_uncorrected           :float[]
        xPeak                       :float[]
        yFitted                     :float[]
        }

    /// Calculates the difference between the scan time used to retreave the peak and the fitted peak midpoint.
    let searchRTMinusFittedRtTarget searchRT (fit:HULQ.QuantifiedPeak) = 
        try
            searchRT - fit.EstimatedParams.[1] 
        with
        | _ -> nan

    /// Calculates the difference between the scan time used to retreave the peak and the fitted peak midpoint.
    let searchRTMinusFittedRtInferred searchRT fit = 
        try
            searchRT - fit.EstimatedParams.[1] 
        with
        | _ -> nan
    
    /// Aims to provide the best scan time estimate when quanitfying a inferred peak.
    let chooseScanTime maxDiff searchRTMinusFittedRT initialScanTime (quantP:HULQ.QuantifiedPeak) = 
        if Array.isEmpty quantP.EstimatedParams then
            initialScanTime 
        elif abs searchRTMinusFittedRT > maxDiff then
            initialScanTime 
        else
            quantP.EstimatedParams.[1]

    ///
    let saveChart sequence globalMod ch (xXic:float[]) (yXic:float[]) ms2s avgScanTime (xToQuantify:float[]) (ypToQuantify:float[]) (fitY:float[])
            (xXicInferred:float[]) (yXicinferred:float[]) (xInferred:float[]) (inferredFit:float[]) (*(xEnvelopeSum:float[]) (yEnvelopeSum:float[])*) (peaks:FSharp.Stats.Signal.PeakDetection.IdentifiedPeak []) (pattern:PeakComparison []) plotDirectory =
        let xic = 
            [
            Chart.Point(xXic, yXic)                     |> Chart.withTraceName "Target XIC"
            peaks
            |> Array.map (fun x -> Chart.Point(x.XData,x.YData)) 
            |> Chart.Combine
            Chart.Point(ms2s)                           |> Chart.withTraceName "MS2s with scores"
            Chart.Point([avgScanTime],[1.])             |> Chart.withTraceName "Weighted Mean of Ms2 scan times"
            Chart.Point((xToQuantify), (ypToQuantify))  |> Chart.withTraceName "Identified Target Peak"
            Chart.Line(xToQuantify,fitY)                |> Chart.withTraceName "Fit of target Peak"
            Chart.Point(xXicInferred, yXicinferred)     |> Chart.withTraceName "Inferred XIC"
            Chart.Line(xInferred,inferredFit)           |> Chart.withTraceName "Fit of inferred Peak"
            //Chart.Point(xEnvelopeSum, yEnvelopeSum)     |> Chart.withTraceName "Target Envelope Sum"
            ]
            |> Chart.Combine
        let pattern = 
            [
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensity))          |> Chart.withTraceName "Measured"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensityCorrected)) |> Chart.withTraceName "Measured Corrected"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.PredictedRelFrequency))      |> Chart.withTraceName "Predicted Relative Frequency"
            //Chart.Point(xEnvelopeSum, yEnvelopeSum)     |> Chart.withTraceName "Target Envelope Sum"
            ]
            |> Chart.Combine
        [xic;pattern]
        |> Chart.Stack(2, 0.1)
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i" sequence globalMod)
        |> Chart.withSize(2500.,800.)
        |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + globalMod.ToString() + "Ch" + ch.ToString())|])

    let saveErrorChart (xXic:float[]) (yXic:float[]) pepIon desc plotDirectory =      
        Chart.Point(xXic, yXic)
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i_%s" pepIon.Sequence pepIon.GlobalMod desc)
        |> Chart.withSize(1500.,800.)
        |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((pepIon.Sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + pepIon.GlobalMod.ToString() + "Ch" + pepIon.Charge.ToString() + "_notQuantified")|])
            
    // Method is based on: https://doi.org/10.1021/ac0600196
    /// Estimates the autocorrelation at lag 1 of a blank signal (containing only noise). Subsequently, the signal of interest is smoothed
    /// several times by a savitzky golay filter using constant polynomial order and variing windowWidth. For each iteration, the deviation
    /// of the smoothed to the original signal is computed and the autocorrelation at lag 1 of this residual noise is computed. The function returns the optimized
    /// window width yielding a autocorrelation at lag 1 closest to the value computed for the blank signal.
    let optimizeWindowWidth polOrder (windowWidthToTest:int[]) noiseAutoCorr (signalOfInterest:float[]) =
        let signalOfInterest' = signalOfInterest |> vector
        //let noiseAutoCorr = Correlation.Vector.autoCorrelation 1 (blankSignal |> vector)
        let filterF w yData = FSharp.Stats.Signal.Filtering.savitzkyGolay w polOrder 0 0 yData |> vector
        let windowWidthToTest' = windowWidthToTest |> Array.filter (fun x -> x%2 <> 0)
        let optimizedWindowWidth =
            windowWidthToTest'
            |> Array.map (fun w ->
                let smoothedY = filterF w signalOfInterest
                let noise = smoothedY - (signalOfInterest')
                w, Correlation.Vector.autoCorrelation 1 noise
                )
            |> Array.minBy (fun (w,ac) -> (ac - noiseAutoCorr) |> abs )
            |> fst
        optimizedWindowWidth
    
    ///
    let initGetWindowWidth (windowEst:Domain.WindowSize) polynomOrder (windowWidthToTest:int[]) =
        match windowEst with
        | Domain.WindowSize.Fixed w  -> fun yData -> w
        | Domain.WindowSize.EstimateUsingAutoCorrelation noiseAutoCorr -> fun yData -> optimizeWindowWidth polynomOrder windowWidthToTest noiseAutoCorr yData

    ///
    let initIdentifyPeaks (peakDetectionParams:Domain.XicProcessing) =
        match peakDetectionParams with 
        | Domain.XicProcessing.SecondDerivative parameters ->
            let getWindowWidth = initGetWindowWidth parameters.WindowSize parameters.PolynomOrder [|5 .. 2 .. 60|] 
            (fun xData yData -> 
                try
                let  windowSize = getWindowWidth yData
                FSharp.Stats.Signal.PeakDetection.SecondDerivative.getPeaks parameters.MinSNR parameters.PolynomOrder windowSize xData yData
                with 
                | ex ->
                    // logger.Trace (sprintf "Quant failed: Peak detection failed with: %A" ex)
                    [||]
                )
        | Domain.XicProcessing.Wavelet parameters ->
            (fun xData yData -> 
                try
                FSharpStats'.Wavelet.identify parameters xData yData
                with 
                | ex ->
                    // logger.Trace (sprintf "Quant failed: Peak detection failed with: %A" ex)
                    [||]
                )
        | Domain.XicProcessing.Gabor3D _ ->
            fun _ _ -> [||]
        
    ///
    let calcCorrelation (xValues:float []) (quantifiedPeak:HULQ.QuantifiedPeak) (inferredPeak:HULQ.QuantifiedPeak) = 
        let getValue (model:HULQ.PeakModel) estParams =
            match model with
            | HULQ.PeakModel.Gaussian m -> 
                m.GetFunctionValue (vector estParams)
            | HULQ.PeakModel.EMG m -> 
                m.GetFunctionValue (vector estParams)
        match quantifiedPeak.Model, inferredPeak.Model with 
        | Some q , Some i -> 
            let xValuesBW = 
                [|for i = 1 to xValues.Length-1 do abs(xValues.[i] - xValues.[i-1])|]
                |> Array.median
            let xValues = [|xValues.[0] .. xValuesBW/2. .. xValues.[xValues.Length-1]|]
            let fQ = getValue q quantifiedPeak.EstimatedParams
            let yQ = xValues |> Array.map fQ
            let fI = getValue i inferredPeak.EstimatedParams
            let yI = xValues |> Array.map fI
            FSharp.Stats.Correlation.Seq.pearson yQ yI
        | _ -> nan

    ///Predicts an isotopic distribution of the given formula at the given charge, normalized by the sum of probabilities, using the MIDAs algorithm
    let generateIsotopicDistributionOfFormulaBySum (charge:int) (seq:AminoAcids.AminoAcid list) =
        seq
        |> BioList.toFormula
        |> Formula.add Formula.Table.H2O
        |> IsotopicDistribution.MIDA.ofFormula IsotopicDistribution.MIDA.normalizeByProbSum 0.01 0.001 charge 
        |> Array.ofList

    ///
    let initComparePredictedAndMeasuredIsotopicCluster inReader ms1s ms1AccuracyEstimate (x_Xic:float[]) (y_Xic:float[]) y_Xic_uncorrected ch peptideSequence tarRt tarMz =    
        /// IsotopicCluster
        let targetIsotopicPattern_predicted = 
            generateIsotopicDistributionOfFormulaBySum ch peptideSequence
        let baseLineCorrectionF = getBaseLineCorrectionOffsetAt tarRt x_Xic y_Xic y_Xic_uncorrected
        let closestMS1 = getClosestMs1 ms1s tarRt
        let peaks' = 
            getSpec inReader closestMS1
            |> Array.filter (fun x -> x.Mz < tarMz + 1. && x.Mz > tarMz - 0.6)
        let recordedVsPredictedPattern = 
            targetIsotopicPattern_predicted
            |> Array.choose (fun (mz,relFreq) ->
                if peaks' |> Array.isEmpty then None
                else
                    let closestRealPeak = peaks' |> Array.minBy (fun peak -> abs(peak.Mz - mz)) 
                    if (abs(closestRealPeak.Mz - mz) < 4. * ms1AccuracyEstimate) then  
                        Some (closestRealPeak,relFreq)
                    else None
                )
            |> Array.groupBy fst
            |> Array.map (fun ((peak),list) -> 
                let (mz,measuredIntensity,predictedRelFrequency) = peak.Mz,peak.Intensity,list |> Array.sumBy snd 
                {Mz=mz;MeasuredIntensity=measuredIntensity;MeasuredIntensityCorrected=measuredIntensity - baseLineCorrectionF;PredictedRelFrequency= predictedRelFrequency}
                )
            |> Array.filter (fun (isoP:PeakComparison) -> isoP.MeasuredIntensityCorrected > 0.)
        let recordedVsPredictedPatternNorm = 
            let sumMeasured = recordedVsPredictedPattern |> Array.sumBy (fun x -> x.MeasuredIntensity)
            let sumMeasuredCorr = recordedVsPredictedPattern |> Array.sumBy (fun x -> x.MeasuredIntensityCorrected)
            let sumMeasuredPred = recordedVsPredictedPattern |> Array.sumBy (fun x -> x.PredictedRelFrequency)
            recordedVsPredictedPattern
            |>Array.map (fun isoP -> 
                {isoP with 
                    MeasuredIntensity = isoP.MeasuredIntensity / sumMeasured; 
                    MeasuredIntensityCorrected = isoP.MeasuredIntensityCorrected / sumMeasuredCorr; 
                    PredictedRelFrequency = isoP.PredictedRelFrequency / sumMeasuredPred }
                )
        let klUnCorr = klDiv (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.MeasuredIntensity)) (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.PredictedRelFrequency))
        let klCorr   = klDiv (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.MeasuredIntensityCorrected)) (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.PredictedRelFrequency))
        {
            PeakComparisons     = recordedVsPredictedPatternNorm
            KLDiv_UnCorrected   = klUnCorr
            KLDiv_Corrected     = klCorr
        }
       
    ///
    let quantifyPeptides diagCharts zipCharts (processParams:Domain.QuantificationParams) (outputDir:string) (cn:SQLiteConnection) (instrumentOutput:string) (scoredPSMs:string)  =
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension scoredPSMs)
        logger.Trace (sprintf "Input file: %s" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        logger.Trace (sprintf "Now performing Quantification using: %s and %s, Results will be written to: %s" instrumentOutput scoredPSMs outputDir)
        // initialize Reader and Transaction
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".quant"
            Path.Combine [|outputDir;fileName|]
        logger.Trace (sprintf "outFilePath:%s" outFilePath)
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
        logger.Trace "Copy peptide DB into Memory"
        let memoryDB = SearchDB.copyDBIntoMemory cn
        logger.Trace "Copy peptide DB into Memory: finished"
        logger.Trace "Get peptide lookUp function"
        let dBParams     = getSDBParams memoryDB
        //let massLookUp = prepareSelectMassByModSequenceAndGlobalMod memoryDB
        let peptideLookUp = getThreadSafePeptideLookUpFromFileBySequenceAndGMod memoryDB dBParams
        let calcIonSeries aal = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal
        logger.Trace "Get peptide lookUp function: finished"
        // initialize Reader and Transaction
        logger.Trace "Init connection to mass spectrum data."
        let inReader = Core.MzIO.Reader.getReader instrumentOutput :?> MzIO.MzSQL.MzSQL
        inReader.Connection.Open()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader       
        let inTr = inReader.BeginTransaction()
        let countMatchedMasses (peptide: LookUpResult<AminoAcids.AminoAcid>)(psms: PSMStatisticsResultFragpipe []) =
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
        logger.Trace "Create RetentionTime index"
        let retTimeIdxed = Query.getMS1RTIdx inReader inRunID
        logger.Trace "Create RetentionTime index:finished"
        
        logger.Trace "Read and sort ms1s"
        /// Returns a Sequence containing all MassSpectra of a single MS-run
        let massSpectra = inReader.ReadMassSpectra(inRunID)
        /// Returns a Array that contains all MS1s sorted by their scanTime
        let ms1SortedByScanTime =
            massSpectra
            |> Seq.filter (fun ms -> MassSpectrum.getMsLevel ms = 1)
            |> Seq.map (fun ms -> MassSpectrum.getScanTime ms, ms)
            |> Seq.sortBy fst
            |> Array.ofSeq
        logger.Trace "Read and sort ms1s:finished"
        
        logger.Trace "Read scored PSMs"
        ///
        let qpsms =
            Csv.CsvReader<PSMStatisticsResultFragpipe>(SchemaMode=Csv.Fill).ReadFile(scoredPSMs,'\t',false,1)
            |> Array.ofSeq
        logger.Trace "Read scored PSMs:finished"
        
        logger.Trace "Estimate precursor mz standard deviation and mz correction."
        ///
        let ms1AccuracyEstimate,scanTimeToMzCorrection = 
            let scanTimeVsDelta = 
                qpsms
                |> Array.map (fun x -> 
                    let precMz = x.PrecursorMZ
                    let theMz  = Mass.toMZ x.TheoMass (float x.Charge)
                    let diff = precMz - theMz 
                    x.ScanTime,diff
                    )
            let borders =
                scanTimeVsDelta
                |> Seq.map snd
                |> Array.ofSeq
                |> FSharp.Stats.Testing.Outliers.tukey 3.
            let filteredValues =
                scanTimeVsDelta
                |> Array.ofSeq
                |> Array.filter (fun (s,d) -> (d <= borders.Upper && d >= borders.Lower) && d < 0.1)
            let scanTimeToMzCorrection =
                let runTime = scanTimeVsDelta |> Array.maxBy fst |> fst
                if runTime > 20. && qpsms.Length > 500 then
                    let binWidth = 
                        System.Math.Min(runTime / 2., 20.)
                    let stabw = filteredValues |> Seq.stDevBy snd
                    let r,f = initSpline binWidth filteredValues
                    f
                else 
                    let m = 
                        filteredValues 
                        |> Seq.map snd 
                        |> Seq.median 
                    fun scanTime -> m
            let stDev = 
                filteredValues 
                |> Seq.map snd 
                |> Seq.stDev
            if diagCharts then 
                [
                Chart.Point(scanTimeVsDelta)
                |> Chart.withTraceName "Raw"
                Chart.Line(scanTimeVsDelta |> Array.sortBy fst |> Array.map (fun (st,d) -> st, scanTimeToMzCorrection st))
                ]
                |> Chart.Combine
                |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; "mzErrorAndCorrection"|])
            stDev, scanTimeToMzCorrection 
        logger.Trace (sprintf "Estimate precursor mz standard deviation and mz correction.:finished, standard deviation: %f" ms1AccuracyEstimate) 

        let qpsmsMzRefined = 
            let refined = 
                qpsms
                |> Array.map (fun x -> 
                    let precMz = x.PrecursorMZ
                    let theoMz  = Mass.toMZ x.TheoMass (float x.Charge)
                    let diffPToTheo = precMz - theoMz 
                    let absDiffToPredDiff = 
                        (abs (scanTimeToMzCorrection x.ScanTime)) - (abs diffPToTheo)
                        |> abs
                    if absDiffToPredDiff > 4.*ms1AccuracyEstimate then 
                        let precMz' = theoMz + scanTimeToMzCorrection x.ScanTime
                        {x with PrecursorMZ = precMz'}
                    else    
                        x
                    )
            if diagCharts then 
                [
                qpsms
                |> Array.map (fun x -> 
                    let precMz = x.PrecursorMZ
                    let theoMz  = Mass.toMZ x.TheoMass (float x.Charge)
                    let diff = precMz - theoMz 
                    x.ScanTime,diff
                    )
                |> Chart.Point
                |> Chart.withTraceName "Raw"
                refined
                |> Array.map (fun x -> 
                    let precMz = x.PrecursorMZ
                    let theoMz  = Mass.toMZ x.TheoMass (float x.Charge)
                    let diff = precMz - theoMz 
                    x.ScanTime,diff
                    )
                |> Chart.Point
                |> Chart.withTraceName "Corrected"
                ]
                |> Chart.Combine
                |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; "precMzCorrected"|])
            refined
                
           
        logger.Trace "init lookup functions"        
        ///
        let comparePredictedAndMeasuredIsotopicCluster = initComparePredictedAndMeasuredIsotopicCluster inReader ms1SortedByScanTime ms1AccuracyEstimate     
        
        let mzWindow = 
            match processParams.XicExtraction.MzWindow_Da with 
            | Domain.Window.Fixed v -> v 
            | Domain.Window.Estimate -> 
                let mzW = ms1AccuracyEstimate*4.
                logger.Trace (sprintf "optimal mz Window for XIC look up found by estimation :%f Da" mzW)  
                mzW

        let readSpecPeaksWithMem = FSharpAux.Memoization.memoize inReader.ReadSpectrumPeaks
        ///
        let getPeaks = Query.initRTProfile readSpecPeaksWithMem
        let getXIC = initGetProcessedXIC logger processParams.BaseLineCorrection getPeaks retTimeIdxed processParams.XicExtraction.ScanTimeWindow mzWindow 0.05        
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
            peptideSequence
            targetMz
            =

            try
                if Array.isEmpty quant.ProjectedRetentionTimes || Array.isEmpty quant.ProjectedIntensities then
                    emptyClusterComparison
                else
                    comparePredictedAndMeasuredIsotopicCluster
                        quant.ProjectedRetentionTimes
                        quant.ProjectedIntensities
                        quant.ProjectedIntensities
                        charge
                        peptideSequence
                        quant.WeightedRetentionTime
                        targetMz
            with
            | ex ->
                logger.Trace (sprintf "Gabor isotope cluster comparison failed: %A" ex)
                emptyClusterComparison

        let quantMzAtScanTime (peptide: LookUpResult<AminoAcids.AminoAcid>) charge scanTime =
            scanTimeToMzCorrection scanTime + Mass.toMZ peptide.Mass (float charge)

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
                    compareGaborIsotopicCluster q charge peptide.BioSequence quantMz

                q.Volume,
                q.ApexIntensity,
                nan,
                ([||] : float[]),
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
            (psms: PSMStatisticsResultFragpipe[])
            =

            try
                let bestQValue, prots =
                    psms
                    |> Array.minBy (fun x -> x.Expectscore)
                    |> fun x -> x.Expectscore, x.ProteinNames

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

                match estimateGaborTargetFromFragPipePsms scanTimeToMzCorrection theoMz psmsWithMatchedSums with
                | None ->
                    logger.Trace (sprintf "Gabor quant failed: No FragPipe seed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                    None
                | Some target ->
                    match getGabor2DQuantification gaborParams target.ExpectedRt target.ExpectedMz target.ExpectedMobility with
                    | None ->
                        logger.Trace (sprintf "Gabor quant failed: No 2D region, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    | Some targetQuant ->
                        let meanScore =
                            psmsWithMatchedSums
                            |> Array.averageBy (fun (psm, _) -> psm.Hyperscore)

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
                                PEPValue                                    = nan
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
                                AlignmentQValue                             = nan
                                IonMobility                                 = targetQuant.WeightedIonMobility
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
                                PEPValue                                    = nan
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
                                AlignmentQValue                             = nan
                                IonMobility                                 = targetQuant.WeightedIonMobility
                            }
                            |> Some
            with
            | ex ->
                logger.Trace (sprintf "Gabor quantfailed: %A" ex)
                None

        let gaborLabelFreeQuantification
            (gaborParams: Domain.Gabor3DParams)
            (pepIon: PeptideIon)
            (psms: PSMStatisticsResultFragpipe[])
            =

            try
                if pepIon.GlobalMod <> 0 then
                    None
                else
                    let bestQValue, prots =
                        psms
                        |> Array.minBy (fun x -> x.Expectscore)
                        |> fun x -> x.Expectscore, x.ProteinNames

                    let unlabledPeptide =
                        peptideLookUp pepIon.Sequence 0

                    let psmsWithMatchedSums =
                        countMatchedMasses unlabledPeptide psms

                    let theoMz =
                        Mass.toMZ unlabledPeptide.Mass (float pepIon.Charge)

                    match estimateGaborTargetFromFragPipePsms scanTimeToMzCorrection theoMz psmsWithMatchedSums with
                    | None ->
                        logger.Trace (sprintf "Gabor quant failed: No FragPipe seed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    | Some target ->
                        match getGabor2DQuantification gaborParams target.ExpectedRt target.ExpectedMz target.ExpectedMobility with
                        | None ->
                            logger.Trace (sprintf "Gabor quant failed: No 2D region, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                            None
                        | Some targetQuant ->
                            let meanScore =
                                psmsWithMatchedSums
                                |> Array.averageBy (fun (psm, _) -> psm.Hyperscore)

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
                                PEPValue                                    = nan
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
                                AlignmentQValue                             = nan
                                IonMobility                                 = targetQuant.WeightedIonMobility
                            }
                            |> Some
            with
            | ex ->
                logger.Trace (sprintf "Gabor quantfailed: %A" ex)
                None
        ///
        let labledQuantification (pepIon:PeptideIon) (psms:PSMStatisticsResultFragpipe []) = 
            try
            let bestQValue,prots = psms |> Array.minBy (fun x -> x.Expectscore) |> fun x -> x.Expectscore,x.ProteinNames
            let unlabledPeptide = peptideLookUp pepIon.Sequence 0
            let labeledPeptide  = peptideLookUp pepIon.Sequence 1
            let targetPeptide = if pepIon.GlobalMod = 0 then unlabledPeptide else labeledPeptide            
            let psmsWithMatchedSums = countMatchedMasses targetPeptide psms 
            let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
            let theoMz = Mass.toMZ targetPeptide.Mass (float pepIon.Charge)
            let averagePSM = average getXIC scanTimeToMzCorrection theoMz psmsWithMatchedSums
            let avgMass = Mass.ofMZ (averagePSM.MeanPrecMz) (pepIon.Charge |> float)
            let peaks = identifyPeaks averagePSM.X_Xic averagePSM.Y_Xic
            if Array.isEmpty peaks then 
                if diagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "noPeaks" plotDirectory
                logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                None
            else
            let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks averagePSM.WeightedAvgScanTime
            let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify                        
            let searchRTMinusFittedRT = searchRTMinusFittedRtTarget averagePSM.WeightedAvgScanTime quantP
            if quantP.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantP.EstimatedParams then 
                if diagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "fittingFailed" plotDirectory    
                logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                None
            else
            let inferredScanTime = chooseScanTime processParams.XicExtraction.ScanTimeWindow searchRTMinusFittedRT averagePSM.WeightedAvgScanTime quantP 
            let clusterComparisonTarget = comparePredictedAndMeasuredIsotopicCluster averagePSM.X_Xic averagePSM.Y_Xic averagePSM.Y_Xic_uncorrected pepIon.Charge targetPeptide.BioSequence quantP.EstimatedParams.[1] averagePSM.MeanPrecMz            
            if pepIon.GlobalMod = 0 then
                let mzHeavy = 
                    let mz = Mass.toMZ (labeledPeptide.Mass) (pepIon.Charge|> float)
                    let correctedMz = scanTimeToMzCorrection inferredScanTime + mz
                    correctedMz 
                let inferredQuant = 
                    let inferredXicHeavy = getInferredXic getXIC averagePSM.WeightedAvgScanTime mzHeavy averagePSM.WeightedAvgIM
                    let inferredPeaksHeavy = identifyPeaks inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic
                    if Array.isEmpty inferredPeaksHeavy then 
                        if diagCharts then saveErrorChart inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic pepIon "noInferredPeaks" plotDirectory
                        logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    else
                    let peakToQuantifyHeavy = BioFSharp.Mz.Quantification.HULQ.getPeakBy inferredPeaksHeavy inferredScanTime
                    let quantPHeavy         = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantifyHeavy          
                    if quantPHeavy.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantPHeavy.EstimatedParams then 
                        if diagCharts then saveErrorChart inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic pepIon "fittingInferredFailed" plotDirectory    
                        logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    else
                    // let inferred_Heavy = quantifyInferredPeak getXIC identifyPeaks mz_Heavy averagePSM.WeightedAvgScanTime inferredScanTime
                    let searchRTMinusFittedRTHeavy = searchRTMinusFittedRtTarget inferredScanTime quantPHeavy
                    let clusterComparisonHeavy = comparePredictedAndMeasuredIsotopicCluster inferredXicHeavy.X_Xic inferredXicHeavy.Y_Xic inferredXicHeavy.Y_Xic_uncorrected pepIon.Charge labeledPeptide.BioSequence quantP.EstimatedParams.[1] mzHeavy
                    let corrLightHeavy  = calcCorrelation averagePSM.X_Xic quantP quantPHeavy
                    {
                    Model                       = quantPHeavy.Model
                    Area                        = quantPHeavy.Area
                    StandardErrorOfPrediction   = quantPHeavy.StandardErrorOfPrediction
                    MeasuredApexIntensity       = quantPHeavy.MeasuredApexIntensity
                    Correlation_Light_Heavy     = corrLightHeavy
                    SearchRTMinusFittedRT       = searchRTMinusFittedRTHeavy 
                    ClusterComparison           = clusterComparisonHeavy     
                    EstimatedParams             = quantPHeavy.EstimatedParams
                    X_Xic                       = inferredXicHeavy.X_Xic
                    Y_Xic                       = inferredXicHeavy.Y_Xic 
                    Y_Xic_uncorrected           = inferredXicHeavy.Y_Xic_uncorrected 
                    xPeak                       = peakToQuantifyHeavy.XData 
                    yFitted                     = quantPHeavy.YPredicted
                    }
                    |> Some                      
                match inferredQuant with 
                | Some successfulQuant -> 
                    if diagCharts then 
                        saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                        peakToQuantify.XData peakToQuantify.YData quantP.YPredicted successfulQuant.X_Xic successfulQuant.Y_Xic successfulQuant.xPeak successfulQuant.yFitted peaks clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = pepIon.Sequence
                    GlobalMod                                   = pepIon.GlobalMod
                    Charge                                      = pepIon.Charge
                    PepSequenceID                               = pepIon.PepSequenceID
                    ModSequenceID                               = pepIon.ModSequenceID
                    PrecursorMZ                                 = averagePSM.MeanPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = unlabledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-unlabledPeptide.Mass)
                    MeanPercolatorScore                         = averagePSM.MeanScore
                    QValue                                      = bestQValue
                    PEPValue                                    = nan
                    ProteinNames                                = prots
                    QuantMz_Light                               = averagePSM.MeanPrecMz
                    Quant_Light                                 = quantP.Area
                    MeasuredApex_Light                          = quantP.MeasuredApexIntensity
                    Seo_Light                                   = quantP.StandardErrorOfPrediction
                    Params_Light                                = quantP.EstimatedParams            
                    Difference_SearchRT_FittedRT_Light          = searchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Light            = clusterComparisonTarget.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Light   = clusterComparisonTarget.KLDiv_Corrected
                    QuantMz_Heavy                               = mzHeavy
                    Quant_Heavy                                 = successfulQuant.Area
                    MeasuredApex_Heavy                          = successfulQuant.MeasuredApexIntensity
                    Seo_Heavy                                   = successfulQuant.StandardErrorOfPrediction
                    Params_Heavy                                = successfulQuant.EstimatedParams 
                    Difference_SearchRT_FittedRT_Heavy          = successfulQuant.SearchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Heavy            = successfulQuant.ClusterComparison.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Heavy   = successfulQuant.ClusterComparison.KLDiv_Corrected
                    Correlation_Light_Heavy                     = successfulQuant.Correlation_Light_Heavy
                    QuantificationSource                        = QuantificationSource.PSM
                    IsotopicPatternMz_Light                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Light     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Light    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Light                               = averagePSM.X_Xic 
                    IntensityTrace_Observed_Light               = averagePSM.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Light              = averagePSM.Y_Xic
                    IsotopicPatternMz_Heavy                     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Heavy     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Heavy    = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Heavy                               = successfulQuant.X_Xic 
                    IntensityTrace_Observed_Heavy               = successfulQuant.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Heavy              = successfulQuant.Y_Xic
                    AlignmentScore                              = nan
                    AlignmentQValue                             = nan
                    IonMobility                                 = averagePSM.WeightedAvgIM
                    }
                    |> Option.Some
                | None -> 
                    if diagCharts then 
                        saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                        peakToQuantify.XData peakToQuantify.YData quantP.YPredicted [||] [||] [||] [||] peaks clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = pepIon.Sequence
                    GlobalMod                                   = pepIon.GlobalMod
                    Charge                                      = pepIon.Charge
                    PepSequenceID                               = pepIon.PepSequenceID
                    ModSequenceID                               = pepIon.ModSequenceID
                    PrecursorMZ                                 = averagePSM.MeanPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = unlabledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-unlabledPeptide.Mass)
                    MeanPercolatorScore                         = averagePSM.MeanScore
                    QValue                                      = bestQValue
                    PEPValue                                    = nan
                    ProteinNames                                = prots
                    QuantMz_Light                               = averagePSM.MeanPrecMz
                    Quant_Light                                 = quantP.Area
                    MeasuredApex_Light                          = quantP.MeasuredApexIntensity
                    Seo_Light                                   = quantP.StandardErrorOfPrediction
                    Params_Light                                = quantP.EstimatedParams            
                    Difference_SearchRT_FittedRT_Light          = searchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Light            = clusterComparisonTarget.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Light   = clusterComparisonTarget.KLDiv_Corrected
                    QuantMz_Heavy                               = mzHeavy
                    Quant_Heavy                                 = nan
                    MeasuredApex_Heavy                          = nan
                    Seo_Heavy                                   = nan
                    Params_Heavy                                = [||]
                    Difference_SearchRT_FittedRT_Heavy          = nan
                    KLDiv_Observed_Theoretical_Heavy            = nan
                    KLDiv_CorrectedObserved_Theoretical_Heavy   = nan
                    Correlation_Light_Heavy                     = nan
                    QuantificationSource                        = QuantificationSource.PSM
                    IsotopicPatternMz_Light                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Light     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Light    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Light                               = averagePSM.X_Xic 
                    IntensityTrace_Observed_Light               = averagePSM.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Light              = averagePSM.Y_Xic
                    IsotopicPatternMz_Heavy                     = [||]
                    IsotopicPatternIntensity_Observed_Heavy     = [||]
                    IsotopicPatternIntensity_Corrected_Heavy    = [||]
                    RtTrace_Heavy                               = [||]
                    IntensityTrace_Observed_Heavy               = [||]
                    IntensityTrace_Corrected_Heavy              = [||]
                    AlignmentScore                              = nan
                    AlignmentQValue                             = nan
                    IonMobility                                 = averagePSM.WeightedAvgIM
                    }
                    |> Option.Some
            else
                let mzLight = 
                    let mz = Mass.toMZ (unlabledPeptide.Mass) (pepIon.Charge|> float)
                    let correctedMz = scanTimeToMzCorrection inferredScanTime + mz
                    correctedMz   
                let inferredQuant = 
                    let inferredXicLight = getInferredXic getXIC averagePSM.WeightedAvgScanTime mzLight averagePSM.WeightedAvgIM
                    let inferredPeaksLight = identifyPeaks inferredXicLight.X_Xic inferredXicLight.Y_Xic
                    if Array.isEmpty inferredPeaksLight then 
                        if diagCharts then saveErrorChart inferredXicLight.X_Xic inferredXicLight.Y_Xic pepIon "noInferredPeaks" plotDirectory
                        logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    else
                    let peakToQuantifyLight = BioFSharp.Mz.Quantification.HULQ.getPeakBy inferredPeaksLight inferredScanTime
                    let quantPLight         = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantifyLight          
                    if quantPLight.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantPLight.EstimatedParams then 
                        if diagCharts then saveErrorChart inferredXicLight.X_Xic inferredXicLight.Y_Xic pepIon "fittingInferredFailed" plotDirectory    
                        logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                        None
                    else
                    // let inferred_Heavy = quantifyInferredPeak getXIC identifyPeaks mz_Heavy averagePSM.WeightedAvgScanTime inferredScanTime
                    let searchRTMinusFittedRTLight = searchRTMinusFittedRtTarget inferredScanTime quantPLight
                    let clusterComparisonLight = comparePredictedAndMeasuredIsotopicCluster inferredXicLight.X_Xic inferredXicLight.Y_Xic inferredXicLight.Y_Xic_uncorrected pepIon.Charge unlabledPeptide.BioSequence quantP.EstimatedParams.[1] mzLight
                    let corrLightHeavy  = calcCorrelation averagePSM.X_Xic quantP quantPLight
                    {
                    Model                       = quantPLight.Model
                    Area                        = quantPLight.Area
                    StandardErrorOfPrediction   = quantPLight.StandardErrorOfPrediction
                    MeasuredApexIntensity       = quantPLight.MeasuredApexIntensity
                    Correlation_Light_Heavy     = corrLightHeavy
                    SearchRTMinusFittedRT       = searchRTMinusFittedRTLight 
                    ClusterComparison           = clusterComparisonLight     
                    EstimatedParams             = quantPLight.EstimatedParams
                    X_Xic                       = inferredXicLight.X_Xic
                    Y_Xic                       = inferredXicLight.Y_Xic 
                    Y_Xic_uncorrected           = inferredXicLight.Y_Xic_uncorrected 
                    xPeak                       = peakToQuantifyLight.XData 
                    yFitted                     = quantPLight.YPredicted
                    }
                    |> Some                      
                match inferredQuant with 
                | Some successfulQuant -> 
                    if diagCharts then 
                        saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                        peakToQuantify.XData peakToQuantify.YData quantP.YPredicted successfulQuant.X_Xic successfulQuant.Y_Xic successfulQuant.xPeak successfulQuant.yFitted peaks clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = pepIon.Sequence
                    GlobalMod                                   = pepIon.GlobalMod
                    Charge                                      = pepIon.Charge
                    PepSequenceID                               = pepIon.PepSequenceID
                    ModSequenceID                               = pepIon.ModSequenceID
                    PrecursorMZ                                 = averagePSM.MeanPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = labeledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-labeledPeptide.Mass)
                    MeanPercolatorScore                         = averagePSM.MeanScore
                    QValue                                      = bestQValue
                    PEPValue                                    = nan
                    ProteinNames                                = prots
                    QuantMz_Light                               = mzLight
                    Quant_Light                                 = successfulQuant.Area
                    MeasuredApex_Light                          = successfulQuant.MeasuredApexIntensity
                    Seo_Light                                   = successfulQuant.StandardErrorOfPrediction
                    Params_Light                                = successfulQuant.EstimatedParams 
                    Difference_SearchRT_FittedRT_Light          = successfulQuant.SearchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Light            = successfulQuant.ClusterComparison.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Light   = successfulQuant.ClusterComparison.KLDiv_Corrected
                    QuantMz_Heavy                               = averagePSM.MeanPrecMz
                    Quant_Heavy                                 = quantP.Area
                    MeasuredApex_Heavy                          = quantP.MeasuredApexIntensity
                    Seo_Heavy                                   = quantP.StandardErrorOfPrediction
                    Params_Heavy                                = quantP.EstimatedParams 
                    Difference_SearchRT_FittedRT_Heavy          = searchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Heavy            = clusterComparisonTarget.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterComparisonTarget.KLDiv_Corrected
                    Correlation_Light_Heavy                     = successfulQuant.Correlation_Light_Heavy
                    QuantificationSource                        = QuantificationSource.PSM
                    IsotopicPatternMz_Light                     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Light     = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Light    = successfulQuant.ClusterComparison.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Light                               = successfulQuant.X_Xic 
                    IntensityTrace_Observed_Light               = successfulQuant.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Light              = successfulQuant.Y_Xic
                    IsotopicPatternMz_Heavy                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Heavy     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Heavy    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Heavy                               = averagePSM.X_Xic 
                    IntensityTrace_Observed_Heavy               = averagePSM.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Heavy              = averagePSM.Y_Xic
                    AlignmentScore                              = nan
                    AlignmentQValue                             = nan
                    IonMobility                                 = averagePSM.WeightedAvgIM
                    }
                    |> Option.Some
                | None ->
                    if diagCharts then 
                        saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                            peakToQuantify.XData peakToQuantify.YData quantP.YPredicted [||] [||] [||] [||] peaks clusterComparisonTarget.PeakComparisons plotDirectory
                    {
                    StringSequence                              = pepIon.Sequence
                    GlobalMod                                   = pepIon.GlobalMod
                    Charge                                      = pepIon.Charge
                    PepSequenceID                               = pepIon.PepSequenceID
                    ModSequenceID                               = pepIon.ModSequenceID
                    PrecursorMZ                                 = averagePSM.MeanPrecMz
                    MeasuredMass                                = avgMass
                    TheoMass                                    = labeledPeptide.Mass
                    AbsDeltaMass                                = abs(avgMass-labeledPeptide.Mass)
                    MeanPercolatorScore                         = averagePSM.MeanScore
                    QValue                                      = bestQValue
                    PEPValue                                    = nan
                    ProteinNames                                = prots
                    QuantMz_Light                               = mzLight
                    Quant_Light                                 = nan
                    MeasuredApex_Light                          = nan
                    Seo_Light                                   = nan
                    Params_Light                                = [||] 
                    Difference_SearchRT_FittedRT_Light          = nan
                    KLDiv_Observed_Theoretical_Light            = nan
                    KLDiv_CorrectedObserved_Theoretical_Light   = nan
                    QuantMz_Heavy                               = averagePSM.MeanPrecMz
                    Quant_Heavy                                 = quantP.Area
                    MeasuredApex_Heavy                          = quantP.MeasuredApexIntensity
                    Seo_Heavy                                   = quantP.StandardErrorOfPrediction
                    Params_Heavy                                = quantP.EstimatedParams 
                    Difference_SearchRT_FittedRT_Heavy          = searchRTMinusFittedRT
                    KLDiv_Observed_Theoretical_Heavy            = clusterComparisonTarget.KLDiv_UnCorrected
                    KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterComparisonTarget.KLDiv_Corrected
                    Correlation_Light_Heavy                     = nan
                    QuantificationSource                        = QuantificationSource.PSM
                    IsotopicPatternMz_Light                     = [||] 
                    IsotopicPatternIntensity_Observed_Light     = [||] 
                    IsotopicPatternIntensity_Corrected_Light    = [||] 
                    RtTrace_Light                               = [||] 
                    IntensityTrace_Observed_Light               = [||] 
                    IntensityTrace_Corrected_Light              = [||] 
                    IsotopicPatternMz_Heavy                     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.Mz)
                    IsotopicPatternIntensity_Observed_Heavy     = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                    IsotopicPatternIntensity_Corrected_Heavy    = clusterComparisonTarget.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                    RtTrace_Heavy                               = averagePSM.X_Xic 
                    IntensityTrace_Observed_Heavy               = averagePSM.Y_Xic_uncorrected
                    IntensityTrace_Corrected_Heavy              = averagePSM.Y_Xic
                    AlignmentScore                              = nan
                    AlignmentQValue                             = nan
                    IonMobility                                 = averagePSM.WeightedAvgIM
                    }
                    |> Option.Some
            with
            | ex ->             
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None

        let lableFreeQuantification (pepIon:PeptideIon)  (psms:PSMStatisticsResultFragpipe []) = 
            try
            if pepIon.GlobalMod <> 0 then None 
            else
            let bestQValue,prots = psms |> Array.minBy (fun x -> x.Expectscore) |> fun x -> x.Expectscore,x.ProteinNames
            let unlabledPeptide = peptideLookUp pepIon.Sequence 0
            let psmsWithMatchedSums = countMatchedMasses unlabledPeptide psms 
            let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
            let theoMz = Mass.toMZ unlabledPeptide.Mass (float pepIon.Charge)
            let averagePSM = average getXIC scanTimeToMzCorrection theoMz psmsWithMatchedSums
            let avgMass = Mass.ofMZ (averagePSM.MeanPrecMz) (pepIon.Charge |> float)      
            let peaks = identifyPeaks averagePSM.X_Xic averagePSM.Y_Xic               
            if Array.isEmpty peaks then 
                if diagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "noPeaks" plotDirectory
                logger.Trace (sprintf "Quant failed: No Peak detected, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                None
            else
            let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks averagePSM.WeightedAvgScanTime
            let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify
            let searchRTMinusFittedRT = searchRTMinusFittedRtTarget averagePSM.WeightedAvgScanTime quantP 
            if quantP.EstimatedParams |> Array.exists nan.Equals || Array.isEmpty quantP.EstimatedParams then 
                if diagCharts then saveErrorChart averagePSM.X_Xic averagePSM.Y_Xic pepIon "fittingFailed" plotDirectory    
                logger.Trace (sprintf "Quant failed: Peak fitting failed, Sequence:%s, GlobalMod:%s, Charge:%s" (pepIon.Sequence |> String.filter (fun x -> x <> '*')) (pepIon.GlobalMod.ToString()) (pepIon.Charge.ToString()))
                None
            else
            let clusterComparisonTarget = comparePredictedAndMeasuredIsotopicCluster averagePSM.X_Xic averagePSM.Y_Xic averagePSM.Y_Xic_uncorrected pepIon.Charge unlabledPeptide.BioSequence quantP.EstimatedParams.[1] averagePSM.MeanPrecMz
            if diagCharts then 
                saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                                    peakToQuantify.XData peakToQuantify.YData quantP.YPredicted  [||] [||] [||] [||] peaks clusterComparisonTarget.PeakComparisons plotDirectory
            {
            StringSequence                              = pepIon.Sequence
            GlobalMod                                   = pepIon.GlobalMod
            Charge                                      = pepIon.Charge
            PepSequenceID                               = pepIon.PepSequenceID
            ModSequenceID                               = pepIon.ModSequenceID
            PrecursorMZ                                 = averagePSM.MeanPrecMz
            MeasuredMass                                = avgMass
            TheoMass                                    = unlabledPeptide.Mass
            AbsDeltaMass                                = abs(avgMass-unlabledPeptide.Mass)
            MeanPercolatorScore                         = averagePSM.MeanScore
            QValue                                      = bestQValue
            PEPValue                                    = nan
            ProteinNames                                = prots
            QuantMz_Light                               = averagePSM.MeanPrecMz
            Quant_Light                                 = quantP.Area
            MeasuredApex_Light                          = quantP.MeasuredApexIntensity
            Seo_Light                                   = quantP.StandardErrorOfPrediction
            Params_Light                                = quantP.EstimatedParams 
            Difference_SearchRT_FittedRT_Light          = searchRTMinusFittedRT
            KLDiv_Observed_Theoretical_Light            = clusterComparisonTarget.KLDiv_UnCorrected
            KLDiv_CorrectedObserved_Theoretical_Light   = clusterComparisonTarget.KLDiv_Corrected
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
            IsotopicPatternMz_Light                     = clusterComparisonTarget.PeakComparisons|> Array.map (fun x -> x.Mz)
            IsotopicPatternIntensity_Observed_Light     = clusterComparisonTarget.PeakComparisons|> Array.map (fun x -> x.MeasuredIntensity)
            IsotopicPatternIntensity_Corrected_Light    = clusterComparisonTarget.PeakComparisons|> Array.map (fun x -> x.MeasuredIntensityCorrected)
            RtTrace_Light                               = averagePSM.X_Xic 
            IntensityTrace_Observed_Light               = averagePSM.Y_Xic_uncorrected
            IntensityTrace_Corrected_Light              = averagePSM.Y_Xic
            IsotopicPatternMz_Heavy                     = [||]
            IsotopicPatternIntensity_Observed_Heavy     = [||]
            IsotopicPatternIntensity_Corrected_Heavy    = [||]
            RtTrace_Heavy                               = [||]
            IntensityTrace_Observed_Heavy               = [||]
            IntensityTrace_Corrected_Heavy              = [||]
            AlignmentScore                              = nan
            AlignmentQValue                             = nan
            IonMobility                                 = averagePSM.WeightedAvgIM
            }
            |> Option.Some
            with
            | ex ->
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None
        logger.Trace "init quantification functions:finished"
        
        logger.Trace "executing quantification"
        let quantResults = 
            qpsmsMzRefined 
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
                        |> Array.sortByDescending (fun x -> x.Hyperscore)
                        |> Array.take x
                    | _ -> 
                        pepIon,
                        psms                  
                )
            |> Array.mapi (fun i (pepIon,psms) -> 
                if i % 100 = 0 then logger.Trace (sprintf "%i peptides quantified" i)
                
                match processParams.XicExtraction.XicProcessing, processParams.PerformLabeledQuantification with 
                | Domain.XicProcessing.Gabor3D gaborParams, (Domain.Labeling.N15Labeling | Domain.Labeling.N15LabelingOnly | Domain.Labeling.Labelshift) ->
                    gaborLabeledQuantification gaborParams pepIon psms
                | Domain.XicProcessing.Gabor3D gaborParams, Domain.Labeling.Unlabeled ->
                    gaborLabelFreeQuantification gaborParams pepIon psms
                | _, (Domain.Labeling.N15Labeling | Domain.Labeling.N15LabelingOnly | Domain.Labeling.Labelshift) ->
                    labledQuantification pepIon psms
                | _, Domain.Labeling.Unlabeled ->
                    lableFreeQuantification pepIon psms
                )
            |> Array.choose id
                    
        let filteredResults = 
            match processParams.PerformLabeledQuantification with 
            |Domain.Labeling.N15Labeling ->
                quantResults
                |> heavyQualityFilter -2. 2.
                |> lightQualityFilter -2. 2.
            |Domain.Labeling.N15LabelingOnly ->
                quantResults
                |> heavyQualityFilter -2. 2.
            |Domain.Labeling.Unlabeled ->
                quantResults
                |> lightQualityFilter -2. 2.
            |Domain.Labeling.Labelshift ->
                quantResults
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
        memoryDB.Dispose()
        logger.Trace "executing quantification:finished"
        
