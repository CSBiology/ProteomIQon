#r "nuget: MzIO, 0.1.7"
#r "nuget: MzIO.Processing, 0.1.7"
#r "nuget: MzLite, 1.1.0"
#r "nuget: MzIO.SQL, 0.1.9"
#r "nuget: Plotly.NET, 5.0.0"
#r "nuget: Plotly.NET.Interactive, 5.0.0"
#r "nuget: MzIO.MzML, 0.1.10"
#r "nuget: Deedle, 2.5.0"
#r "nuget: BioFSharp.ImgP, 2.0.0-beta5"
#r "nuget: BioFSharp.Mz, 0.1.5-beta"
#r "nuget: FsMath, 0.0.2"
#r "nuget: FSharpAux, 2.1.0"
#r "nuget: ProteomIQon, 0.0.10"
#r "nuget: BioFSharp.Mz, 0.1.5-beta"
#r "nuget: FSharp.Stats, 0.4.0"
#r "nuget: BioFSharp, 2.0.0-beta4"
#r "nuget: BioFSharp.Mz, 0.1.5-beta"
#r "nuget: FSharpAux, 1.0.0"
#r "nuget: FSharpAux.IO, 1.0.0"
#r "nuget: ProteomIQon"

open System.IO
open System.Data.SQLite
open ProteomIQon
open ProteomIQon.Core
open ProteomIQon.Core.MzIO
open ProteomIQon.Dto
open FSharp.Stats
open BioFSharp.Mz.Quantification
open BioFSharp.Mz
open FSharpAux.IO.SchemaReader
open Plotly.NET
open BioFSharp
open MzIO.Processing
open BioFSharp.Mz.SearchDB

open Plotly.NET
open System
open FSharpAux 
open BioFSharp.ImgP
open BioFSharp.Mz.Quantification
open BioFSharp.Mz
open System
open FsMath
open FSharpAux
open FSharp.Stats
open MzIO
open MzIO.Processing
open MzLite
open MzIO.MzSQL
open System.Data
// open Plotly.NET
// open Plotly.NET.TraceObjects
// open Plotly.NET.LayoutObjects
// open Plotly.NET.Interactive
open MzIO.Commons.Arrays
open MzIO.Processing.MzIOLinq
open MzIO.Binary
open Deedle
open System.Data.SQLite
open System.Linq
open MzIO.IO.MzML
open FSharp.Stats.Signal
open System.IO

do fsi.AddPrinter(fun (printer:Deedle.Internal.IFsiFormattable) -> "\n" + (printer.Format()))



let inReader = new MzSQL(@"/home/paulinehans/Dokumente/testRunGabor/binned_spectra_1.000.mzlite")
inReader.Connection.Open()
let inRunID  = "sample=0"
let inTr = inReader.BeginTransaction()

// let ms = inReader.ReadMassSpectrum("merged=475568 frame=53975 scanStart=1 scanEnd=918")
// ms |> MzIO.Processing.MassSpectrum.getScanTime

let initRTProfile (readspecPeaks:string -> Peak1DArray)  (rtIndex: IMzIOArray<RtIndexEntry>) (rtRange: RangeQuery) (mzRange: RangeQuery) =
    let entries = RtIndexEntry.Search(rtIndex, rtRange).ToArray()
    //printfn "RtProfile %i" entries.Length
    [|
        for rtIdx = 0 to entries.Length-1 do
            let entry = entries.[rtIdx]
            let peaks = (readspecPeaks entry.SpectrumID).Peaks
            let p = 
                (RtIndexEntry.MzSearch (peaks, mzRange)).DefaultIfEmpty(Peak1D(0., mzRange.LockValue))
                |> Seq.map (fun x -> RtIndexEntry.AsPeak2D (x, entry.Rt))
                |> Seq.toArray
                
            p      
    |]


let filePathQuant = ("/home/paulinehans/Dokumente/testRunGabor/GaborOutput/binned_spectra_1.000_quant")
let quantFile = Frame.ReadCsv (filePathQuant, hasHeaders =true, separators = "\t")
let quantifiedRange = 
    quantFile
    |> Frame.mapRows (fun rk s -> 
        {|
            MassLight = s.GetAs<float>("QuantMz_Light");
            RTTraceLight = s.GetAs<string>("RtTrace_Light").Split";" |> Array.map float |> Array.median;
            QuantLight = s.GetAs<float>("Quant_Light")
            ParamsLight = s.GetAs<string>("Params_Light")
        |}
    )
    |> Series.values
    |> Seq.toArray
    |> Array.sortByDescending (fun r -> r.QuantLight)
    |> Array.take 1
let retTimeIdxed = Query.getMS1RTIdx inReader inRunID
let allQuerys =
    quantifiedRange 
    |> Array.map (fun r -> 
        let rt = RangeQuery (r.RTTraceLight, 2.0)
        let mz = RangeQuery (r.MassLight, 0.1)
        initRTProfile inReader.ReadSpectrumPeaks retTimeIdxed  rt mz
    )  

let transform  = 
    let waveletData = 
        allQuerys 
        |> Array.map (fun x -> 
            x    
            |> Array.collect (fun x -> 
                x 
                |> Array.filter (fun y -> y.IonMobility.IsSome)
            ) 
            |> Array.map (fun x -> x.Rt, x.IonMobility.Value, x.Intensity)  
        )
    waveletData
    |> Array.concat
transform



//let inReaderMS = new MzMLReaderMIRIM("D:\Testfiles\MS292NT_Slot1-25_1_1375.mzML")
//let inReaderPeaks = new MzMLReaderMIRIM("D:\Testfiles\MS292NT_Slot1-25_1_1375.mzML")
//let inTrMS = inReaderMS.BeginTransaction()
//let inTrPeaks = inReaderPeaks.BeginTransaction()

//let ms2 = inReaderMS.ReadMassSpectrum("merged=377907 frame=45978 scanStart=1 scanEnd=918")

//let unzipIMzliteArray (a:IMzIOArray<Peak1D>) = 
//    let mzData = Array.zeroCreate a.Length
//    let intensityData = Array.zeroCreate a.Length
//    for i = 0 to a.Length-1 do 
//        let peak = a.[i]
//        mzData.[i] <- peak.Mz
//        intensityData.[i] <- peak.Intensity
//    mzData,intensityData

//let old =
//    inReaderPeaks.getSpecificPeak1DArraySequentialWithMIRIM(ms2.ID)
//    |> fun pa ->
//        let mzData,intensityData = unzipIMzliteArray pa.Peaks
//        let ionMobilityData = (pa?Mirim |> unbox<float array>)
//        mzData,intensityData,ionMobilityData

let new1 =
    inReader.SelectPeak1DArray(ms.ID, true)
    |> fun pa ->
        PeakArray.mzIntensityIonMobilityArrayOf pa

//Chart.Point3D(old |> fun (x,y,z) -> Array.zip3 x y z)
//|> Chart.withXAxisStyle("Mz")
//|> Chart.withYAxisStyle("Intensity")
//|> Chart.withZAxisStyle("Ion Mobility")
//|> Chart.show
//Chart.Point3D(new1 |> fun (x,y,z) -> Array.zip3 x y z)
//|> Chart.withXAxisStyle("Mz")
//|> Chart.withYAxisStyle("Intensity")
//|> Chart.withZAxisStyle("Ion Mobility")
//|> Chart.show

type SeqWrapper<'T>(sequence:'T seq) =
            
    interface System.Collections.Generic.IEnumerable<'T> with
        member this.GetEnumerator() =
            sequence.AsEnumerable().GetEnumerator()
 
    interface System.Collections.IEnumerable with
        member this.GetEnumerator() = 
            sequence.GetEnumerator()

    interface IMzIOArray<'T> with
        member this.Length = 
            sequence.Count()

        member this.Item with get idx = sequence |>Seq.item idx

new1
|> fun (x,y,z) -> Array.zip3 x y z
|> Array.filter (fun (x,y,z) -> z < 0.9607 + 0.05 && z > 0.9607 - 0.05)
|> Array.length
|> Array.filter (fun (x,y,z) -> x < 631.3194 + 2. && x > 631.3194 - 2.)

let mzRange2 = RangeQuery(631.3194, 2.)
let imRange2 = RangeQuery(0.9607, 0.05)
let peaks = 
    (inReader.ReadSpectrumPeaks "merged=475568 frame=53975 scanStart=1 scanEnd=918").Peaks
    |> Seq.filter (fun x -> x.Mz > 1090 && x.Mz < 1095)
    //|> Seq.sortBy (fun x -> x.Mz)
    //|> SeqWrapper<Peak1D>

Chart.Point(peaks |> Seq.map (fun x -> x.IonMobility.Value,x.Intensity))
|> Chart.show

#time
let t = 
    (RtIndexEntry.IonMobilitySearch (peaks, imRange2)).DefaultIfEmpty(new Peak1D(0., mzRange2.LockValue, imRange2.LockValue))
    |> Seq.filter (fun x -> x.Mz < mzRange2.HighValue && x.Mz > mzRange2.LowValue)
t
|> Seq.toArray
|> Array.groupBy (fun x -> x.Mz)
|> Array.map (fun (mz,peaks) -> 
    let i = peaks |> Array.sumBy (fun x -> x.Intensity)
    new Peak1D(i, mz)
)

//old = new1
let retTimeIdxed = Query.getMS1RTIdx inReader inRunID

/// Extract a rt profile for specified target mass and rt range.
/// Mz range peak aggregation is closest lock mz.
/// Profile array with index corresponding to continous mass spectra over rt range and mz range given.
let initRTProfile3 (rtIndex: IMzIOArray<RtIndexEntry>) (rtRange: RangeQuery) (mzRange: RangeQuery) (ionMobilityRange: RangeQuery)=
    let entries = RtIndexEntry.Search(rtIndex, rtRange)
    entries
    |>Seq.map (fun entry ->
        let peaks = (inReader.ReadSpectrumPeaks entry).Peaks
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
    //profile


let rtRange = RangeQuery(95.6818666666666, 10.)
let mzRange = RangeQuery(1093.0861, 2.)
let imRange = RangeQuery(1.29, 0.05)

//let woIM =
//    inReader.RtProfile(retTimeIdxed, rtRange, mzRange)

let wIM =
    initRTProfile3 retTimeIdxed rtRange mzRange imRange
wIM |> Array.filter (fun x -> x.[0].Intensity > 0)
wIM.[0]

//woIM.Length
//wIM.Length

//woIM |> Array.map (fun x -> x.Rt,x.Intensity)
//|> Chart.Point
//|> Chart.show
1
wIM |> Array.map (fun x -> x.Rt,x.Intensity)
|> Chart.Point
|> Chart.show





let cn = SearchDB.getDBConnection @"C:\Users\jonat\OneDrive\Desktop\ChlamyTest.db"
let memoryDB = SearchDB.copyDBIntoMemory cn 
let pepDBTr = memoryDB.BeginTransaction()


open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute
type PSMStatisticsResultFragpipe = {
    // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space separated by '_'
    [<FieldAttribute(0)>]
    PSMId                        : string
    [<FieldAttribute(1)>]
    PepSequenceID                : int
    [<FieldAttribute(2)>]
    ModSequenceID                : int
    [<FieldAttribute(3)>]
    ScanTime                     : float
    [<FieldAttribute(4)>]
    Charge                       : int
    [<FieldAttribute(5)>]
    PrecursorMZ                  : float
    [<FieldAttribute(6)>]
    TheoMass                     : float
    [<FieldAttribute(7)>]
    AbsDeltaMass                 : float
    [<FieldAttribute(8)>]
    IonMobility                  : float
    [<FieldAttribute(9)>]
    PeptideLength                : int
    [<FieldAttribute(10)>]
    MissCleavages                : int
    [<FieldAttribute(11)>]
    Expectscore                  : float
    [<FieldAttribute(12)>]
    Hyperscore                   : float
    [<FieldAttribute(13)>]
    Probability                  : float
    [<FieldAttribute(14)>]
    StringSequence               : string
    [<FieldAttribute(15)>]
    ProteinNames                 : string
    [<FieldAttribute(16)>]
    GlobalMod                    : int
}

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
let substractBaseLine (baseLineParams:Domain.BaseLineCorrection) (yData:float []) =
    if yData.Length > 500 then
        yData
    else
        let baseLine = FSharp.Stats.Signal.Baseline.baselineAls' baseLineParams.MaxIterations baseLineParams.Lambda baseLineParams.P yData |> Array.ofSeq
        Array.map2 (fun y b ->
                        let c = y - b
                        if c < 0. then 0. else c
                    ) yData baseLine

///
let initGetProcessedXIC (baseLineCorrection:Domain.BaseLineCorrection option) getPeaks idx scanTimeWindow mzWindow_Da meanScanTime meanPrecMz =
    let rtQuery = Query.createRangeQuery meanScanTime scanTimeWindow
    let mzQuery = Query.createRangeQuery meanPrecMz mzWindow_Da
    let retData',itzData' =
        let tmp =
            getPeaks idx rtQuery mzQuery
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
        retData', substractBaseLine baseLineParams itzData', itzData'
    | None ->
        retData',itzData',itzData'

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
let weightedMean (weights:seq<'T>) (items:seq<'T>) =
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
    
let getInferredXic getXic targetScanTime targetMz =
    let (retData,itzData,uncorrectedItzData)   =
            getXic targetScanTime targetMz 
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

let dBParams     = SearchDB'.getSDBParams memoryDB
//let massLookUp = prepareSelectMassByModSequenceAndGlobalMod memoryDB
let peptideLookUp = SearchDB'.getThreadSafePeptideLookUpFromFileBySequenceAndGMod memoryDB dBParams
let calcIonSeries aal = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal


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

let qpsms =
    Csv.CsvReader<PSMStatisticsResultFragpipe>(SchemaMode=Csv.Fill).ReadFile(@"D:\Testfiles\binned_spectra_1.000.mzlite.qpsm",'\t',false,1)
    |> Array.ofSeq

let massSpectra = inReader.ReadMassSpectra(inRunID)
/// Returns a Array that contains all MS1s sorted by their scanTime
let ms1SortedByScanTime =
    massSpectra
    |> Seq.filter (fun ms -> MassSpectrum.getMsLevel ms = 1)
    |> Seq.map (fun ms -> MassSpectrum.getScanTime ms, ms)
    |> Seq.sortBy fst
    |> Array.ofSeq

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
    stDev, scanTimeToMzCorrection
1
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
    refined

let comparePredictedAndMeasuredIsotopicCluster = initComparePredictedAndMeasuredIsotopicCluster inReader ms1SortedByScanTime ms1AccuracyEstimate

let mzWindow = 
    let mzW = ms1AccuracyEstimate*4.
    mzW

let initGetProcessedXIC_IM scanTimeWindow mzWindow_Da imWindow meanScanTime meanPrecMz im =
    let rtQuery = Query.createRangeQuery meanScanTime scanTimeWindow
    let mzQuery = Query.createRangeQuery meanPrecMz mzWindow_Da
    let imQuery = Query.createRangeQuery im imWindow
    let retData',itzData' =
        let tmp =
            inReader.RtProfile (retTimeIdxed, rtQuery, mzQuery, imQuery)
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

    retData',itzData',itzData'

let waveletParams :Domain.WaveletParameters = 
    {
        Borderpadding           = None    
        BorderPadMethod         = FSharp.Stats.Signal.Padding.BorderPaddingMethod.Random 
        InternalPaddingMethod   = FSharp.Stats.Signal.Padding.InternalPaddingMethod.LinearInterpolation 
        HugeGapPaddingMethod    = FSharp.Stats.Signal.Padding.HugeGapPaddingMethod.Zero
        HugeGapPaddingDistance  = 100.
        MinPeakDistance         = None
        MinPeakLength           = Some 0.1
        MaxPeakLength           = 1.5 
        NoiseQuantile           = 0.01 
        MinSNR                  = 0.01  
    }

let XicExtraction: Domain.XicExtraction = 
    {
        TopKPSMs                     = None
        ScanTimeWindow               = 2.
        MzWindow_Da                  = Domain.Window.Estimate
        XicProcessing                = ProteomIQon.Domain.XicProcessing.Wavelet waveletParams
    }


let grouped =
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

let getXIC = initGetProcessedXIC_IM 2. 0.01 0.05
1
let topXICs =
    grouped
    |> Array.sortBy (fun (_,psms) -> psms |> Array.minBy (fun x -> x.Expectscore))
    |> Array.take 10
    |> Array.mapi (fun i (pepIon, psms) ->
        printfn "%A" i
        let bestQValue,prots = psms |> Array.minBy (fun x -> x.Expectscore) |> fun x -> x.Expectscore,x.ProteinNames
        let unlabledPeptide = peptideLookUp pepIon.Sequence 0
        let psmsWithMatchedSums = countMatchedMasses unlabledPeptide psms
        let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
        let theoMz = Mass.toMZ unlabledPeptide.Mass (float pepIon.Charge)
        let averagePSM = average getXIC id theoMz psmsWithMatchedSums
        averagePSM
    )

let getXIC2 scanTimeWindow mzWindow_Da imWindow meanScanTime meanPrecMz im=
    let rtQuery = Query.createRangeQuery meanScanTime scanTimeWindow
    let mzQuery = Query.createRangeQuery meanPrecMz mzWindow_Da
    let imQuery = Query.createRangeQuery im imWindow
    let retData',itzData' =
        let tmp =
            inReader.RtProfile (retTimeIdxed, rtQuery, mzQuery, imQuery)
            |> Array.map (fun (p:MzIO.Binary.Peak2D) -> p.Rt , p.Intensity)
        tmp
        |> Array.unzip
    retData',itzData'

let topXICsDirect =
    qpsms
    |> Array.sortByDescending (fun x -> x.Hyperscore)
    |> Array.mapi (fun i x ->
        printfn "%A" i
        let rt,intensity = getXIC2 2. 0.01 0.1 x.ScanTime x.PrecursorMZ x.IonMobility
        if intensity |> Array.sum > 0. then
            Chart.Point (rt, intensity)
            |> Chart.saveHtml $@"D:\Testfiles\Charts\{x.PSMId}_{x.StringSequence}_{x.Charge}_{intensity |> Array.max}"
    )

topXICsDirect
|> Array.map (fun (x,y,_) ->
    if y |> Array.sum > 0. then
        Chart.Point (x, y)
        |> Some
    else
        None
)
|> Array.choose id
|> fun x -> 
    printfn "%A"x.Length
    x
|> Array.map Chart.show