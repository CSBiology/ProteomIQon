#r "nuget: ProteomIQon, 0.0.8"
#r "nuget: MzIO, 0.1.1"
#r "nuget: MzIO.Processing, 0.1.2"
#r "nuget: MzIO.SQL, 0.1.4"
#r "nuget: MzIO.MzML, 0.1.6"
#r "nuget: System.Data.SQLite.Core, 1.0.118"

open ProteomIQon
open ProteomIQon.Domain
open ProteomIQon.Core
open System.IO
open System.Collections.Generic
open BioFSharp.Mz
open MzIO
open MzIO.Binary
open MzIO.Model
open MzIO.Processing
open MzIO.IO
open MzIO.MzSQL
open MzIO.IO.MzML
open ProteomIQon.Core.MzIO.Processing
open ProteomIQon.Dto
open ProteomIQon.Domain


    ///
let private initPeakPicking (*(logger: NLog.Logger)*) (peakPickingParams:PeakPicking) (peaks: Peak1DArray) =

    match peakPickingParams with
    | PeakPicking.ProfilePeaks ->
        peaks.Peaks
        |> Core.MzIO.Peaks.unzipIMzliteArray
    | PeakPicking.Centroid (CentroidizationMode.Wavelet waveletParams) ->
        match waveletParams.PaddingParams with
        | Some pParams ->
            let initPaddingParameters yThreshold =
                SignalDetection.Padding.createPaddingParameters
                    yThreshold
                    pParams.MaximumPaddingPoints
                    pParams.Padding_MzTolerance
                    pParams.WindowSize
                    pParams.SpacingPerc

            let initwaveletParameters yThreshold =
                SignalDetection.Wavelet.createWaveletParameters
                    waveletParams.NumberOfScales
                    yThreshold
                    waveletParams.Centroid_MzTolerance
                    waveletParams.SNRS_Percentile
                    waveletParams.MinSNR
                    waveletParams.RefineMZ
                    waveletParams.SumIntensities
            match waveletParams.YThreshold with
            | YThreshold.Fixed yThreshold ->
                let paddingParams = initPaddingParameters yThreshold
                let waveletParameters = initwaveletParameters yThreshold
                let mzData, intensityData =
                    peaks.Peaks
                    |> Core.MzIO.Peaks.unzipIMzliteArray
                let paddedMz,paddedIntensity =
                    SignalDetection.Padding.paddDataBy paddingParams mzData intensityData
                BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters paddedMz paddedIntensity
            | YThreshold.MinSpectrumIntensity ->
                //logger.Trace (sprintf "ms2Centroidization with: %A" waveletParams)
                let mzData, intensityData =
                    peaks.Peaks
                    |> Core.MzIO.Peaks.unzipIMzliteArray
                let yThreshold = Array.min intensityData
                let paddedMz,paddedIntensity =
                    let paddingParams = initPaddingParameters yThreshold
                    SignalDetection.Padding.paddDataBy paddingParams mzData intensityData
                let waveletParameters = initwaveletParameters yThreshold
                BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters paddedMz paddedIntensity
        | None ->
            let initwaveletParameters yThreshold =
                SignalDetection.Wavelet.createWaveletParameters
                    waveletParams.NumberOfScales
                    yThreshold
                    waveletParams.Centroid_MzTolerance
                    waveletParams.SNRS_Percentile
                    waveletParams.MinSNR
                    waveletParams.RefineMZ
                    waveletParams.SumIntensities
            match waveletParams.YThreshold with
            | YThreshold.Fixed yThreshold ->
                //logger.Trace (sprintf "ms1Centroidization with: %A" waveletParams)
                let waveletParameters = initwaveletParameters yThreshold
                let mzData, intensityData =
                    peaks.Peaks
                    |> Core.MzIO.Peaks.unzipIMzliteArray
                BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters mzData intensityData
            | YThreshold.MinSpectrumIntensity ->
                let mzData, intensityData =
                    peaks.Peaks
                    |> Core.MzIO.Peaks.unzipIMzliteArray
                let yThreshold = Array.min intensityData
                let waveletParameters = initwaveletParameters yThreshold
                BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters mzData intensityData
    | PeakPicking.Centroid CentroidizationMode.Manufacturer ->
        failwith "Manufacturer peak picking is only supported for .baf (Bruker) files."
    
let binBy (projection: 'a -> float) bandwidth (data: seq<'a>) =
    if bandwidth = 0. then raise (System.DivideByZeroException("Bandwidth cannot be 0."))
    let halfBw = bandwidth / 2.0
    let decBandwidth = decimal bandwidth
    let tmp = 
        data
        |> Seq.groupBy (fun x -> (decimal (projection x) / decBandwidth) |> float |> floor) 
        |> Seq.map (fun (k,values) -> 
            let count = (Seq.length(values)) |> float
            if k < 0. then
                ((k  * bandwidth) + halfBw, values)   
            else
                ((k + 1.) * bandwidth) - halfBw, values)
        |> Seq.sortBy fst
    tmp

let fixSpectrum (m:MzIO.Model.MassSpectrum) =
    if isNull(m.Precursors) then
        m.Precursors <- new MzIO.Model.PrecursorList()
    if isNull(m.Scans) then
        m.Scans <- new MzIO.Model.ScanList()
    if isNull(m.Products) then
        m.Products <- new MzIO.Model.ProductList()
    m

let createPeak1DArrayCopy (source: MzIO.Binary.Peak1DArray) =
    let pa = MzIO.Binary.Peak1DArray()
    pa.CompressionType <- source.CompressionType
    pa.IntensityDataType <- source.IntensityDataType
    pa.MzDataType <- source.MzDataType
    pa

let createBinnedPeaks (copyMirim: bool) (binSize: float) (peakArray: MzIO.Binary.Peak1DArray) = 

    let zippedPeaks = peakArray.Peaks |> Seq.zip (peakArray?Mirim |> unbox<float array>)

    let binnedPeakData =
        zippedPeaks
        |> binBy (fun (mirim, peak) -> mirim) binSize

    binnedPeakData
    |> Seq.map(fun (bin, binnedData) ->
        let pa = createPeak1DArrayCopy peakArray
        pa.Peaks <-
            MzIO.Commons.Arrays.ArrayWrapper(
                binnedData
                |> Seq.map snd
                |> Seq.toArray
            )
        if copyMirim then
            pa?Mirim <-
                binnedData
                |> Seq.map fst
                |> Array.ofSeq
        bin,pa
    )

let insertSpectrum (compress:BinaryDataCompressionType) (outReader: MzSQL.MzSQL) (runID:string)
    (ms1PeakPicking: Peak1DArray -> float [] * float []) (ms2PeakPicking: Peak1DArray -> float [] * float [])
        (spectrum: MassSpectrum) (peaks: Peak1DArray) =
    match MassSpectrum.getMsLevel spectrum with
    | 1 ->
        let mzData,intensityData =
            try
            ms1PeakPicking peaks
            with
            | _ -> [||],[||]
        if Array.isEmpty mzData || Array.isEmpty intensityData then ()
        else
            let peaks' = PeakArray.createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData 
            outReader.Insert (runID, spectrum, peaks')
    | 2 ->
        let mzData,intensityData =
            try
            ms2PeakPicking peaks
            with
            | _ -> [||],[||]
        if Array.isEmpty mzData || Array.isEmpty intensityData then ()
        else
            let peaks' = PeakArray.createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData
            outReader.Insert (runID, spectrum, peaks')
    | _ ->
        failwith "Only mass spectra of level 1 and 2 are supported."

/// Returns the default runID used by manufacturers
let getDefaultRunID (mzReader:IMzIODataReader) = 
    match mzReader with
    | :? MzSQL as r                 -> "sample=0"
    | :? MzML.MzMLReader as r       -> "sample=0"
    | :? MzMLReaderMIRIM as r       -> "sample=0"

let processFile (processParams:MzMLtoMzLiteParams) (outputDir:string) (instrumentOutput:string) =

    //let logger = Logging.createLogger (Path.GetFileNameWithoutExtension instrumentOutput)

    //logger.Trace (sprintf "Input file: %s" instrumentOutput)
    //logger.Trace (sprintf "Output directory: %s" outputDir)
    //logger.Trace (sprintf "Parameters: %A" processParams)

    
    //let tmp = File.ReadAllText instrumentOutput
    //File.WriteAllText(instrumentOutput, tmp.Replace("&quot;", ""))

    let inReaderMS = new MzMLReaderMIRIM(instrumentOutput)
    let inReaderPeaks = new MzMLReaderMIRIM(instrumentOutput)
    let inRunID  = getDefaultRunID inReaderMS
    let inTrMS = inReaderMS.BeginTransaction()
    let inTrPeaks = inReaderPeaks.BeginTransaction()

    let ms1PeakPicking = initPeakPicking (*logger*) processParams.MS1PeakPicking
    let ms2PeakPicking = initPeakPicking (*logger*) processParams.MS2PeakPicking

    let outDirPath =
        let fileName = Path.GetFileNameWithoutExtension instrumentOutput
        Path.Combine(outputDir, fileName)
            
    //logger.Trace $"Creating directory {outDirPath} for binned results of {instrumentOutput}"
        
    Directory.CreateDirectory outDirPath |> ignore

    //logger.Trace $"Reading spectra from {instrumentOutput}"

    let spectra = inReaderMS.ReadMassSpectra(inRunID)
    //inReaderMS.ResetReader()

    //logger.Trace "Done reading spectra"
    //logger.Trace $"Reading model from {instrumentOutput}"
        
    //let model = inReaderMS.Model
    //inReaderMS.ResetReader()
        
    //logger.Trace "Done reading model"

    //logger.Trace $"Start writing binned mzlite files"
    //logger.Trace $"Total number of binned files: {spectrumMap.Count}"
    let connectionMap = new Dictionary<string, MzSQL.MzSQL*System.Data.SQLite.SQLiteTransaction>()
    spectra
    |> Seq.take 10000
    |> Seq.iteri (fun i spectrum ->
        if i % 1000 = 0 then printfn "%i" i
        let data = inReaderPeaks.getSpecificPeak1DArraySequentialWithMIRIM(spectrum.ID)
        let binResult = createBinnedPeaks false 0.002 data
        binResult
        |> Seq.iter(fun (bin, peaks) ->
            let outFile = Path.Combine(outputDir, $"binned_spectra_%.3f{bin}.mzlite")
            let outReader,outTr =
                if connectionMap.ContainsKey(outFile) then
                    connectionMap.[outFile]
                else
                    let outReader = new MzSQL(outFile)
                    let _ = outReader.Open()
                    let outTr = outReader.BeginTransaction()
                    //try
                    //    outReader.InsertModel model
                    //    //logger.Trace "Model inserted."
                    //with
                    //    | ex -> failwith $"Inserting model failed: {ex}"
                    connectionMap.Add(outFile, (outReader, outTr))
                    outReader, outTr
            let outRunID  = getDefaultRunID outReader
            insertSpectrum processParams.Compress outReader outRunID ms1PeakPicking ms2PeakPicking (changeScanTimeToMinutes (fixSpectrum spectrum)) peaks
        )
    )
    for x in connectionMap do
        let outReader = fst x.Value
        let outTr = snd x.Value
        outTr.Commit()
        outTr.Dispose()
        outReader.Dispose()
    inTrPeaks.Commit()
    inTrPeaks.Dispose()
    inTrMS.Commit()
    inTrMS.Dispose()
    inReaderMS.Dispose()
    inReaderPeaks.Dispose()
    //logger.Trace "Done."

let deserialized = 
    System.IO.File.ReadAllText(@"C:\Users\jonat\source\repos\ProteomIQon\src\ProteomIQon\defaultParams\TIMsMzMLtoMzLiteParams.json")
    |> Json.deserialize<Dto.MzMLtoMzLiteParams>
    |> PreprocessingParams.toDomain
#time
processFile deserialized @"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\outTest" "C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\Ara_60min_wTrap_Aurora_DDA_Slot1-4_113.mzML"

let reader = new MzMLReaderMIRIM("C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\Ara_60min_wTrap_Aurora_DDA_Slot1-4_113.mzML")
//let reader2 = new MzMLReaderMIRIM("C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\Ara_60min_wTrap_Aurora_DDA_Slot1-4_113.mzML")

let inTrMS = reader.BeginTransaction()
reader.Model
//let inTrMS2 = reader.BeginTransaction()
let spectra = reader.ReadMassSpectra("sample=0")
spectra |> Array.ofSeq
//reader.ResetReader()
//spectra
//|> Seq.mapi (fun i s -> 
//    if i % 1000 = 0 then printfn "%i" i
//    reader2.getSpecificPeak1DArraySequentialWithMIRIM(s.ID)
//)
//|> Array.ofSeq

//spectra.[0].ID
//reader.ResetReader()
//let data = reader.getSpecificPeak1DArraySequentialWithMIRIM(spectra.[50].ID)
//data.Peaks.Length
//(data?Mirim :?> float array).Length

//let binResult = createBinnedPeaks false 0.002 data

//let mz,intensity = data.Peaks |> Core.MzIO.Peaks.unzipIMzliteArray
//let ionMobility = (data?Mirim :?> float array)
//let a = Array.zip3 mz intensity ionMobility
//open Plotly.NET
//1
//Chart.Scatter3d(a, StyleParam.Mode.Markers)
//|> Chart.withX_AxisStyle "m/z"
//|> Chart.withY_AxisStyle "Intensity"
//|> Chart.withZ_AxisStyle "Voltage"
//|> Chart.withSize (1200., 900.)
//|> Chart.Show

//Chart.Point(mz,intensity)
//|> Chart.withX_AxisStyle "m/z"
//|> Chart.withY_AxisStyle "Intensity"
//|> Chart.withSize (1200., 900.)
//|> Chart.Show