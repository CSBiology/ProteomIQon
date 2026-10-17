namespace ProteomIQon

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


module MzMLIonMobilityToMzLite =
    ///
    let private initPeakPicking (logger: NLog.Logger) (peakPickingParams:PeakPicking) (peaks: Peak1DArray) =

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

    let toIMPeaks (peakArray: MzIO.Binary.Peak1DArray) = 
        let imArray = (peakArray?Mirim |> unbox<float array>)
        let peaks = peakArray.Peaks
        let zipped = Seq.zip imArray peaks
        let pa = createPeak1DArrayCopy peakArray
        pa.Peaks <-
            MzIO.Commons.Arrays.ArrayWrapper(
                zipped
                |> Seq.map(fun (ionMobility, peak) ->
                    Peak1D(peak.Intensity,peak.Mz,ionMobility)
                )
                |> Seq.toArray
            )
        pa
        
    let unzipIonMobilityMzliteArray (a:Commons.Arrays.IMzIOArray<Peak1D>) = 
        let mzData = Array.zeroCreate a.Length
        let intensityData = Array.zeroCreate a.Length
        let ionMobilityData = Array.zeroCreate a.Length
        for i = 0 to a.Length-1 do 
            let peak = a.[i]
            mzData.[i] <- peak.Mz
            intensityData.[i] <- peak.Intensity
            ionMobilityData.[i] <- peak.IonMobility.Value
        mzData,intensityData,ionMobilityData

    /// Creates Peak1DArray of mzData array and intensityData Array
    let createPeak1DArray compression mzBinaryDataType intensityBinaryDataType ionMobilityDataType (mzData:float []) (intensityData:float []) (ionMobilityData:float []) =
        let zipedData = Array.map3 (fun mz intz im -> Peak1D(intz,mz,im)) mzData intensityData ionMobilityData
        let newPeakA = Commons.Arrays.MzIOArray.ToMzIOArray zipedData
        let peak1DArray = new Peak1DArray(compression,intensityBinaryDataType, mzBinaryDataType,newPeakA, ionMobilityDataType)
        peak1DArray

    let insertSpectrum (compress:BinaryDataCompressionType) (outReader: MzSQL.MzSQL) (runID:string)
        (ms1PeakPicking: Peak1DArray -> float [] * float []) (ms2PeakPicking: Peak1DArray -> float [] * float [])
            (spectrum: MassSpectrum) (peaks: Peak1DArray) =
        match MassSpectrum.getMsLevel spectrum with
        | 1 ->
            //let mzData,intensityData =
            //    try
            //    ms1PeakPicking peaks
            //    with
            //    | _ -> [||],[||]
            //if Array.isEmpty mzData || Array.isEmpty intensityData then ()

            let mzData, intensityData, ionMobilityData =
                unzipIonMobilityMzliteArray peaks.Peaks

            //else
            let peaks' = createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData ionMobilityData
            outReader.Insert (runID, spectrum, peaks')
        | 2 ->
            //let mzData,intensityData =
            //    try
            //    ms2PeakPicking peaks
            //    with
            //    | _ -> [||],[||]
            //if Array.isEmpty mzData || Array.isEmpty intensityData then ()
            //else
            let mzData, intensityData, ionMobilityData =
                unzipIonMobilityMzliteArray peaks.Peaks

            let peaks' = createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData ionMobilityData
            outReader.Insert (runID, spectrum, peaks')
        | _ ->
            failwith "Only mass spectra of level 1 and 2 are supported."

    /// Returns the default runID used by manufacturers
    let getDefaultRunID (mzReader:IMzIODataReader) = 
        match mzReader with
        | :? MzSQL as r                 -> "sample=0"
        | :? MzML.MzMLReader as r       -> "sample=0"
        | :? MzMLReaderMIRIM as r       -> "sample=0"

    let processFile (processParams:MzMLtoMzLiteParams) (fixFile: bool) (outputDir:string) (instrumentOutput:string) =

        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension instrumentOutput)

        logger.Trace (sprintf "Input file: %s" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)

        if fixFile then
            logger.Trace "Fixing file."
            File.ReadAllLines instrumentOutput
            |> Array.map (fun s -> s.Replace ("&quot",""))
            |> fun c -> File.WriteAllLines (instrumentOutput,c)

        //let tmp = File.ReadAllText instrumentOutput
        //File.WriteAllText(instrumentOutput, tmp.Replace("&quot;", ""))

        let inReaderMS = new MzMLReaderMIRIM(instrumentOutput)
        let inReaderMS = new MzMLReaderMIRIM(instrumentOutput)
        let inRunID  = getDefaultRunID inReaderMS
        let inTrMS = inReaderMS.BeginTransaction()
        // the peaks of the spectra of the usual layout are read directly from the file, the others with the MIRIM reader
        let rawPeakReader = new MzMLRawPeakReaderMIRIM(instrumentOutput)
        let fallbackPeakReader = lazy (new MzMLReaderMIRIM(instrumentOutput))
        let useRawPeakReader = processParams.Compress = BinaryDataCompressionType.ZLib
        let peakDataEncoder = TimsPeaks.PeakDataEncoder()

        let ms1PeakPicking = initPeakPicking logger processParams.MS1PeakPicking
        let ms2PeakPicking = initPeakPicking logger processParams.MS2PeakPicking

        let outFilePath =
            let fileName = Path.GetFileNameWithoutExtension instrumentOutput
            Path.Combine(outputDir, (fileName + ".mzlite"))
            
        logger.Trace $"Creating directory {outFilePath} for binned results of {instrumentOutput}"
        
        Directory.CreateDirectory outputDir |> ignore

        logger.Trace $"Reading model from {instrumentOutput}"

        let model = inReaderMS.Model
        inReaderMS.ResetReader()

        logger.Trace "Done reading model"

        logger.Trace $"Reading spectra from {instrumentOutput}"

        let spectra = inReaderMS.ReadMassSpectra(inRunID)

        logger.Trace "Done reading spectra"

        logger.Trace $"Start writing binned mzlite files"
        let outReader = new MzSQL(outFilePath)
        let _ = outReader.Open()
        let outTr = outReader.BeginTransaction()
        let outRunID  = getDefaultRunID outReader
        // the description of the peak arrays is the same for all spectra
        let peakArrayJson =
            lazy (
                createPeak1DArray processParams.Compress BinaryDataType.Float64 BinaryDataType.Float64 BinaryDataType.Float64 [||] [||] [||]
                |> MzIO.Json.MzIOJson.ToJson
            )
        let insertCommand =
            new System.Data.SQLite.SQLiteCommand(
                "INSERT INTO Spectrum (RunID, SpectrumID, Description, PeakArray, PeakData) VALUES(@runID, @spectrumID, @description, @peakArray, @peakData)",
                outReader.Connection)
        for (name, dbType) in [ "@runID", System.Data.DbType.String; "@spectrumID", System.Data.DbType.String; "@description", System.Data.DbType.String; "@peakArray", System.Data.DbType.String; "@peakData", System.Data.DbType.Binary ] do
            insertCommand.Parameters.Add(name, dbType) |> ignore
        spectra
        |> Seq.iteri (fun i spectrum ->
            if i % 1000 = 0 then logger.Trace $"{i}"
            let spectrum' = changeScanTimeToMinutes (fixSpectrum spectrum)
            if useRawPeakReader && rawPeakReader.TryReadPeaks spectrum.ID then
                match MassSpectrum.getMsLevel spectrum' with
                | 1 | 2 -> ()
                | _ -> failwith "Only mass spectra of level 1 and 2 are supported."
                insertCommand.Parameters.["@runID"].Value <- outRunID
                insertCommand.Parameters.["@spectrumID"].Value <- spectrum'.ID
                insertCommand.Parameters.["@description"].Value <- MzIO.Json.MzIOJson.MassSpectrumToJson spectrum'
                insertCommand.Parameters.["@peakArray"].Value <- peakArrayJson.Value
                insertCommand.Parameters.["@peakData"].Value <- peakDataEncoder.Encode(rawPeakReader.Intensity, rawPeakReader.Mz, rawPeakReader.IonMobility, rawPeakReader.Count)
                insertCommand.ExecuteNonQuery() |> ignore
            else
                let data = fallbackPeakReader.Value.getSpecificPeak1DArraySequentialWithMIRIM(spectrum.ID)
                let dataIM = toIMPeaks data
                insertSpectrum processParams.Compress outReader outRunID ms1PeakPicking ms2PeakPicking spectrum' dataIM
        )
        try
            outReader.InsertModel model
            logger.Trace $"Model inserted for {outFilePath}."
        with
            | ex -> failwith $"Inserting model failed: {ex}"
        outTr.Commit()
        outTr.Dispose()
        outReader.Dispose()
        inTrMS.Commit()
        inTrMS.Dispose()
        inReaderMS.Dispose()
        rawPeakReader.Dispose()
        if fallbackPeakReader.IsValueCreated then fallbackPeakReader.Value.Dispose()
        logger.Trace "Done."
