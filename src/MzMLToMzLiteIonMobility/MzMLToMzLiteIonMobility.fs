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
        let inRunID  = getDefaultRunID inReaderMS
        let inTrMS = inReaderMS.BeginTransaction()
        // the peaks are read in file order by a second reader while inReaderMS reads the spectra
        let inReaderPeaks = new MzMLReaderMIRIM(instrumentOutput)

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
        spectra
        |> Seq.iteri (fun i spectrum ->
            if i % 1000 = 0 then logger.Trace $"{i}"
            let spectrum' = changeScanTimeToMinutes (fixSpectrum spectrum)
            let peaks = inReaderPeaks.ReadPeakArraysSequential spectrum.ID
            match MassSpectrum.getMsLevel spectrum' with
            | 1 | 2 -> outReader.Insert(outRunID, spectrum', processParams.Compress, peaks.Mz, peaks.Intensity, peaks.IonMobility, peaks.Count)
            | _ -> failwith "Only mass spectra of level 1 and 2 are supported."
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
        inReaderPeaks.Dispose()
        logger.Trace "Done."
