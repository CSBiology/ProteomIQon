namespace ProteomIQon

open Domain
open Core
open System.IO
open BioFSharp.Mz
open MzIO
open MzIO.Binary
open MzIO.Model
open MzIO.Processing
open MzIO.IO
open MzIO.MzSQL
open MzIO.IO.MzML
open ProteomIQon.Core.MzIO.Processing

module MzMLToMzLite =

    ///
    let private initPeakPicking (reader: MzMLReader) (peakPickingParams:PeakPicking) (outputDir:string) (instrumentOutput:string)=

        //outputDir and instrumentOutput only added for logger
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension instrumentOutput)
        let getP1D (id: string)=
            reader.getSpecificPeak1DArraySequential(id)
        match peakPickingParams with
        | PeakPicking.ProfilePeaks ->
            fun (massSpec:MassSpectrum) ->
                (getP1D massSpec.ID).Peaks
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
                    fun (massSpec:MassSpectrum) ->
                        let mzData, intensityData =
                            (getP1D massSpec.ID).Peaks
                            |> Core.MzIO.Peaks.unzipIMzliteArray
                        let paddedMz,paddedIntensity =
                            SignalDetection.Padding.paddDataBy paddingParams mzData intensityData
                        BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters paddedMz paddedIntensity
                | YThreshold.MinSpectrumIntensity ->
                    logger.Trace (sprintf "ms2Centroidization with: %A" waveletParams)
                    fun (massSpec:MassSpectrum) ->
                        let mzData, intensityData =
                            (getP1D massSpec.ID).Peaks
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
                    logger.Trace (sprintf "ms1Centroidization with: %A" waveletParams)
                    let waveletParameters = initwaveletParameters yThreshold
                    fun (massSpec:MassSpectrum) ->
                        let mzData, intensityData =
                            (getP1D massSpec.ID).Peaks
                            |> Core.MzIO.Peaks.unzipIMzliteArray
                        BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters mzData intensityData
                | YThreshold.MinSpectrumIntensity ->
                    fun (massSpec:MassSpectrum) ->
                        let mzData, intensityData =
                            (getP1D massSpec.ID).Peaks
                            |> Core.MzIO.Peaks.unzipIMzliteArray
                        let yThreshold = Array.min intensityData
                        let waveletParameters = initwaveletParameters yThreshold
                        BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters mzData intensityData
        | PeakPicking.Centroid CentroidizationMode.Manufacturer ->
            failwith "Manufacturer peak picking is only supported for .baf (Bruker) files."

    ///
    let insertSprectrum (compress:BinaryDataCompressionType) (outReader: MzSQL.MzSQL) (runID:string)
        (ms1PeakPicking: MassSpectrum -> float [] * float []) (ms2PeakPicking: MassSpectrum -> float [] * float [])
            (spectrum: MassSpectrum) =
        match MassSpectrum.getMsLevel spectrum with
        | 1 ->
            let mzData,intensityData =
                try
                ms1PeakPicking spectrum
                with
                | _ -> [||],[||]
            if Array.isEmpty mzData || Array.isEmpty intensityData then ()
            else
            let peaks = PeakArray.createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData 
            outReader.Insert(runID, spectrum, peaks)
        | 2 ->
            let mzData,intensityData =
                try
                ms2PeakPicking spectrum
                with
                | _ -> [||],[||]
            if Array.isEmpty mzData || Array.isEmpty intensityData then ()
            else
            let peaks = PeakArray.createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData
            outReader.Insert(runID, spectrum, peaks)
        | _ ->
            failwith "Only mass spectra of level 1 and 2 are supported."

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

        logger.Trace "Init connection to input data base."
        // initialize Reader and Transaction
        let inReaderMS = new MzMLReader(instrumentOutput)
        let inReaderPeaks = new MzMLReader(instrumentOutput)
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReaderMS
        let inTrMS = inReaderMS.BeginTransaction()
        let inTrPeaks = inReaderPeaks.BeginTransaction()

        logger.Trace "Creating mzlite file."
        // initialize Reader and Transaction
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".mzlite"
            Path.Combine [|outputDir;fileName|]

        let outReader = new MzSQL(outFilePath)
        let cn = outReader.Open()
        /// All files created by this application will have a unified runID.
        let outRunID  = Core.MzIO.Reader.getDefaultRunID outReader
        let outTr = outReader.BeginTransaction()
        logger.Trace "Try inserting Model."
        try
            outReader.InsertModel inReaderMS.Model
            logger.Trace "Model inserted."
        with
        | ex -> logger.Trace $"Inserting model failed: {ex}"
        inReaderMS.ResetReader()
        logger.Trace "Initiating peak picking functions."
        // Initialize PeakPickingFunctions
        let ms1PeakPicking = initPeakPicking inReaderPeaks processParams.MS1PeakPicking outputDir instrumentOutput
        let ms2PeakPicking = initPeakPicking inReaderPeaks processParams.MS2PeakPicking outputDir instrumentOutput

        logger.Trace "Getting mass spectra."
        // Get all mass spectra
        let massSpectra = 
            inReaderMS.ReadMassSpectra(inRunID)
        logger.Trace "Filtering mass spectra according to retention time."
        // Filter mass spectra by minimum or maximum scan time time.
        let massSpectraF =
            match processParams.StartRetentionTime, processParams.EndRetentionTime with
            | Some s, Some e ->
                massSpectra
                |> Seq.filter (fun ms ->
                                let scanTime = MassSpectrum.getScanTime ms
                                scanTime > s && scanTime < e
                              )
            | Some s, None ->
                massSpectra
                |> Seq.filter (fun ms ->
                                let scanTime = MassSpectrum.getScanTime ms
                                scanTime > s
                              )
            | None, Some e ->
                massSpectra
                |> Seq.filter (fun ms ->
                                let scanTime = MassSpectrum.getScanTime ms
                                scanTime < e
                              )
            | None, None ->
                massSpectra
        //logger.Trace (sprintf "Copying %i mass spectra to output data base." (Seq.length massSpectraF))
        ///
        let insert =
            massSpectraF
            |> Seq.filter (fun ms ->
                                let level = MassSpectrum.getMsLevel ms
                                level = 1 || level = 2
                          )
            |> Seq.iteri (fun i ms ->
                            if i%1000 = 0 then logger.Trace (sprintf "%i" i)
                            try
                                insertSprectrum processParams.Compress outReader outRunID ms1PeakPicking ms2PeakPicking (changeScanTimeToMinutes ms)
                            with
                            | ex ->
                                logger.Trace (sprintf "File:%s ID: %s could not be inserted. Exeption:%A" instrumentOutput ms.ID ex)
                        )
        insert
        inTrMS.Commit()
        inTrMS.Dispose()
        inTrPeaks.Commit()
        inTrPeaks.Dispose()
        inReaderMS.Dispose()
        inReaderPeaks.Dispose()
        outTr.Commit()
        outTr.Dispose()
        outReader.Dispose()
        logger.Trace "Done."