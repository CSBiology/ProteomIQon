namespace ProteomIQon

open Domain
open Core 
open Logary
open System.IO
open MzLite 
open BioFSharp.Mz
open MzLite.Model
open MzLite.Binary
open Core.MzLite.Reader
open Core.MzLite.Peaks

module Preprocessing = 

    /// 
    let private initPeakPicking (reader:MzDataReader) (peakPickingParams:PeakPicking) = 
        match reader, peakPickingParams with
        | MzDataReader.Baf r, PeakPicking.Centroid CentroidizationMode.Manufacturer -> 
            fun (massSpec:MassSpectrum) -> 
                r.ReadSpectrumPeaks(massSpec.ID,true).Peaks
                |> unzipIMzliteArray
        | _ as r , PeakPicking.ProfilePeaks ->
            fun (massSpec:MassSpectrum) -> 
                (MzDataReader.toIMzliteDataReader r).ReadSpectrumPeaks(massSpec.ID).Peaks
                |> unzipIMzliteArray
        | _ as r , PeakPicking.Centroid (CentroidizationMode.Wavelet waveletParams) -> 
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
                match waveletParams.YThreshold with 
                | YThreshold.Fixed yThreshold -> 
                    let paddingParams = initPaddingParameters yThreshold
                    let waveletParameters = initwaveletParameters yThreshold
                    fun (massSpec:MassSpectrum) -> 
                        let mzData, intensityData = 
                            (MzDataReader.toIMzliteDataReader r).ReadSpectrumPeaks(massSpec.ID).Peaks
                            |> unzipIMzliteArray
                        let paddedMz,paddedIntensity = 
                            let paddingParams = initPaddingParameters yThreshold
                            SignalDetection.Padding.paddDataBy paddingParams mzData intensityData
                        BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters paddedMz paddedIntensity
                | YThreshold.MinSpectrumIntensity -> 
                    fun (massSpec:MassSpectrum) -> 
                        let mzData, intensityData = 
                            (MzDataReader.toIMzliteDataReader r).ReadSpectrumPeaks(massSpec.ID).Peaks
                            |> unzipIMzliteArray
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
                match waveletParams.YThreshold with 
                | YThreshold.Fixed yThreshold -> 
                    let waveletParameters = initwaveletParameters yThreshold
                    fun (massSpec:MassSpectrum) -> 
                        let mzData, intensityData = 
                            (MzDataReader.toIMzliteDataReader r).ReadSpectrumPeaks(massSpec.ID).Peaks
                            |> unzipIMzliteArray
                        BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters mzData intensityData
                | YThreshold.MinSpectrumIntensity -> 
                    fun (massSpec:MassSpectrum) -> 
                        let mzData, intensityData = 
                            (MzDataReader.toIMzliteDataReader r).ReadSpectrumPeaks(massSpec.ID).Peaks
                            |> unzipIMzliteArray
                        let yThreshold = Array.min intensityData 
                        let waveletParameters = initwaveletParameters yThreshold
                        BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters mzData intensityData
        | _ as r , PeakPicking.Centroid CentroidizationMode.Manufacturer ->
            failwith "Manufacturer peak picking is only supported for .baf (Bruker) files."               
                    
    ///                    
    let insertSprectrum (compress:bool) (outReader: MzLite.SQL.MzLiteSQL) (runID:string) 
        (ms1PeakPicking: MassSpectrum -> float [] * float []) (ms2PeakPicking: MassSpectrum -> float [] * float [])
            (spectrum: MassSpectrum) =
        match MassSpectrum.getMsLevel spectrum with 
        | 1 -> 
            let mzData,intensityData = 
                ms1PeakPicking spectrum
            let peaks = createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData
            outReader.Insert(runID, spectrum, peaks)           
        | 2 -> 
            let mzData,intensityData = 
                ms2PeakPicking spectrum
            let peaks = createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData
            outReader.Insert(runID, spectrum, peaks)
        | _ -> 
            failwith "Only mass spectra of level 1 and 2 are supported."

    let processFile (processParams:PreprocessingParams) (outputDir:string) (instrumentOutput:string) =
        printfn "Hallo"
        let logger = Log.create()
        Message.event Info (sprintf "Now preprocessing: %s /nResults will be written to: %s" instrumentOutput outputDir) |> logger.logSimple
        printfn "Hallo2"
        // initialize Reader and Transaction
        let inReader = 
            match Core.MzLite.Reader.getReader instrumentOutput with 
            | Result.Ok reader -> 
                Message.event Info ("Reader initialized successfully") |> logger.logSimple                
                reader
            | Result.Error ex -> 
                failwith ex

        let inRunID = Core.MzLite.Reader.getDefaultRunID inReader

        let inTr = 
            let tmp = Core.MzLite.Reader.beginTransaction inReader                    
            Message.event Info ("Transaction Scope initialized") |> logger.logSimple
            tmp 
        printfn "Hallo3"

        // initialize Reader and Transaction
        let outFilePath = 
            let fileName = Path.GetFileName instrumentOutput
            Path.Combine [|outputDir;fileName|]
            
        let outReader = 
            match Core.MzLite.Reader.createMzLiteSQLDB outFilePath with 
            | Result.Ok reader -> 
                Message.event Info ("Output database initialized successfully") |> logger.logSimple                
                reader
            | Result.Error ex -> 
                failwith ex
        printfn "Hallo4"
        // TODO: Factor this out.
        /// All files created by this application will have a unified runID.
        let outRunID = "sample=0"

        let outTr = 
            let tmp = outReader.BeginTransaction()                    
            Message.event Info ("Transaction Scope initialized") |> logger.logSimple
            tmp 
        printfn "Hallo5"
        // Initialize PeakPickingFunctions
        let ms1PeakPicking = initPeakPicking inReader processParams.MS1PeakPicking
        let ms2PeakPicking = initPeakPicking inReader processParams.MS2PeakPicking

        Message.event Info ("Getting all mass spectra") |> logger.logSimple
        // Get all mass spectra 
        let massSpectra = 
            getMassSpectra inReader inRunID

        printfn "Hallo6"
        Message.event Info ("Filtering mass spectra by scan time.") |> logger.logSimple
        // Filter mass spectra by minimum or maximum scan time time.            
        let massSpectraF = 
            match processParams.StartRetentionTime, processParams.EndRetentionTime with
            | Some s, Some e ->
                massSpectra
                |> Seq.filter (fun ms -> 
                                let scanTime = Core.MzLite.MassSpectrum.getScanTime ms
                                scanTime > s && scanTime < e
                              )
            | Some s, None ->
                massSpectra
                |> Seq.filter (fun ms -> 
                                let scanTime = Core.MzLite.MassSpectrum.getScanTime ms
                                scanTime > s 
                              )
            | None, Some e ->
                massSpectra
                |> Seq.filter (fun ms -> 
                                let scanTime = Core.MzLite.MassSpectrum.getScanTime ms
                                scanTime < e
                              )
            | None, None ->
                massSpectra
        
        printfn "Hallo7"
        Message.event Info ("processing and copying mass spectra into output data base.") |> logger.logSimple
        ///  
        massSpectraF
        |> Seq.iter (insertSprectrum processParams.Compress outReader outRunID ms1PeakPicking ms2PeakPicking)
        
        Message.event Info ("Done") |> logger.logSimple
        
