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
open MzLite.Bruker

module Preprocessing = 
    open MzLite.IO
    open MzLite.SQL

    /// 
    let private initPeakPicking (reader:IMzLiteDataReader) (peakPickingParams:PeakPicking) = 
        match reader, peakPickingParams with
        | :? BafFileReader as r, PeakPicking.Centroid CentroidizationMode.Manufacturer -> 
            fun (massSpec:MassSpectrum) -> 
                r.ReadSpectrumPeaks(massSpec.ID,true).Peaks
                |> unzipIMzliteArray
        | _ as r , PeakPicking.ProfilePeaks ->
            fun (massSpec:MassSpectrum) -> 
                r.ReadSpectrumPeaks(massSpec.ID).Peaks
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
                            r.ReadSpectrumPeaks(massSpec.ID).Peaks
                            |> unzipIMzliteArray
                        let paddedMz,paddedIntensity = 
                            SignalDetection.Padding.paddDataBy paddingParams mzData intensityData
                        BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters paddedMz paddedIntensity
                | YThreshold.MinSpectrumIntensity -> 
                    printfn "ms2Centroidization with: %A" waveletParams
                    fun (massSpec:MassSpectrum) -> 
                        let mzData, intensityData = 
                            r.ReadSpectrumPeaks(massSpec.ID).Peaks
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
                    printfn "ms1Centroidization with: %A" waveletParams
                    let waveletParameters = initwaveletParameters yThreshold
                    fun (massSpec:MassSpectrum) -> 
                        let mzData, intensityData = 
                            r.ReadSpectrumPeaks(massSpec.ID).Peaks
                            |> unzipIMzliteArray
                        BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters mzData intensityData
                | YThreshold.MinSpectrumIntensity -> 
                    fun (massSpec:MassSpectrum) -> 
                        let mzData, intensityData = 
                            r.ReadSpectrumPeaks(massSpec.ID).Peaks
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
                try
                ms1PeakPicking spectrum
                with 
                | _ -> [||],[||]
            let peaks = createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData
            outReader.Insert(runID, spectrum, peaks)           
        | 2 -> 
            let mzData,intensityData = 
                try
                ms2PeakPicking spectrum
                with 
                | _ -> [||],[||]
            let peaks = createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData
            outReader.Insert(runID, spectrum, peaks)
        | _ -> 
            failwith "Only mass spectra of level 1 and 2 are supported."

    let processFile (processParams:PreprocessingParams) (outputDir:string) (instrumentOutput:string) =
        printfn "Now preprocessing: %s Results will be written to: %s" instrumentOutput outputDir
        
        printfn "Init connection to input data base." 
        // initialize Reader and Transaction
        let inReader = Core.MzLite.Reader.getReader instrumentOutput  
        let inRunID  = Core.MzLite.Reader.getDefaultRunID inReader
        let inTr = inReader.BeginTransaction()                    

        printfn "Creating mzlite file." 
        // initialize Reader and Transaction
        let outFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".mzlite"
            Path.Combine [|outputDir;fileName|]
            
        let outReader = new MzLiteSQL(outFilePath) 
        /// All files created by this application will have a unified runID.
        let outRunID  = Core.MzLite.Reader.getDefaultRunID outReader    
        let outTr = outReader.BeginTransaction()                    

        printfn "Initiating peak picking functions."
        // Initialize PeakPickingFunctions
        let ms1PeakPicking = initPeakPicking inReader processParams.MS1PeakPicking
        let ms2PeakPicking = initPeakPicking inReader processParams.MS2PeakPicking
 
        printfn "Getting mass spectra."
        // Get all mass spectra 
        let massSpectra = inReader.ReadMassSpectra(inRunID)

        printfn "Filtering mass spectra according to retention time."
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
        
        printfn "Copying %i mass spectra to output data base." (Seq.length massSpectraF) 
        ///  
        massSpectraF
        |> Seq.filter (fun ms -> 
                            let level = Core.MzLite.MassSpectrum.getMsLevel ms 
                            level = 1 || level = 2
                      )
        |> Seq.iter (fun ms -> 
                        try
                            insertSprectrum processParams.Compress outReader outRunID ms1PeakPicking ms2PeakPicking ms
                        
                        with 
                        | ex -> 
                            printfn "File:%s ID: %s could not be inserted. Exeption:%A" instrumentOutput ms.ID ex
                    )
        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()
        outTr.Commit()
        outTr.Dispose()
        outReader.Dispose()
        printfn "Done."
