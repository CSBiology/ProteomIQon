namespace ProteomIQon

open Domain
open Dto
open MzLite
open MzLite.IO
open MzLite.Commons.Arrays
open MzLite.Binary

module Core =
    
    module MzLite = 

        module Reader = 
    
            type MzDataReader =
                | Wiff of Wiff.WiffFileReader 
                | Baf  of Bruker.BafFileReader
                | Raw  of Thermo.ThermoRawFileReader
                | MzLite of SQL.MzLiteSQL

            module MzDataReader =

                let toIMzliteDataReader reader =
                    match reader with 
                    | Wiff   r -> r :> IMzLiteDataReader 
                    | Baf    r -> r :> IMzLiteDataReader
                    | Raw    r -> r :> IMzLiteDataReader
                    | MzLite r -> r :> IMzLiteDataReader

            /// Initializes a transaction scope.
            let beginTransaction mzReader =
                mzReader 
                |> MzDataReader.toIMzliteDataReader
                |> fun x -> x.BeginTransaction()

            /// Returns the default runID used by manufacturers
            let getDefaultRunID mzReader = 
                match mzReader with
                | Wiff r -> "sample=0" 
                | Baf  r -> "run_1" 
                | Raw  r -> "run_1"
                | MzLite r -> "sample=0"
        
            let getReader (instrumentOutput:string) : Result<MzDataReader,string> = 
                match System.IO.Path.GetExtension instrumentOutput with 
                | ".wiff" -> 
                    try
                    printfn "hallo 2.1"
                    let wiffReader = new Wiff.WiffFileReader(instrumentOutput) 
                    Ok (Wiff wiffReader)
                    with ex ->
                    printfn "hallo 2.2"
                    raise ex
                    Error ex.Message 
                | ".d"    -> 
                    try
                    let bafReader = new Bruker.BafFileReader(instrumentOutput)
                    Ok (Baf bafReader)
                    with ex ->
                    Error ex.Message
                | ".mzlite" -> 
                    try
                    let mzLiteReader = new SQL.MzLiteSQL(instrumentOutput)
                    Ok (MzLite mzLiteReader)
                    with ex ->
                    Error ex.Message
                | ".RAW" -> 
                    let rawReader = new Thermo.ThermoRawFileReader(instrumentOutput)
                    Ok (Raw rawReader)
                | _       -> Error "Reader could not be opened. Only the formats .wiff (ABSciex), baf (Bruker), Raw (Thermo) or .mzlite (CSBiology) are supported." 
            
            let createMzLiteSQLDB (filePath:string)= 
                try
                let mzLiteReader = new SQL.MzLiteSQL(filePath)
                Ok mzLiteReader
                with ex ->
                Error ex.Message

            /// Returns the default runID used by manufacturers
            let getMassSpectrum mzReader specID = 
                mzReader
                |> MzDataReader.toIMzliteDataReader
                |> fun x -> x.ReadMassSpectrum(specID)

            /// Returns the default runID used by manufacturers
            let getMassSpectra mzReader runID = 
                mzReader
                |> MzDataReader.toIMzliteDataReader
                |> fun x -> x.ReadMassSpectra(runID)

        module MassSpectrum = 
            open MzLite.Model
            open System

                        /// Returns the ID of the MassSpectrum
            let getID (massSpectrum: MassSpectrum) =
                massSpectrum.ID  

            /// Returns the MsLevel of the MassSpectrum 
            let getMsLevel (massSpectrum: MassSpectrum) = 
                if massSpectrum.CvParams.Contains("MS:1000511") then 
                    (massSpectrum.CvParams.["MS:1000511"].Value) |> Convert.ToInt32
                else 
                    -1

            /// Returns the ScanTime (formerly: RetentionTime) of the MassSpectrum
            let getScanTime (massSpectrum: MassSpectrum) =  
                if massSpectrum.Scans.[0].CvParams.Contains("MS:1000016") then
                    massSpectrum.Scans.[0].CvParams.["MS:1000016"].Value |> Convert.ToDouble        
                else 
                    -1.    
    
            /// Returns PrecursorMZ of MS2 spectrum
            let getPrecursorMZ (massSpectrum: MassSpectrum) =
                if massSpectrum.Precursors.[0].SelectedIons.[0].CvParams.Contains("MS:1002234") then
                    massSpectrum.Precursors.[0].SelectedIons.[0].CvParams.["MS:1002234"].Value:?> float  // |> Convert.ToInt32        
                else 
                    -1.  

        module Peaks = 
            
            ///
            let unzipIMzliteArray (a:IMzLiteArray<Peak1D>) = 
                let mzData = Array.zeroCreate a.Length
                let intensityData = Array.zeroCreate a.Length
                for i = 0 to a.Length-1 do 
                    let peak = a.[i]
                    mzData.[i] <- peak.Mz
                    intensityData.[i] <- peak.Intensity
                mzData,intensityData

            ///
            let createIMzLiteArray mzData intensityData = 
                Array.map2 (fun mz intz -> MzLite.Binary.Peak1D(intz,mz)) mzData intensityData 
                |> MzLite.Commons.Arrays.MzLiteArray.ToMzLiteArray

            /// Creates Peak1DArray of mzData array and intensityData Array
            let createPeak1DArray compression mzBinaryDataType intensityBinaryDataType (mzData:float []) (intensityData:float []) =
                match compression with
                | true -> 
                    let peak1DArray = new Peak1DArray(BinaryDataCompressionType.ZLib,intensityBinaryDataType, mzBinaryDataType)
                    let zipedData = Array.map2 (fun mz intz -> MzLite.Binary.Peak1D(intz,mz)) mzData intensityData 
                    let newPeakA = MzLite.Commons.Arrays.MzLiteArray.ToMzLiteArray zipedData
                    peak1DArray.Peaks <- newPeakA
                    peak1DArray
                | false -> 
                    let peak1DArray = new Peak1DArray(BinaryDataCompressionType.NoCompression,intensityBinaryDataType, mzBinaryDataType)
                    let zipedData = Array.map2 (fun mz intz -> MzLite.Binary.Peak1D(intz,mz)) mzData intensityData 
                    let newPeakA = MzLite.Commons.Arrays.MzLiteArray.ToMzLiteArray zipedData
                    peak1DArray.Peaks <- newPeakA
                    peak1DArray