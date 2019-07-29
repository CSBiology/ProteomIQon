namespace ProteomIQon

open Domain
open Dto
open System
open System.IO
open MzLite
open MzLite.IO
open MzLite.Commons.Arrays
open MzLite.Binary
open MzLite.Bruker
open MzLite.Wiff
open MzLite.Thermo
open MzLite.SQL

module Core =
    
    module MzLite = 

        module Reader = 
            let getMzLiteFiles directoryPath = 
                Directory.GetFiles(directoryPath,("*.mzlite"))   

            let getThermoRawFiles directoryPath = 
                Array.append (Directory.GetFiles(directoryPath,("*.raw"))) (Directory.GetFiles(directoryPath,("*.RAW")))   
                |> Array.distinct
            let getWiffFiles directoryPath = 
                Directory.GetFiles(directoryPath,("*.wiff"))   

            let getBrukerFiles directoryPath = 
                Directory.GetDirectories(directoryPath,("*.d"))   

            let getMSFilePaths directoryPath =
                [|
                    getMzLiteFiles directoryPath
                    getThermoRawFiles directoryPath
                    getWiffFiles directoryPath
                    getBrukerFiles directoryPath
                |]
                |> Array.concat
            let getReader (instrumentOutput:string) = 
                match System.IO.Path.GetExtension instrumentOutput with 
                | ".wiff" -> 
                    let wiffReader = new Wiff.WiffFileReader(instrumentOutput) 
                    wiffReader :> IMzLiteDataReader
                | ".d"    -> 
                    let bafPath = Path.Combine[|instrumentOutput;"analysis.baf"|]
                    let bafReader = new Bruker.BafFileReader(bafPath)
                    bafReader :> IMzLiteDataReader
                | ".mzlite" -> 
                    let mzLiteReader = new SQL.MzLiteSQL(instrumentOutput)
                    mzLiteReader :> IMzLiteDataReader
                | ".raw" -> 
                    let rawReader = new Thermo.ThermoRawFileReader(instrumentOutput)
                    rawReader :> IMzLiteDataReader
                | ".RAW" -> 
                    let rawReader = new Thermo.ThermoRawFileReader(instrumentOutput)
                    rawReader :> IMzLiteDataReader
                | _       ->  
                    failwith "Reader could not be opened. Only the formats .wiff (ABSciex), baf (Bruker), .raw (Thermo) or .mzlite (CSBiology) are supported." 
            
            /// Returns the default runID used by manufacturers
            let getDefaultRunID (mzReader:IMzLiteDataReader) = 
                match mzReader with
                | :? WiffFileReader as r        -> "sample=0" 
                | :? BafFileReader as r         -> "run_1" 
                | :? ThermoRawFileReader as r   -> "run_1"
                | :? MzLiteSQL as r             -> "sample=0"
                | _                             -> failwith "There is no known default RunID for the given Reader"
            /// Initializes a transaction scope.
            let beginTransaction (mzReader:IMzLiteDataReader) =
                mzReader.BeginTransaction()

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
        module Query = 
            
            open MzLite.Processing
            open MzLite.Wiff
            open MzLite.IO
            
            ///
            let createRangeQuery v offset =
                new RangeQuery(v, offset)
            
            ///
            let getMS1RTIdx (reader:IMzLiteDataReader) runId = 
                reader.BuildRtIndex(runId)

            /// 
            let getXIC (reader:IMzLiteDataReader) (rtIdx:MzLite.Commons.Arrays.IMzLiteArray<MzLiteLinq.RtIndexEntry>) (rtQuery:RangeQuery) (mzQuery:RangeQuery) = 
                reader.RtProfile(rtIdx, rtQuery, mzQuery) 

            ///
            let getXICs (reader:IMzLiteDataReader) (rtIdx:MzLite.Commons.Arrays.IMzLiteArray<MzLiteLinq.RtIndexEntry>) (rtQuery:RangeQuery) (mzQueries:RangeQuery []) = 
                reader.RtProfiles(rtIdx, rtQuery, mzQueries) 
               
            ///       
            let createSwathQuery targetMz rtQuery ms2MzQueries =
                new SwathQuery(targetMz, rtQuery, ms2MzQueries)

            ///
            let getSwathIdx (reader:IMzLiteDataReader) runId =
                SwathIndexer.Create(reader, runId)

            ///
            let getSwathXics (reader:IMzLiteDataReader) (swathIdx:SwathIndexer) swathQuery = 
                swathIdx.GetMS2(reader, swathQuery)

            ///        
            let getSwathXICsBy (reader:IMzLiteDataReader) (swathIdx:SwathIndexer) (rtQuery:RangeQuery) (ms2MzQueries:RangeQuery []) tarMz = 
                let swathQ = createSwathQuery tarMz rtQuery ms2MzQueries
                getSwathXics reader swathIdx swathQ
