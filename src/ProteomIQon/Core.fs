namespace ProteomIQon

open Domain
open Dto
open System
open System.IO
open MzIO
open MzIO.IO
open MzIO.Commons.Arrays
open MzIO.Binary
//open MzIO.Bruker
open MzIO.Wiff
//open MzIO.(*Thermo*)
open MzIO.MzSQL

module Core =
    
    module MzIO = 

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
                    wiffReader :> IMzIODataReader
                //| ".d"    -> 
                //    let bafPath = Path.Combine[|instrumentOutput;"analysis.baf"|]
                //    let bafReader = new Bruker.BafFileReader(bafPath)
                //    bafReader :> IMzIODataReader
                | ".mzlite" -> 
                    let mzLiteReader = new MzSQL.MzSQL(instrumentOutput)
                    mzLiteReader :> IMzIODataReader
                //| ".raw" -> 
                //    let rawReader = new Thermo.ThermoRawFileReader(instrumentOutput)
                //    rawReader :> IMzIODataReader
                //| ".RAW" -> 
                //    let rawReader = new Thermo.ThermoRawFileReader(instrumentOutput)
                //    rawReader :> IMzIODataReader
                | _       ->  
                    failwith "Reader could not be opened. Only the formats .wiff (ABSciex), baf (Bruker), .raw (Thermo) or .mzlite (CSBiology) are supported." 
            
            /// Returns the default runID used by manufacturers
            let getDefaultRunID (mzReader:IMzIODataReader) = 
                match mzReader with
                | :? WiffFileReader as r        -> "sample=0" 
                //| :? BafFileReader as r         -> "run_1" 
                //| :? ThermoRawFileReader as r   -> "run_1"
                | :? MzSQL as r                 -> "sample=0"
                //| :? Thermo.ThermoRawFileReader as r -> "run_1"
            /// Initializes a transaction scope.
            let beginTransaction (mzReader:IMzIODataReader) =
                mzReader.BeginTransaction()


        module Peaks = 

            ///
            let unzipIMzliteArray (a:IMzIOArray<Peak1D>) = 
                let mzData = Array.zeroCreate a.Length
                let intensityData = Array.zeroCreate a.Length
                for i = 0 to a.Length-1 do 
                    let peak = a.[i]
                    mzData.[i] <- peak.Mz
                    intensityData.[i] <- peak.Intensity
                mzData,intensityData


