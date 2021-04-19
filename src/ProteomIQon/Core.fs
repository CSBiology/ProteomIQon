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
//open MzIO.Wiff
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
                //| ".mzml" -> 
                    //let mzmlReader = new mzmlReader.WiffFileReader(instrumentOutput) 
                    //wiffReader :> IMzIODataReader
                | ".mzlite" -> 
                    let mzLiteReader = new MzSQL.MzSQL(instrumentOutput)
                    mzLiteReader :> IMzIODataReader
                | _       ->  
                    failwith "Reader could not be opened. Only the formats .mzml or .mzlite (CSBiology) are supported." 
            
            /// Returns the default runID used by manufacturers
            let getDefaultRunID (mzReader:IMzIODataReader) = 
                match mzReader with
                | :? MzSQL as r                 -> "sample=0"

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

    module InputPaths =
        
        open FSharpAux

        let getRelativePath (reference: string) (path: string) =
            Path.Combine(reference, path)
        
        let parsePath (extension: string) fileOrDirectoryPath =
            if File.Exists fileOrDirectoryPath then
                [|fileOrDirectoryPath|]
            elif Directory.Exists fileOrDirectoryPath then
                let fps = Directory.GetFiles(fileOrDirectoryPath ,(sprintf "*.%s" extension))
                fps
            else 
                [||]

        let parsePaths (extension: string) (fileOrDirectoryPaths:seq<string>) =
            fileOrDirectoryPaths
            |> Seq.collect (parsePath extension)
