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
open MzIO.Model
open MzIO.MetaData.PSIMSExtension
open MzIO.Model.CvParam

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

            let getMzMLFiles directoryPath = 
                [|
                    Directory.GetFiles(directoryPath,("*.mzML"))
                    Directory.GetFiles(directoryPath,("*.mzml"))
                |]
                |> Array.concat
                |> Array.distinct

            let getMSFilePaths directoryPath =
                [|
                    getMzLiteFiles directoryPath
                    getThermoRawFiles directoryPath
                    getWiffFiles directoryPath
                    getBrukerFiles directoryPath
                    getMzMLFiles directoryPath
                |]
                |> Array.concat

            let getMzLiteMzMLPaths directoryPath =
                [|
                    getMzLiteFiles directoryPath
                    getMzMLFiles directoryPath
                |]
                |> Array.concat

            let getReader (instrumentOutput:string) = 
                match System.IO.Path.GetExtension instrumentOutput with 
                | ".mzML" -> 
                    let mzmlReader = new MzIO.IO.MzML.MzMLReader(instrumentOutput) 
                    mzmlReader :> IMzIODataReader
                | ".mzlite" -> 
                    let mzLiteReader = new MzSQL.MzSQL(instrumentOutput)
                    mzLiteReader :> IMzIODataReader
                | _       ->  
                    failwith "Reader could not be opened. Only the formats .mzml or .mzlite (CSBiology) are supported." 
            
            /// Returns the default runID used by manufacturers
            let getDefaultRunID (mzReader:IMzIODataReader) = 
                match mzReader with
                | :? MzSQL as r                 -> "sample=0"
                | :? MzML.MzMLReader as r       -> "sample=0"

            /// Initializes a transaction scope.
            let beginTransaction (mzReader:IMzIODataReader) =
                mzReader.BeginTransaction()

            /// Opens a connection.
            let openConnection (mzReader: IMzIODataReader) =
                match mzReader with
                | :? MzIO.MzSQL.MzSQL as r -> r.Connection.Open()
                | _ -> ()

        module Processing =
            
            /// Changes the unit for ScanTime (formerly: RetentionTime) of the MassSpectrum to minutes
            let changeScanTimeToMinutes (massSpectrum: MassSpectrum) =
                let newScanList = new ScanList()
                let scans =  
                    massSpectrum.Scans.GetProperties false
                    |> Array.ofSeq
                    |> Array.map (fun scan -> 
                        match scan.Value with
                        | :? Scan -> 
                            let s = (scan.Value :?> Scan)
                            let rt = s.TryGetValue(PSIMS_Scan.ScanStartTime)
                            match rt with
                            | Some t -> 
                                let oldParam = t :?> IParamBase<IConvertible>
                                let oldParamValue =
                                    Convert.ToDouble((tryGetValue oldParam).Value)
                                let oldParamUnit =
                                    (tryGetCvUnitAccession oldParam).Value
                                match oldParamUnit with
                                | "UO:0000010" ->
                                    let newParam =
                                        let timeInMin = oldParamValue / 60.
                                        let cvParam = new CvParam<IConvertible>(PSIMS_Scan.ScanStartTime, ParamValue.WithCvUnitAccession((timeInMin :> IConvertible), "UO:0000031"))
                                        cvParam
                                    s.RemoveItem(PSIMS_Scan.ScanStartTime)
                                    s.AddCvParam(newParam)
                                    newScanList.Add(Guid.NewGuid().ToString(),s)
                                | "UO:0000031" -> newScanList.Add(Guid.NewGuid().ToString(),s)
                                | "UO:0000032" ->
                                    let newParam =
                                        let timeInMin = oldParamValue * 60.
                                        let cvParam = new CvParam<IConvertible>(PSIMS_Scan.ScanStartTime, ParamValue.WithCvUnitAccession((timeInMin :> IConvertible), "UO:0000031"))
                                        cvParam
                                    s.RemoveItem(PSIMS_Scan.ScanStartTime)
                                    s.AddCvParam(newParam)
                                    newScanList.Add(Guid.NewGuid().ToString(),s)
                                | _ -> failwith "Unknown Time Unit in Scans"

                            | None -> newScanList.Add(Guid.NewGuid().ToString(),s)
                        | :? UserParam<IConvertible> ->
                            newScanList.AddUserParam(scan.Value :?> UserParam<IConvertible>)
                        | :? CvParam<IConvertible> ->
                            newScanList.AddCvParam(scan.Value :?> CvParam<IConvertible>)
                        | _ -> failwith "unexpected input type"
                    )
                massSpectrum.Scans <- newScanList
                massSpectrum


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
        
        let parsePath (getFiles: string -> string[]) fileOrDirectoryPath =
            if File.Exists fileOrDirectoryPath then
                [|fileOrDirectoryPath|]
            elif Directory.Exists fileOrDirectoryPath then
                let fps = getFiles fileOrDirectoryPath
                fps
            else 
                [||]

        let parsePaths (getFiles: string -> string[]) (fileOrDirectoryPaths:seq<string>) =
            fileOrDirectoryPaths
            |> Seq.collect (parsePath getFiles)

    module Zipping =
    
        open System.IO
        open System.IO.Compression

        let zipDirectory (globPattern) (logger: NLog.Logger) (directoryName: string) =
            let files = Directory.GetFiles(directoryName, globPattern)
            try
                use ms = new MemoryStream()
                (
                    use archive = new ZipArchive(ms, ZipArchiveMode.Create)
                    files
                    |> Array.iter (fun fileName ->
                        let data = File.ReadAllBytes fileName
                        let entry = archive.CreateEntry(Path.GetFileName fileName)
                        use entryStream = entry.Open()
                        use bw = new BinaryWriter(entryStream)
                        bw.Write(data)
                    )
                )
                Ok (ms.ToArray())
            with e ->
                Error <| logger.Trace (sprintf "Cannot zip stream %s: %s" directoryName e.Message)

        let saveZippedDirectory path (logger: NLog.Logger) fileName (data: byte []) =
            try
                let path = Path.Combine(path, fileName + ".zip")
                if File.Exists path then
                    use fs = File.Open(Path.Combine(path), FileMode.Append)
                    fs.Write(data, 0, data.Length)
                    Ok ()
                else
                    use fs = File.Open(Path.Combine(path), FileMode.OpenOrCreate)
                    fs.Write(data, 0, data.Length)
                    Ok ()
            with e ->
                Error <| logger.Trace (sprintf "Cannot save directory %s: %s" fileName e.Message)
