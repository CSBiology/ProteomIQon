namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PSMStatistics
open System.Reflection
open ProteomIQon.Core.InputPaths
open System

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv = 
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)       
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse(argv)
        //let results = parser.ParseCommandLine(inputs = argv, raiseOnUsage = false)
        let i = results.GetResult PSMs            |> List.map getPathRelativeToDir
        let o = results.GetResult OutputDirectory |> getPathRelativeToDir
        let p = results.GetResult ParamFile       |> getPathRelativeToDir
        let d = results.GetResult PeptideDataBase |> getPathRelativeToDir
        Logging.generateConfig o
        let logger = Logging.createLogger "PSMStatistics"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Info (sprintf "Peptide data base -d = %s" d)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let p = 
            Json.ReadAndDeserialize<Dto.PSMStatisticsParams> p
            |> Dto.PSMStatisticsParams.toDomain
        let dbConnection =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        let files = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.psm"))) i
            |> Array.ofSeq
        if files.Length = 1 then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Processing %s" files.[0])
            psmStats p o dbConnection files.[0]
        else
            logger.Info (sprintf "multiple files")
            logger.Trace (sprintf "Processing multiple files: %A" files)
            let c = 
                match results.TryGetResult Parallelism_Level with 
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            files 
            |> FSharpAux.PSeq.map (psmStats p o dbConnection) 
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        dbConnection.Dispose()
        logger.Info "Done"
        0
