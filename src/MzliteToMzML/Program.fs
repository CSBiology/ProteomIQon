namespace ProteomIQon

open System
open System.IO
open System.Reflection
open Argu
open ProteomIQon.Core
open ProteomIQon.Core.InputPaths
open CLIArgumentParsing

module console1 =

    [<EntryPoint>]
    let main argv = 
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult InstrumentOutput |> getPathRelativeToDir
        let o = results.GetResult OutputDirectory  |> getPathRelativeToDir
        let p = results.GetResult ParamFile        |> getPathRelativeToDir
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "mzMLConverter"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let converterParams =
                Json.ReadAndDeserialize<Dto.MzMLConverterParams> p
                |> Dto.MzMLConverterParams.toDomain
        if File.Exists i || (Directory.Exists i && i.EndsWith(".d"))  then
            logger.Info "single file"
            MzliteToMzML.convertFile converterParams o i
            logger.Trace (sprintf "Converting %s" i)
        elif Directory.Exists i then
            logger.Info "multiple files"
            let files =
                Core.MzIO.Reader.getMSFilePaths i
                |> Array.sort
            logger.Trace (sprintf "Converting multiple files: %A" files)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            try
            let partitionedFiles =
                files
                |> Array.splitInto c
            [for i in partitionedFiles do yield async { return i |> Array.map (MzliteToMzML.convertFile converterParams o)}]
            |> Async.Parallel
            |> Async.RunSynchronously
            |> ignore
            with
            | ex -> printfn "%A" ex
        else
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        logger.Info "Done"
        0
