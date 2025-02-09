namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
open System.Reflection
open ProteomIQon.Core.InputPaths
open ProteomIQon.Core
open ProteomIQon

module console1 =

    [<EntryPoint>]
    let main argv =
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult InstrumentOutput |> List.map getPathRelativeToDir
        let o = results.GetResult OutputDirectory  |> getPathRelativeToDir
        let p = results.GetResult ParamFile        |> getPathRelativeToDir
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "MzMLToMzLite"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let processParams =
            Json.ReadAndDeserialize<Dto.PreprocessingParams> p
            |> Dto.PreprocessingParams.toDomain

        let files = 
            parsePaths (fun path -> MzIO.Reader.getMzMLFiles path) i
            |> Array.ofSeq

        if files.Length = 1  then
            logger.Info "single file"
            logger.Trace (sprintf "Preprocessing %s" files.[0])
            MzMLIonMobilityToMzLite.processFile processParams o files.[0]
        else
            logger.Info "multiple files"
            logger.Trace (sprintf "Preprocessing multiple files: %A" files)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            try
            let partitionedFiles =
                files
                |> Array.splitInto c
            [for i in partitionedFiles do yield async { return i |> Array.map (MzMLIonMobilityToMzLite.processFile processParams o)}]
            |> Async.Parallel
            |> Async.RunSynchronously
            |> ignore
            with
            | ex -> printfn "%A" ex
        logger.Info "Done"
        0