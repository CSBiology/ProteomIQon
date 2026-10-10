namespace ProteomIQon

open System
open System.IO
open Core.InputPaths
open QuantBasedAlignment_v2
open CLIArgumentParsing
open Argu

module console1 =

    /// Logging.generateConfig sends Trace and Warn to the log file and Info to the console.
    /// Also send Debug, Info, Error and Fatal to the log file and warnings and errors to the console.
    let extendLogging () =
        let config = NLog.LogManager.Configuration
        config.AddRule(NLog.LogLevel.Debug, NLog.LogLevel.Info, config.FindTargetByName "file")
        config.AddRule(NLog.LogLevel.Error, NLog.LogLevel.Fatal, config.FindTargetByName "file")
        config.AddRule(NLog.LogLevel.Warn, NLog.LogLevel.Fatal, config.FindTargetByName "console")
        NLog.LogManager.Configuration <- config
        NLog.LogManager.ReconfigExistingLoggers()

    [<EntryPoint>]
    let main argv =
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i   = results.GetResult TargetFiles |> List.map getPathRelativeToDir
        let ii  = results.GetResult SourceFiles |> List.map getPathRelativeToDir
        let o   = results.GetResult OutputDirectory    |> getPathRelativeToDir
        let dc  = results.Contains DiagnosticCharts
        let seed = results.TryGetResult Seed |> Option.defaultValue 1
        let c = results.TryGetResult Parallelism_Level |> Option.defaultValue 1
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        extendLogging ()
        let logger = Logging.createLogger "QuantBasedAlignment_v2"
        logger.Info (sprintf "TargetFiles -i = %A" i)
        logger.Info (sprintf "SourceFiles -ii = %A" ii)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "Seed = %i" seed)
        logger.Trace (sprintf "CLIArguments: %A" results)
        try
            let quantFiles (paths: string list) =
                for path in paths do
                    if not (File.Exists path || Directory.Exists path) then
                        failwithf "%s does not exist." path
                let files =
                    parsePaths (fun path -> Directory.GetFiles(path,("*.quant")) |> Array.sort) paths
                    |> Array.ofSeq
                    |> Array.map Path.GetFullPath
                if files.Length = 0 then
                    failwithf "No .quant files in %s." (String.concat ", " paths)
                files
            let targetFiles = quantFiles i
            let sourceFiles = quantFiles ii
            // the outputs and the source names are named by the file name
            let sameName =
                Array.append targetFiles sourceFiles
                |> Array.distinct
                |> Array.groupBy Path.GetFileNameWithoutExtension
                |> Array.filter (fun (_, paths) -> paths.Length > 1)
            if sameName.Length > 0 then
                failwithf "Input files from different directories have the same name: %s" (sameName |> Array.map (snd >> String.concat ", ") |> String.concat "; ")
            logger.Info (sprintf "%i target and %i source files, running on %i cores" targetFiles.Length sourceFiles.Length c)
            alignRuns dc o seed c targetFiles sourceFiles
            logger.Info "Done"
            0
        with
        | ex ->
            let ex = match ex with | :? AggregateException as agg when agg.InnerExceptions.Count = 1 -> agg.InnerExceptions.[0] | _ -> ex
            logger.Error (sprintf "QuantBasedAlignment_v2 failed: %s" ex.Message)
            logger.Trace (sprintf "%O" ex)
            NLog.LogManager.Flush()
            1
