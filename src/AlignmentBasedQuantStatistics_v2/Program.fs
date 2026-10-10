namespace ProteomIQon

open System
open System.IO
open Core.InputPaths
open AlignmentBasedQuantStatistics_v2
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

    /// Matches the quant, align and aligned quant files by file name (matchFiles) or by position.
    let matchFiles (logger: NLog.Logger) (matchByName: bool) (quantFiles: string []) (alignFiles: string []) (alignQuantFiles: string []) =
        if matchByName then
            let byName (files: string []) (path: string) =
                files |> Array.tryFind (fun f -> Path.GetFileNameWithoutExtension f = Path.GetFileNameWithoutExtension path)
            quantFiles
            |> Array.choose (fun quantFile ->
                match byName alignFiles quantFile, byName alignQuantFiles quantFile with
                | Some alignFile, Some alignQuantFile -> Some (quantFile, alignFile, alignQuantFile)
                | None, _ ->
                    logger.Warn (sprintf "No align file for %s, the run is not used." quantFile)
                    None
                | _, None ->
                    logger.Warn (sprintf "No aligned quant file for %s, the run is not used." quantFile)
                    None
                )
        elif quantFiles.Length <> alignFiles.Length || quantFiles.Length <> alignQuantFiles.Length then
            failwithf "%i quant, %i align and %i aligned quant files cannot be matched by position; use -mf to match them by file name." quantFiles.Length alignFiles.Length alignQuantFiles.Length
        else
            Array.zip3 quantFiles alignFiles alignQuantFiles

    [<EntryPoint>]
    let main argv =
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i   = results.GetResult Quant |> List.map getPathRelativeToDir
        let ii  = results.GetResult Align |> List.map getPathRelativeToDir
        let iii = results.GetResult AlignedQuant |> List.map getPathRelativeToDir
        let l   = results.GetResult QuantLearn |> List.map getPathRelativeToDir
        let ll  = results.GetResult AlignLearn |> List.map getPathRelativeToDir
        let lll = results.GetResult AlignedQuantLearn |> List.map getPathRelativeToDir
        let outputDir = results.GetResult OutputDirectory |> getPathRelativeToDir
        let paramPath = results.GetResult ParamFile |> getPathRelativeToDir
        let dc  = results.Contains DiagnosticCharts
        let ts  = results.Contains TestSet
        let mf  = results.Contains MatchFiles
        let c   = results.TryGetResult Parallelism_Level |> Option.defaultValue 1
        Directory.CreateDirectory(outputDir) |> ignore
        Logging.generateConfig outputDir
        extendLogging ()
        let logger = Logging.createLogger "AlignmentBasedQuantStatistics_v2"
        logger.Info (sprintf "QuantFilePath -i = %A" i)
        logger.Info (sprintf "AlignFilePath -ii = %A" ii)
        logger.Info (sprintf "AlignedQuantFilePath -iii = %A" iii)
        logger.Info (sprintf "QuantFilePathForLearning -l = %A" l)
        logger.Info (sprintf "AlignFilePathForLearning -ll = %A" ll)
        logger.Info (sprintf "AlignedQuantFilePathForLearning -lll = %A" lll)
        logger.Info (sprintf "OutputFilePath -o = %s" outputDir)
        logger.Info (sprintf "ParamFilePath -p = %s" paramPath)
        logger.Trace (sprintf "CLIArguments: %A" results)
        try
            let p =
                Json.ReadAndDeserialize<Dto.AlignmentBasedQuantStatisticsV2Params> paramPath
                |> Dto.AlignmentBasedQuantStatisticsV2Params.toDomain
            logger.Info (sprintf "Parameters: %A" p)
            if not (p.PositiveRtSdFactor > 0. && p.NegativeRtSdFactor >= p.PositiveRtSdFactor && p.MaxAbsLog2QuantRatio > 0.) then
                failwith "The parameters need 0 < PositiveRtSdFactor <= NegativeRtSdFactor and MaxAbsLog2QuantRatio > 0."
            if p.NumberOfFolds < 2 then
                failwith "NumberOfFolds must be at least 2."
            let files (extension: string) (paths: string list) =
                for path in paths do
                    if not (File.Exists path || Directory.Exists path) then
                        failwithf "%s does not exist." path
                parsePaths (fun path -> Directory.GetFiles(path, "*" + extension) |> Array.sort) paths
                |> Array.ofSeq
                |> Array.map Path.GetFullPath
            let scoringRuns = matchFiles logger mf (files ".quant" i) (files ".align" ii) (files ".quant" iii)
            let learningRuns = matchFiles logger mf (files ".quant" l) (files ".align" ll) (files ".quant" lll)
            // the outputs, test sets and logs are named by the file name of the PSM .quant file
            for runs in [ scoringRuns; learningRuns ] do
                let sameName =
                    runs
                    |> Array.distinct
                    |> Array.groupBy (fun (quantPath, _, _) -> Path.GetFileNameWithoutExtension quantPath)
                    |> Array.filter (fun (_, files) -> files.Length > 1)
                if sameName.Length > 0 then
                    failwithf "Input files from different directories have the same name: %s" (sameName |> Array.map (fun (_, files) -> files |> Array.map (fun (quantPath, _, _) -> quantPath) |> String.concat ", ") |> String.concat "; ")
            if scoringRuns.Length = 0 then failwith "No runs to score."
            if learningRuns.Length = 0 then failwith "No runs to learn from."
            logger.Info (sprintf "%i runs to score, %i runs to learn from, reading on %i cores" scoringRuns.Length learningRuns.Length c)
            scoreRuns dc ts outputDir c p scoringRuns learningRuns
            logger.Info "Done"
            0
        with
        | ex ->
            let ex = match ex with | :? AggregateException as agg when agg.InnerExceptions.Count = 1 -> agg.InnerExceptions.[0] | _ -> ex
            logger.Error (sprintf "AlignmentBasedQuantStatistics_v2 failed: %s" ex.Message)
            logger.Trace (sprintf "%O" ex)
            NLog.LogManager.Flush()
            1
