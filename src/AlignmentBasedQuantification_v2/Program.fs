namespace ProteomIQon

open System
open System.IO
open BioFSharp.Mz
open Core.InputPaths
open PeakSelection
open AlignmentBasedQuantification_v2
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
        let i  = results.GetResult InstrumentOutput    |> getPathRelativeToDir
        let ii = results.GetResult AlignedPeptides     |> getPathRelativeToDir
        let iv = results.GetResult QuantifiedPeptides  |> getPathRelativeToDir
        let o  = results.GetResult OutputDirectory     |> getPathRelativeToDir
        let p  = results.GetResult ParamFile           |> getPathRelativeToDir
        let d  = results.GetResult PeptideDataBase     |> getPathRelativeToDir
        let mf = results.Contains MatchFiles
        let dc = results.Contains DiagnosticCharts
        let c  = results.TryGetResult Parallelism_Level |> Option.defaultValue 1
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        extendLogging ()
        let logger = Logging.createLogger "AlignmentBasedQuantification_v2"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "AlignedPeptides -ii = %s" ii)
        logger.Info (sprintf "QuantifiedPeptides -iv = %s" iv)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        try
            let processParams =
                Json.ReadAndDeserialize<Dto.AlignmentBasedQuantificationV2Params> p
                |> Dto.AlignmentBasedQuantificationV2Params.toDomain
            if not (processParams.XicExtraction.ScanTimeWindow > 0. && processParams.MinScanTimeWindow > 0. && processParams.ScanTimeWindowSdFactor > 0.) then
                failwith "The parameters need XicExtraction.ScanTimeWindow, MinScanTimeWindow and ScanTimeWindowSdFactor > 0."
            if not (processParams.DecoyFraction >= 0. && processParams.DecoyFraction <= 1.) then
                failwith "DecoyFraction must be between 0 and 1."
            match processParams.XicExtraction.IonMobilityWindow with
            | Domain.Window.Fixed v when not (v > 0.) -> failwith "IonMobilityWindow Fixed needs a value greater than 0."
            | _ -> ()
            // warn before reading any run
            match processParams.XicExtraction.XicProcessing with
            | Domain.XicProcessing.Wavelet w when w.BorderPadMethod = FSharp.Stats.Signal.Padding.BorderPaddingMethod.Random ->
                logger.Warn "Wavelet BorderPadMethod Random pads the XIC borders with unseeded random values, so the results differ slightly between runs; use Zero for reproducible results."
            | _ -> ()
            if not (File.Exists d) then
                failwithf "The peptide data base %s does not exist." d
            for path in [i; ii; iv] do
                if not (File.Exists path || Directory.Exists path) then
                    failwithf "%s does not exist." path
            do
                use dbConnection = SearchDB.getDBConnection d
                BioFSharp.Mz.SearchDB.Db.SQLiteQuery.setIndexOnModSequenceAndGlobalMod dbConnection |> ignore
            let runs =
                if File.Exists i then
                    [|i, iv, ii|]
                elif Directory.Exists i && Directory.Exists ii && Directory.Exists iv then
                    let mzFiles = Directory.GetFiles(i, "*.mzlite") |> Array.sort
                    let alignFiles = Directory.GetFiles(ii, "*.align") |> Array.sort
                    let quantFiles = Directory.GetFiles(iv, "*.quant") |> Array.sort
                    if mf then
                        let byName (files: string []) (name: string) = files |> Array.tryFind (fun f -> Path.GetFileNameWithoutExtension f = name)
                        mzFiles
                        |> Array.choose (fun mzFile ->
                            let name = Path.GetFileNameWithoutExtension mzFile
                            match byName alignFiles name, byName quantFiles name with
                            | Some alignFile, Some quantFile -> Some (mzFile, quantFile, alignFile)
                            | alignFile, quantFile ->
                                logger.Warn (sprintf "%s is not quantified: %s" mzFile (if alignFile.IsNone then "no .align file" else "no .quant file"))
                                None
                            )
                    else
                        if mzFiles.Length <> alignFiles.Length || mzFiles.Length <> quantFiles.Length then
                            failwithf "%i .mzlite, %i .align and %i .quant files cannot be matched by position; use -mf to match them by name." mzFiles.Length alignFiles.Length quantFiles.Length
                        Array.zip3 mzFiles quantFiles alignFiles
                else
                    failwith "-i, -ii and -iv have to be either files or directories."
            logger.Info (sprintf "%i runs, running on %i cores" runs.Length c)
            let failed =
                runs
                |> FSharpAux.PSeq.map (fun (mzFile, quantFile, alignFile) ->
                    try
                        quantifyPeptides dc processParams o d mzFile quantFile alignFile
                        None
                    with
                    | ex ->
                        logger.Error (sprintf "%s failed: %s" mzFile ex.Message)
                        logger.Trace (sprintf "%O" ex)
                        Some mzFile
                    )
                |> FSharpAux.PSeq.withDegreeOfParallelism c
                |> Array.ofSeq
                |> Array.choose id
            if failed.Length > 0 then
                logger.Error (sprintf "%i of %i runs failed." failed.Length runs.Length)
                NLog.LogManager.Flush()
                1
            else
                logger.Info "Done"
                0
        with
        | ex ->
            let ex = match ex with | :? AggregateException as agg when agg.InnerExceptions.Count = 1 -> agg.InnerExceptions.[0] | _ -> ex
            logger.Error (sprintf "AlignmentBasedQuantification_v2 failed: %s" ex.Message)
            logger.Trace (sprintf "%O" ex)
            NLog.LogManager.Flush()
            1
