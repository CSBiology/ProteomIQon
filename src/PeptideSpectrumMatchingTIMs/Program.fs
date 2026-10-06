namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
open PeptideSpectrumMatchingTIMs
open ProteomIQon.Core
open ProteomIQon.Core.InputPaths

module console1 =
    open BioFSharp.Mz

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
        let d = results.GetResult PeptideDataBase  |> getPathRelativeToDir
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "PeptideSpectrumMatchingTIMs"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Info (sprintf "Peptide data base -d = %s" d)
        logger.Trace (sprintf "CLIArguments: %A" results)
        use dbConnection =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the peptide data base is not a valid file path."
        let p = 
            Json.ReadAndDeserialize<Dto.PeptideSpectrumMatchingTIMsParams> p
            |> Dto.PeptideSpectrumMatchingTIMsParams.toDomain
        let files = 
            parsePaths MzIO.Reader.getMzLiteFiles i
            |> Array.ofSeq
        match files |> Array.tryFind (fun f -> not (String.Equals(Path.GetExtension f, ".mzlite", StringComparison.OrdinalIgnoreCase))) with
        | Some other -> failwithf "Only .mzlite files are supported, cannot search %s." other
        | None -> ()
        if files.Length = 0 then failwith "No .mzlite file found under the given input paths."
        let c =
            match results.TryGetResult Parallelism_Level with
            | Some c when c > 0 -> c
            | Some _ -> invalidArg "-c" "Must be positive."
            | None -> 1
        logger.Trace (sprintf "Program is running on %i cores" c)
        logger.Trace "Loading the peptide table."
        // 0 means: derive the ceiling from the heaviest precursor of every run
        let maxFragmentMass =
            match results.TryGetResult Max_Fragment_Mass with
            | Some m when Double.IsFinite m && m > 0. -> m
            | Some _ -> invalidArg "-f" "Must be finite and positive."
            | None -> 0.
        let slices =
            match results.TryGetResult Slices with
            | Some s when s > 0 -> s
            | Some _ -> invalidArg "-s" "Must be positive."
            | None -> 1
        let keepSpectra = results.Contains Keep_Spectra
        let mobilityMergePpm = results.TryGetResult Mobility_Merge_Ppm
        mobilityMergePpm |> Option.iter validateMobilityMergePpm
        validateParameters p maxFragmentMass slices
        validateOutputPaths o files
        let sdbParams, table = prepareTable p dbConnection (fun msg -> logger.Trace msg)
        logger.Trace "Loading the peptide table: finished."
        // as many files at once as -c says, each on one thread, against index slices built once
        // for the group; the next group starts when the group is done
        files
        |> Array.chunkBySize c
        |> Array.iter (fun group ->
            logger.Trace (sprintf "Scoring %i files together: %A" group.Length group)
            match mobilityMergePpm with
            | Some ppm -> scoreRunsWithMobilityTolerance ppm p o (fun msg -> logger.Trace msg) maxFragmentMass slices keepSpectra sdbParams table group
            | None -> scoreRuns p o (fun msg -> logger.Trace msg) maxFragmentMass slices keepSpectra sdbParams table group)
        logger.Info "Done"
        0
