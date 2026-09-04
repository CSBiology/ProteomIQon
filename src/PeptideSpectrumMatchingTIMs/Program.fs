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
        let dbConnection =
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
            | _ -> 1
        logger.Trace (sprintf "Program is running on %i cores" c)
        logger.Trace "Building the fragment index."
        let sdbParams, table, index = prepareIndex p dbConnection c (fun msg -> logger.Trace msg)
        logger.Trace "Building the fragment index: finished."
        if files.Length = 1 then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Scoring spectra for %s" files.[0])
            scoreSpectra p o sdbParams table index files.[0]
        else
            logger.Info (sprintf "multiple files")
            logger.Trace (sprintf "Scoring multiple files: %A" files)
            files 
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> FSharpAux.PSeq.iter (scoreSpectra p o sdbParams table index)
        dbConnection.Dispose()
        logger.Info "Done"
        0
