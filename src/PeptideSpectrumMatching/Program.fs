namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PeptideSpectrumMatching
open System.Reflection
open System
open ProteomIQon.Core.InputPaths

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv = 
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult InstrumentOutput |> List.map getPathRelativeToDir
        let o = results.GetResult OutputDirectory  |> getPathRelativeToDir
        let p = results.GetResult ParamFile        |> getPathRelativeToDir
        let d = results.GetResult PeptideDataBase  |> getPathRelativeToDir
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "PeptideSpectrumMatching"
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
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        let p = 
            Json.ReadAndDeserialize<Dto.PeptideSpectrumMatchingParams> p
            |> Dto.PeptideSpectrumMatchingParams.toDomain
        let files = 
            parsePaths "mzlite" i
            |> Array.ofSeq
        if files.Length = 1 then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Scoring spectra for %s" files.[0])
            scoreSpectra p o dbConnection files.[0]
        else
            logger.Info (sprintf "multiple files")
            logger.Trace (sprintf "Scoring multiple files: %A" files)
            let c = 
                match results.TryGetResult Parallelism_Level with 
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            files 
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> FSharpAux.PSeq.iter (scoreSpectra p o dbConnection)
        dbConnection.Dispose()
        logger.Info "Done"
        0