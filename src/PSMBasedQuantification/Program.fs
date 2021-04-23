namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PSMBasedQuantification
open System.Reflection
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
        let ii = results.GetResult ScoredPSMs      |> List.map getPathRelativeToDir
        let p = results.GetResult ParamFile        |> getPathRelativeToDir
        let o = results.GetResult OutputDirectory  |> getPathRelativeToDir
        let d = results.GetResult PeptideDataBase  |> getPathRelativeToDir
        Logging.generateConfig o
        let logger = Logging.createLogger "PSMBasedQuantification"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "ScoredPSMsPath -ii = %A" ii)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let p =
            Json.ReadAndDeserialize<Dto.QuantificationParams> p
            |> Dto.QuantificationParams.toDomain
        let dbConnection =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        logger.Trace "Set Index on data base if not present."
        SearchDB'.setIndexOnModSequenceAndGlobalMod dbConnection |> ignore
        logger.Trace "Set Index on data base if not present: finished"
        let mzfiles = 
            parsePaths "mzlite" i
            |> Array.ofSeq
        let pepfiles = 
            parsePaths "qpsm" ii
            |> Array.ofSeq
        if mzfiles.Length = 1 && pepfiles.Length = 1 then
            logger.Info "single file"
            logger.Trace (sprintf "mz files : %A" mzfiles)
            logger.Trace (sprintf "PEP files : %A" pepfiles)
            quantifyPeptides p o dbConnection mzfiles.[0] pepfiles.[0]
        else
            logger.Info "multiple files"
            logger.Trace (sprintf "mz files : %A" mzfiles)
            logger.Trace (sprintf "PEP files : %A" pepfiles)
            let mzFilesAndPepFiles =
                mzfiles
                |> Array.choose (fun mzFilePath ->
                    match pepfiles  |> Array.tryFind (fun pepFilePath -> (Path.GetFileNameWithoutExtension pepFilePath) = (Path.GetFileNameWithoutExtension mzFilePath)) with
                    | Some pepFile -> Some(mzFilePath,pepFile)
                    | None ->
                        logger.Trace (sprintf "no qpsmFileFor %s" mzFilePath)
                        None
                )
            if mzfiles.Length <> mzFilesAndPepFiles.Length then
                logger.Info (sprintf "There are %i mzFiles but %i files containing scored PSMs" mzfiles.Length pepfiles.Length)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            mzFilesAndPepFiles
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> FSharpAux.PSeq.iter (fun (i,ii) -> quantifyPeptides p o dbConnection i ii)
        dbConnection.Dispose()

        logger.Info "Done"
        0