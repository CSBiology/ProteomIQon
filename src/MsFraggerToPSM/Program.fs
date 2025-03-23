namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
open System.Reflection
open ProteomIQon.Core.InputPaths
open ProteomIQon.Core
open ProteomIQon
open MsFraggerToPSM
open BioFSharp.Mz

module console1 =

    [<EntryPoint>]
    let main argv =
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult MSFraggerOutput  |> List.map getPathRelativeToDir
        let ii = results.GetResult MzLiteOutput     |> List.map getPathRelativeToDir
        let o = results.GetResult OutputDirectory  |> getPathRelativeToDir
        let d = results.GetResult PeptideDataBase        |> getPathRelativeToDir
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "MSFraggerToPSM"
        logger.Info (sprintf "MSFraggerOutput -i = %A" i)
        logger.Info (sprintf "MzLiteOutput -i = %A" ii)
        logger.Info (sprintf "OutputFolderPath -o = %s" o)
        logger.Info (sprintf "DBPath -p = %s" d)
        logger.Trace (sprintf "CLIArguments: %A" results)

        let dbConnection =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."

        let msFraggerValidatedFiles = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.tsv"))) i
            |> Array.ofSeq

        let mzLiteFiles = 
            parsePaths (fun path -> MzIO.Reader.getMzMLFiles path) i
            |> Array.ofSeq

        if mzLiteFiles.Length = 1  then
            logger.Info "single file"
            logger.Trace (sprintf "processing %s" mzLiteFiles.[0])
            convertToPSM msFraggerValidatedFiles.[0] mzLiteFiles.[0] o dbConnection
        else
            logger.Info "multiple files"
            logger.Trace (sprintf "processing multiple files: %A" mzLiteFiles)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            try
            let partitionedFiles =
                Array.zip msFraggerValidatedFiles mzLiteFiles
                |> Array.splitInto c
            [for i in partitionedFiles do yield async { return i |> Array.map (fun (msFragger, mzLite) -> (convertToPSM msFragger mzLite o dbConnection))}]
            |> Async.Parallel
            |> Async.RunSynchronously
            |> ignore
            with
            | ex -> printfn "%A" ex
        logger.Info "Done"
        0