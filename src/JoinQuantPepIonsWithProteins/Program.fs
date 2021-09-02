namespace ProteomIQon

open System
open System.IO
open System.Reflection
open Argu
open ProteomIQon.Core
open ProteomIQon.Core.InputPaths
open JoinQuantPepIonsWithProteins
open CLIArgumentParsing

module console1 =
    [<EntryPoint>]
    let main argv = 
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult QuantifiedPeptides |> List.map getPathRelativeToDir
        let ii = results.GetResult InferredProteins  |> List.map getPathRelativeToDir
        let o = results.GetResult OutputDirectory    |> getPathRelativeToDir
        let mf = results.Contains MatchFiles
        let c = 
            match results.TryGetResult Parallelism_Level with 
            | Some c    -> c
            | None      -> 1
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "JoinQuantPepIonsWithProteins"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let quantPepIonPaths = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.quant"))) i
            |> Array.ofSeq
        let protPaths = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.prot"))) ii
            |> Array.ofSeq
        let quantAndProtFiles =
            if mf then 
                quantPepIonPaths
                |> Array.choose (fun quantPepIonPath ->
                    match protPaths  |> Array.tryFind (fun protPath -> (Path.GetFileNameWithoutExtension protPath) = (Path.GetFileNameWithoutExtension quantPepIonPath)) with
                    | Some protPath -> Some(quantPepIonPath,protPath)
                    | None ->
                        logger.Trace (sprintf "no .prot File for %s" quantPepIonPath)
                        None
                )
            else
                Array.zip quantPepIonPaths protPaths
        quantAndProtFiles 
        |> FSharpAux.PSeq.withDegreeOfParallelism c
        |> FSharpAux.PSeq.iter (fun (quant,prot) -> joinQuantPepIonsWithProteins o quant prot)
        logger.Info "Done"
        0
