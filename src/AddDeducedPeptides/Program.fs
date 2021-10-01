namespace ProteomIQon

open System
open System.IO
open System.Reflection
open Argu
open ProteomIQon.Core
open ProteomIQon.Core.InputPaths
open AddDeducedPeptides
open CLIArgumentParsing

module console1 =

    [<EntryPoint>]
    let main argv = 
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i  = results.GetResult QuantFiles |> List.map getPathRelativeToDir
        let ii = results.GetResult ProtFiles |> List.map getPathRelativeToDir
        let o  = results.GetResult OutputDirectory |> getPathRelativeToDir
        Logging.generateConfig o
        let logger = Logging.createLogger "AddDeducedPeptides"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "InputFilePath -ii = %A" ii)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let quantFiles = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.quant"))) i
            |> Array.ofSeq
        let protFiles = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.prot"))) ii
            |> Array.ofSeq
        logger.Trace (sprintf "Inputfiles: \n%A\n%A" quantFiles protFiles)
        addDeducedPeptides quantFiles protFiles o
        logger.Info "Done"
        0
