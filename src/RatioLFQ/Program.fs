namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
open System.Reflection
open ProteomIQon.Core.InputPaths
open ProteomIQon.Core
open ProteomIQon

module console1 =

    [<EntryPoint>]
    let main argv =
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult InputFile   |> getPathRelativeToDir
        let o = results.GetResult OutputFile  |> getPathRelativeToDir
        let pc = results.GetResult ProteinColumn
        let rac = results.GetResult RatioColumnEnding
        let rfc = results.GetResult ReferenceColumnEnding
        let resultDirectory = 
            (Directory.GetParent o).FullName
            |> getPathRelativeToDir
        Directory.CreateDirectory(resultDirectory) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "MzMLToMzLite"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Trace (sprintf "CLIArguments: %A" results)
        RatioLFQ.lfq i o pc rac rfc
        logger.Info "Done"
        0