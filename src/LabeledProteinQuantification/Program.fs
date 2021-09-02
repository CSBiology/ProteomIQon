namespace ProteomIQon

open System
open System.IO
open System.Reflection
open Argu
open ProteomIQon.Core
open ProteomIQon.Core.InputPaths
open CLIArgumentParsing

module console1 =

    [<EntryPoint>]
    let main argv = 
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult ProteinAssignedQuantPepIons |> List.map getPathRelativeToDir
        let o = results.GetResult OutputDirectory  |> getPathRelativeToDir
        let p = results.GetResult ParamFile        |> getPathRelativeToDir
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "LabeledProteinQuantification"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let p = 
            Json.ReadAndDeserialize<Dto.LabeledQuantificationParams> p
            |> Dto.LabeledQuantificationParams.toDomain
        let files = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.quantAndProt"))) i
            |> Array.ofSeq
        LabeledProteinQuantification.labeledQuantification p o files
        logger.Info "Done"
        0
