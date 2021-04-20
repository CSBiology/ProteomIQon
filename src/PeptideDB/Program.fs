namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open System.Reflection
open ProteomIQon.Core.InputPaths

module console1 =

    [<EntryPoint>]
    let main argv =
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult FastaPath |> getPathRelativeToDir
        let o = results.GetResult OutputDirectory |> getPathRelativeToDir
        let p = results.GetResult ParamFile |> getPathRelativeToDir
        Logging.generateConfig o
        let logger = Logging.createLogger "PeptideDB"
        logger.Trace (sprintf "CLIArguments %A" results)
        logger.Info (sprintf "OutputFilePath %s" o)
        logger.Info (sprintf "ParamFilePath %s" p)
        Directory.CreateDirectory(o) |> ignore
        let processParams =
            Json.ReadAndDeserialize<Dto.PeptideDBParams> p
            |> Dto.PeptideDBParams.toDomain           
        PeptideDB.createPeptideDB processParams o i
        0