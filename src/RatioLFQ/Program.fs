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
        RatioLFQ.lfq i o pc rac rfc
        0