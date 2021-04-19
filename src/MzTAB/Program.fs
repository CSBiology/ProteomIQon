namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
open System.Reflection
open ProteomIQon.Core.InputPaths

module console1 =

    [<EntryPoint>]
    let main argv =
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult TableSortFile |> getPathRelativeToDir
        let ii = results.GetResult ProtFile     |> getPathRelativeToDir
        let iii = results.GetResult QuantFile   |> getPathRelativeToDir
        let iv = results.GetResult QpsmFile     |> getPathRelativeToDir
        let o = results.GetResult OutputFile    |> getPathRelativeToDir
        let p = results.GetResult ParamFile     |> getPathRelativeToDir
        Logging.generateConfig (Path.GetDirectoryName o)
        if System.IO.File.Exists o then
            System.IO.File.Delete o
        let logger = Logging.createLogger "mzTab"
        logger.Info (sprintf "InputFile -i = %s" i)
        logger.Info (sprintf "InputFolder -ii = %s" ii)
        logger.Info (sprintf "InputFolder -iii = %s" iii)
        logger.Info (sprintf "InputFolder -iiii = %s" iv)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let processParams =
                Json.ReadAndDeserialize<Dto.MzTABParams> p
                |> Dto.MzTABParams.toDomain
        MzTAB.createMzTab o i ii iii iv processParams
        logger.Info "Done"
        0
