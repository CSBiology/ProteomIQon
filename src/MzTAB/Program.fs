namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
open System.Reflection

module console1 =

    [<EntryPoint>]
    let main argv =
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let results = parser.Parse argv
        let i' = results.GetResult TableSortFile
        let ii' = results.GetResult ProtFile
        let iii' = results.GetResult QuantFile
        let iv' = results.GetResult QpsmFile
        let o' = results.GetResult OutputFile
        let p' = results.GetResult ParamFile
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let i = Path.Combine(directory, i')
        let ii = Path.Combine(directory, ii')
        let iii = Path.Combine(directory, iii')
        let iv = Path.Combine(directory, iv')
        let o = Path.Combine(directory, o')
        let p = Path.Combine(directory, p')
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
