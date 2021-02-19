namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
module console1 =

    [<EntryPoint>]
    let main argv =
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let results = parser.Parse argv
        let iTab = results.GetResult TableSortFile
        let iProt = results.GetResult ProtFile
        let iQuant = results.GetResult QuantFile
        let iQpsm = results.GetResult QpsmFile
        let o = results.GetResult OutputFile
        let p = results.GetResult ParamFile
        Logging.generateConfig (Path.GetDirectoryName o)
        if System.IO.File.Exists o then
            System.IO.File.Delete o
        let logger = Logging.createLogger "mzTab"
        logger.Info (sprintf "InputFile -i = %s" iTab)
        logger.Info (sprintf "InputFolder -ii = %s" iProt)
        logger.Info (sprintf "InputFolder -iii = %s" iQuant)
        logger.Info (sprintf "InputFolder -iiii = %s" iQpsm)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let processParams =
                Json.ReadAndDeserialize<Dto.MzTABParams> p
                |> Dto.MzTABParams.toDomain
        MzTAB.createMzTab o iTab iProt iQuant iQpsm processParams
        logger.Info "Done"
        0
