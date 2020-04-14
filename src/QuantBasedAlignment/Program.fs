namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open QuantBasedAlignment

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv =
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let results = parser.Parse argv
        let i = results.GetResult QuantifiedPeptides
        let o = results.GetResult OutputDirectory
        let p = results.GetResult ParamFile
        Logging.generateConfig o
        let logger = Logging.createLogger "PSMBasedQuantification"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let p =
            Json.ReadAndDeserialize<Dto.QuantificationParams> p
            |> Dto.QuantificationParams.toDomain
        if File.Exists i then
            logger.Info "single file detected"
            failwithf "%s is a file path, please specify a directory containing .quant files" i
        elif Directory.Exists i then
            logger.Info "directory found"
            let quantfiles =
                Directory.GetFiles(i,("*.quant"))
            logger.Trace (sprintf "PEP files : %A" quantfiles)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            printParams 
            |> ignore
        else
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        logger.Info "Done"
        0
