namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
module console1 =

    [<EntryPoint>]
    let main argv = 
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let results = parser.Parse argv
        let i = results.GetResult InstrumentOutput
        let o = results.GetResult OutputDirectory
        let p = results.GetResult ParamFile
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "mzMLConverter"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let converterParams =
                Json.ReadAndDeserialize<Dto.MzMLConverterParams> p
                |> Dto.MzMLConverterParams.toDomain

        if File.Exists i || (Directory.Exists i && i.EndsWith(".d"))  then
            logger.Info "single file"
            MzMLConverter.convertFile converterParams o i
            logger.Trace (sprintf "Converting %s" i)
        elif Directory.Exists i then
            logger.Info "multiple files"
            let files =
                Core.MzIO.Reader.getMSFilePaths i
                |> Array.sort
            logger.Trace (sprintf "Converting multiple files: %A" files)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            try
            let partitionedFiles =
                files
                |> Array.splitInto c
            [for i in partitionedFiles do yield async { return i |> Array.map (MzMLConverter.convertFile converterParams o)}]
            |> Async.Parallel
            |> Async.RunSynchronously
            |> ignore
            with
            | ex -> printfn "%A" ex
        else
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        logger.Info "Done"
        0