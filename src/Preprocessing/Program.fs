namespace ProteomIQon

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
        let i' = results.GetResult InstrumentOutput
        let o' = results.GetResult OutputDirectory
        let p' = results.GetResult ParamFile
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let i = Path.Combine(directory, i')
        let o = Path.Combine(directory, o')
        let p = Path.Combine(directory, p')
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "Preprocessing"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let processParams =
                Json.ReadAndDeserialize<Dto.PreprocessingParams> p
                |> Dto.PreprocessingParams.toDomain

        if File.Exists i || (Directory.Exists i && i.EndsWith(".d"))  then
            logger.Info "single file"
            Preprocessing.processFile processParams o i
            logger.Trace (sprintf "Preprocessing %s" i)
        elif Directory.Exists i then
            logger.Info "multiple files"
            let files =
                Core.MzIO.Reader.getMSFilePaths i
                |> Array.sort
            logger.Trace (sprintf "Preprocessing multiple files: %A" files)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            try
            let partitionedFiles =
                files
                |> Array.splitInto c
            [for i in partitionedFiles do yield async { return i |> Array.map (Preprocessing.processFile processParams o)}]
            |> Async.Parallel
            |> Async.RunSynchronously
            |> ignore
            with
            | ex -> printfn "%A" ex
        else
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        logger.Info "Done"
        0