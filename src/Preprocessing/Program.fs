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
                Core.MzLite.Reader.getMSFilePaths i
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

            //files
            //|> FSharpAux.PSeq.map (fun inp ->
            //                        try
            //                            Preprocessing.processFile processParams o inp
            //                        with
            //                        | ex ->
            //                            printfn "%A" ex

            //                      )
            //|> FSharpAux.PSeq.withDegreeOfParallelism c
            //|> Array.ofSeq
            //|> ignore
            with
            | ex -> printfn "%A" ex
        else
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        logger.Info "Hit any key to exit."
        System.Console.ReadKey() |> ignore
        0