namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open SWATHAnalysis
module console1 =

    [<EntryPoint>]
    let main argv = 
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let results = parser.Parse argv
        let i = results.GetResult InstrumentOutput 
        let o = results.GetResult OutputDirectory
        let p = results.GetResult ParamFile
        let l = results.GetResult Library
        Logging.generateConfig o
        let logger = Logging.createLogger "SWATHAnalysis"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Info (sprintf "LibraryFilePath -l = %s" l)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let p =
            Json.ReadAndDeserialize<Dto.SWATHAnalysisParams> p
            |> Dto.SWATHAnalysisParams.toDomain
        if File.Exists i then
            logger.Info "single file"
            quantify l p i o
        elif Directory.Exists i then
            logger.Info "multiple files"
            let files = 
                Directory.GetFiles(i,("*.mzlite"))
            logger.Trace (sprintf "Quantifying multiple swath files: %A" files)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            try
                let partitionedFiles =
                    files
                    |> Array.splitInto c
                [for i in partitionedFiles do yield async { return i |> Array.map (fun file -> (quantify l p file o))}]
                |> Async.Parallel
                |> Async.RunSynchronously
                |> ignore
            with | ex -> printfn "%A" ex
        else
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        printfn "Hit any key to exit."
        System.Console.ReadKey() |> ignore
        0
