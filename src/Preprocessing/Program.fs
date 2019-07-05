namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu

module console1 =

    [<EntryPoint>]
    let main argv = 
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name)) 
        let usage  = parser.PrintUsage()
        printfn "%s" usage
        let results = parser.Parse argv
        let manufacturerOutput = results.TryGetResult InstrumentOutput
        let outputDir = results.TryGetResult OutputDirectory
        let paramF = results.TryGetResult ParamFile
        match manufacturerOutput, outputDir, paramF with 
        | Some i, Some o , Some p -> 
            printfn "InputFilePath -i = %s" i
            printfn "InputFilePath -o = %s" o
            printfn "InputFilePath -p = %s" p
            Directory.CreateDirectory(o) |> ignore
            let processParams = 
                Json.ReadAndDeserialize<Dto.PreprocessingParams> p
                |> Dto.PreprocessingParams.toDomain
            Preprocessing.processFile processParams o i
        | _ -> failwith "params are not guut"
        System.Console.ReadKey() |> ignore
        printfn "Hit any key to exit."
        0
