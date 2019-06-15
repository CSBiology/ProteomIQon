namespace ProteomIQon

open System
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
        let outputDir = results.TryGetResult OutputDirectory
        let paramF = results.TryGetResult ParamFile
        match outputDir, paramF with 
        | Some o , Some p -> 
            printfn "Output directory -o = %s" o
            printfn "Parameter file path -p = %s" p
            let processParams = 
                Json.ReadAndDeserialize<Dto.PeptideDBParams> p
                |> Dto.PeptideDBParams.toDomain
            PeptideDB.createPeptideDB processParams o 
        | _ -> failwith "params are not guut"
        System.Console.ReadKey() |> ignore
        printfn "Hit any key to exit."
        0
