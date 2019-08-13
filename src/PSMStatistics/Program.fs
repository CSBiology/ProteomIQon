namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PSMStatistics

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv = 
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name)) 
        let usage  = parser.PrintUsage()
        printfn "%s" usage
        let results = parser.Parse argv
        let i = results.GetResult PSMs
        let o = results.GetResult OutputDirectory
        let p = results.GetResult ParamFile
        let d = results.GetResult PeptideDataBase
        printfn "InputFilePath -i = %s" i
        printfn "InputFilePath -o = %s" o
        printfn "InputFilePath -p = %s" p
        printfn "InputFilePath -p = %s" d
        Directory.CreateDirectory(o) |> ignore
        let p = 
            Json.ReadAndDeserialize<Dto.PSMStatisticsParams> p
            |> Dto.PSMStatisticsParams.toDomain
        let dbConnection = 
            if File.Exists d then
                SearchDB.getDBConnection d
            else
                failwith "The given path to the peptide database is neither a valid file path nor a valid directory path."

        if File.Exists i then
            printfn "singleFile"
            pepValueCalcAndProteinInference p o dbConnection i
        elif Directory.Exists i then 
            printfn "multiple files"
            let files = 
                Directory.GetFiles(i,("*.psm"))                
            let c = 
                match results.TryGetResult Parallelism_Level with 
                | Some c    -> c
                | None      -> 1
            files 
            |> FSharpAux.PSeq.map (pepValueCalcAndProteinInference p o dbConnection) 
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        else 
            failwith "The given path to the PSMs is neither a valid file path nor a valid directory path."

        printfn "Hit any key to exit."
        System.Console.ReadKey() |> ignore
        0