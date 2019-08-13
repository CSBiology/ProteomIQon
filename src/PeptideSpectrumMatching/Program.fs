namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PeptideSpectrumMatching
module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv = 
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name)) 
        let usage  = parser.PrintUsage()
        printfn "%s" usage
        let results = parser.Parse argv
        let i = results.GetResult InstrumentOutput
        let o = results.GetResult OutputDirectory
        let p = results.GetResult ParamFile
        let d = results.GetResult PeptideDataBase
        printfn "InputFilePath -i = %s" i
        printfn "InputFilePath -o = %s" o
        printfn "InputFilePath -p = %s" p
        Directory.CreateDirectory(o) |> ignore
        let p = 
            Json.ReadAndDeserialize<Dto.PeptideSpectrumMatchingParams> p
            |> Dto.PeptideSpectrumMatchingParams.toDomain
        let dbConnection = 
            if File.Exists d then
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."

        if File.Exists i then
            printfn "singleFile"
            scoreSpectra p o dbConnection i
        elif Directory.Exists i then 
            printfn "multiple files"
            let files = 
                Directory.GetFiles(i,("*.mzlite"))                
            let c = 
                match results.TryGetResult Parallelism_Level with 
                | Some c    -> c
                | None      -> 1
            files 
            |> FSharpAux.PSeq.map (scoreSpectra p o dbConnection) 
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        else 
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."

        printfn "Hit any key to exit."
        System.Console.ReadKey() |> ignore
        0