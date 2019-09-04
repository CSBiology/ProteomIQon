namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PSMBasedQuantification

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
        let ii = results.GetResult ScoredPSMs
        let o = results.GetResult OutputDirectory
        let p = results.GetResult ParamFile
        let d = results.GetResult PeptideDataBase
        printfn "InputFilePath -i = %s" i
        printfn "InputFilePath -ii = %s" ii
        printfn "InputFilePath -o = %s" o
        printfn "InputFilePath -p = %s" p
        Directory.CreateDirectory(o) |> ignore
        let p = 
            Json.ReadAndDeserialize<Dto.QuantificationParams> p
            |> Dto.QuantificationParams.toDomain
        
        let dbConnection = 
            if File.Exists d then
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        printfn "Set Index on data base if not present."
        setIndexOnModSequenceAndGlobalMod dbConnection |> ignore
        printfn "Set Index on data base if not present: finished"

        if File.Exists i then
            printfn "singleFile"
            quantifyPeptides p o dbConnection i ii
        elif Directory.Exists i && Directory.Exists ii then 
            printfn "multiple files"
            let mzfiles = 
                Directory.GetFiles(i,("*.mzlite"))                
            let pepfiles = 
                Directory.GetFiles(ii,("*.qpsm"))                
            let mzFilesAndPepFiles = 
                mzfiles
                |> Array.choose (fun mzFilePath -> 
                                match pepfiles  |> Array.tryFind (fun pepFilePath -> (Path.GetFileNameWithoutExtension pepFilePath) = (Path.GetFileNameWithoutExtension mzFilePath)) with
                                | Some pepFile -> Some(mzFilePath,pepFile)
                                | None -> 
                                    printfn "no qpsmFileFor %s" mzFilePath
                                    None

                              )
            if mzfiles.Length <> mzFilesAndPepFiles.Length then 
                printfn "There are %i mzFiles but %i files containing scored PSMs" mzfiles.Length pepfiles.Length
            let c = 
                match results.TryGetResult Parallelism_Level with 
                | Some c    -> c
                | None      -> 1
            mzFilesAndPepFiles 
            |> FSharpAux.PSeq.map (fun (i,ii) -> quantifyPeptides p o dbConnection i ii) 
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        else 
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."

        printfn "Hit any key to exit."
        System.Console.ReadKey() |> ignore
        0