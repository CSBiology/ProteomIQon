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
        let i     = results.GetResult InputFolder
        let fastA = results.GetResult FastA
        let gff3  = results.GetResult GFF3
        let o     = results.GetResult OutputDirectory
        let p     = results.GetResult ParamFile
        printfn "InputFilePath -i = %s" i
        printfn "InputFilePath -f = %s" fastA
        printfn "InputFilePath -g = %s" gff3
        printfn "InputFilePath -o = %s" o
        printfn "InputFilePath -p = %s" p
        Directory.CreateDirectory(o) |> ignore
        let subFolder = (i.Split ([|"\\"|], System.StringSplitOptions.None))
        Directory.CreateDirectory(o + "\\" + subFolder.[subFolder.Length - 1]) |> ignore
        let proteinInferenceParams = 
                Json.ReadAndDeserialize<Dto.ProteinInferenceParams> p
                |> Dto.ProteinInferenceParams.toDomain
        ProteinInference.inferProteins gff3 fastA proteinInferenceParams o i
        printfn "Hit any key to exit."
        System.Console.ReadKey() |> ignore
        0