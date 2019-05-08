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
        let input = 
            if results.Contains InstrumentOutput then 
                results.GetResults InstrumentOutput
            else [""]
        Library.printParams (results.GetAllResults())
        printfn "Hit any key to exit."
        System.Console.ReadKey() |> ignore
        0
