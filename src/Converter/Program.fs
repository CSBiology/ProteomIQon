namespace ProteomIQon

open System
open CLIArgumentParsing
open Argu
module console1 =

    [<EntryPoint>]
    let main argv = 
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let results = parser.Parse argv
        Library.printParams (results.GetAllResults())
        printfn "Hit any key to exit."
        System.Console.ReadKey() |> ignore
        0
