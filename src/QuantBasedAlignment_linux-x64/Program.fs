namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
open QuantBasedAlignment
open System.Reflection
open ProteomIQon.Core.InputPaths
open QuantBasedAlignment

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv =
        printfn "Jou"
        QuantBasedAlignment.execute argv
        0
