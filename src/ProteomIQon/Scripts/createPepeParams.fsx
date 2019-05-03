// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.dll"

open ProteomIQon
open ProteomIQon.Domain
open ProteomIQon.Dto
//open ProteomIQon.Library

let defaultPepeParams = 
    {
        QValueThreshold     = 0.01
        PepValueThreshold   = 0.05
        ProteinIDRegex      = (@"\w*")
    }
  
let serialized = 
    defaultPepeParams
    |> Json.serialize

System.IO.File.WriteAllText( @"C:\Users\david\Source\Repos\netCoreRepos\ProteomIQon\src\ProteomIQon\defaultParams\pepeParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(@"C:\Users\david\Source\Repos\netCoreRepos\ProteomIQon\src\ProteomIQon\defaultParams\pepeParams.json")
    |> Json.deserialize<PEPEParams>
    |> PEPEParams.toDomain

