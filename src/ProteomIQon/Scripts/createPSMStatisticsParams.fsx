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

let defaultPSMStatistics : Dto.PSMStatisticsParams = 
    {
        QValueThreshold             = 0.01
        PepValueThreshold           = 0.05
        ParseProteinIDRegexPattern  = (@"id")
        KeepTemporaryFiles          = true
    }
  
let serialized = 
    defaultPSMStatistics
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\pSMStatisticsParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\pSMStatisticsParams.json")
    |> Json.deserialize<Dto.PSMStatisticsParams>
    |> PSMStatisticsParams.toDomain
 
