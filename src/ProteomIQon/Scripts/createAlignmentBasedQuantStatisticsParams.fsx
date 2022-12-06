#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
#r @"../../../bin\alignmentbasedquantification\net5.0\FSharpAux.IO.dll"
#r "nuget: Newtonsoft.Json, 13.0.2"

open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain


let defaultAlignmentBasedQuantStatisticsParams :Dto.AlignmentBasedQuantStatisticsParams = 
    {
        PositiveQuantMzCutoff = 0.01
        PositiveQuantCutoff = 0.1
        NegativeQuantMzCutoff = 0.99
        NegativeQuantCutoff = 0.9
    }
 

let serialized = 
    defaultAlignmentBasedQuantStatisticsParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\AlignmentBasedQuantStatisticsParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\AlignmentBasedQuantStatisticsParams.json")
    |> Json.deserialize<Dto.QuantificationParams>
    |> QuantificationParams.toDomain




