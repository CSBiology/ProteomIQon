#r "netstandard"
#r @"../../../bin\MzMLToMzLite\net5.0\Newtonsoft.Json.dll"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"

open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Dto.RatioLFQParams
open ProteomIQon.Domain

let defaultRatioLFQParams : Dto.RatioLFQParams =
    {
        Reference     = ReferenceChannel.Heavy
        MinSharedIons = 1
        CentreRuns    = false
        PairWeighting = PairWeighting.SharedIonCount
    }

let serialized =
    defaultRatioLFQParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\RatioLFQParams.json",serialized)

let deserialized =
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\RatioLFQParams.json")
    |> Json.deserialize<Dto.RatioLFQParams>
    |> RatioLFQParams.toDomain
