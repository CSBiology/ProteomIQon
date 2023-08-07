#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
#r @"../../../bin\alignmentbasedquantification\net5.0\FSharpAux.IO.dll"
#r "nuget: Newtonsoft.Json, 13.0.2"

open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain


let defaultLabelEfficiencyCalculationParams :Dto.LabelEfficiencyCalculationParams = 
    {
        LowerBound = 0.8
        UpperBound = 1.
    }
 

let serialized = 
    defaultLabelEfficiencyCalculationParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\LabelEfficiencyCalculationParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\LabelEfficiencyCalculationParams.json")
    |> Json.deserialize<Dto.LabelEfficiencyCalculationParams>
    |> LabelEfficiencyCalculationParams.toDomain




