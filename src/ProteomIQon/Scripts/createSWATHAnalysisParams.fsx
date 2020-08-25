// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"

open ProteomIQon
open ProteomIQon.Dto

let defaultSWATHAnalysisParams: Dto.SWATHAnalysisParams =
    {
        PeptideList          = None
        MatchingTolerancePPM = 100.
    }

let serialized = 
    defaultSWATHAnalysisParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\SWATHAnalysisParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\SWATHAnalysisParams.json")
    |> Json.deserialize<Dto.SWATHAnalysisParams>
    |> SWATHAnalysisParams.toDomain