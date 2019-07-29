// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
#r @"../../../packages\BioFSharp.Mz\lib\netstandard2.0\BioFSharp.Mz.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"

open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain

let defaultPreprocessingParams :Dto.PreprocessingParams = 

         
    {
        Compress                    = true
        StartRetentionTime          = None
        EndRetentionTime            = None 
        MS1PeakPicking              = PeakPicking.Centroid CentroidizationMode.Manufacturer
        MS2PeakPicking              = PeakPicking.Centroid CentroidizationMode.Manufacturer
    }


let serialized = 
    defaultPreprocessingParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\preprocessingParams_Bruker.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\preprocessingParams_Bruker.json")
    |> Json.deserialize<Dto.PreprocessingParams>
    |> PreprocessingParams.toDomain
