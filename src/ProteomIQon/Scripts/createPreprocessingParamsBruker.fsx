// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r @"../../../bin/preprocessing\net47\BioFSharp.Mz.dll"
#r @"../../../bin/preprocessing\net47\FSharpAux.IO.dll"
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain

let defaultPreprocessingParams :Dto.PreprocessingParams = 

         
    {
        Compress                    = Compression.NoCompression
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
