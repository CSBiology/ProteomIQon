// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r @"../../../bin/preprocessing\net47\MzIO.dll"
#r @"../../../bin/preprocessing\net47\BioFSharp.Mz.dll"
#r @"../../../bin/preprocessing\net47\FSharpAux.IO.dll"
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain

let defaultMzMLtoMzLiteParams :Dto.MzMLtoMzLiteParams = 

    {
        Compress                    = MzIO.Binary.BinaryDataCompressionType.NoCompression
        StartRetentionTime          = None
        EndRetentionTime            = None 
        MS1PeakPicking              = PeakPicking.ProfilePeaks
        MS2PeakPicking              = PeakPicking.ProfilePeaks 
    }


let serialized = 
    defaultMzMLtoMzLiteParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\MzMLtoMzLiteParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\MzMLtoMzLiteParams.json")
    |> Json.deserialize<Dto.MzMLtoMzLiteParams>
    |> PreprocessingParams.toDomain
