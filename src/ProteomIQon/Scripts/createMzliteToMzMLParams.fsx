
// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r @"../../../bin\MzMLToMzLite\net5.0\MzIO.dll"
#r @"../../../bin\MzMLToMzLite\net5.0\Newtonsoft.Json.dll"
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"

open ProteomIQon
open ProteomIQon.Dto
open MzIO

let defaultMzliteToMzMLParams :Dto.MzliteToMzMLParams = 

    {
        Compress                    = MzIO.Binary.BinaryDataCompressionType.NoCompression
        StartRetentionTime          = None
        EndRetentionTime            = None
    }


let serialized = 
    defaultMzliteToMzMLParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\MzliteToMzMLParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\MzliteToMzMLParams.json")
    |> Json.deserialize<Dto.MzMLtoMzLiteParams>
    |> MzMLtoMzLiteParams.toDomain
