
// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r @"../../../bin\mzmlconverter\net5.0\MzIO.dll"
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"

open ProteomIQon
open ProteomIQon.Dto
open MzIO

let defaultMzMLConverterParams :Dto.MzMLConverterParams = 

    {
        Compress                    = MzIO.Binary.BinaryDataCompressionType.NoCompression
        StartRetentionTime          = None
        EndRetentionTime            = None
    }


let serialized = 
    defaultMzMLConverterParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\mzMLConverterParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\mzMLConverterParams.json")
    |> Json.deserialize<Dto.MzMLConverterParams>
    |> MzMLConverterParams.toDomain
