// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"

open ProteomIQon
open ProteomIQon.Dto

let defaultConsensusSpectralLibraryParams: Dto.ConsensusSpectralLibraryParams =
    {
        RTTolerance = 10.
        iRTPeptides = []
    }

let serialized = 
    defaultConsensusSpectralLibraryParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\SpectralLibraryParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\SpectralLibraryParams.json")
    |> Json.deserialize<Dto.ConsensusSpectralLibraryParams>
    |> ConsensusSpectralLibraryParams.toDomain