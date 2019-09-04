// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../packages\BioFSharp\lib\netstandard2.0\BioFSharp.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"
#r @"../../../packages\BioFSharp.Mz\lib\netstandard2.0\BioFSharp.Mz.dll"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"

open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain
open BioFSharp.Mz

let defaultProteinInferenceParams: Dto.ProteinInferenceParams =
    {
        ProteinIdentifierRegex = "Cre\\S+"
        Protein                = ProteinInference.IntegrationStrictness.Maximal
        Peptide                = ProteinInference.PeptideUsageForQuantification.Minimal
        GroupFiles             = true
    }

let serialized = 
    defaultProteinInferenceParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\ProteinInferenceParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\ProteinInferenceParams.json")
    |> Json.deserialize<Dto.ProteinInferenceParams>
    |> ProteinInferenceParams.toDomain

