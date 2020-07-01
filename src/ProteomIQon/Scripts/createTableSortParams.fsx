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

let defaultTableSortParams: Dto.TableSortParams =
    {
        SeparatorIn                 = "\t"
        SeparatorOut                = '\t'
        QuantFieldsToFilterOn       = [|(FilterOnField.create "Quant_Light" (None) (Some 0.)); (FilterOnField.create "Quant_Heavy" (None) (Some 0.))|]
        ProtFieldsToFilterOn        = [|(*(FilterOnField.create "EvidenceClass" None (Some 1.))*)|]
        QuantColumnsOfInterest      = [|"Quant_Light";"Quant_Heavy"|]
        ProtColumnsOfInterest       = [|"DistinctPeptideCount"; "QValue"; "TargetScore"|]
        AggregatorFunction          = AggregationMethod.Median
        AggregatorFunctionIntensity = AggregationMethod.Median
        AggregatorPepToProt         = AggregationMethod.Median
        TukeyQuant                  = [|("Quant_Light", 2.); ("Quant_Heavy", 2.)|]
        TukeyProt                   = [||]
        Labeled                     = true
    }

let serialized = 
    defaultTableSortParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\TableSortParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\TableSortParams.json")
    |> Json.deserialize<Dto.TableSortParams>
    |> TableSortParams.toDomain
