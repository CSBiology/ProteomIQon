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
        // light = Name of the unlabeled peptide quantification column in the quant file
        // heavy = Name of the labeled peptide quantification column in the quant file (optional, if present, "Ratio" column is automatically added)
        // proteinIDs = Name of the column containing the protein names in the prot file
        // pepSequence = Name of the column containing the peptide sequence in the quant file
        // pepSequences = Name of the column containing the peptide sequences in the prot file
        EssentialFields             = EssentialFields.create "Quant_Light" (Some "Quant_Heavy") "GroupOfProteinIDs" "StringSequence" "PeptideSequence"
        QuantFieldsToFilterOn       = [|(FilterOnField.create "Quant_Light" (None) (Some 0.)); (FilterOnField.create "Quant_Heavy" (None) (Some 0.))|]
        ProtFieldsToFilterOn        = [|(*(FilterOnField.create "Class" None (Some 1.))*)|]
        // Columns of interest can be any column present in the original tables. The column "Ratio" is automatically added
        // when heavy labeled peptides are present. "DistinctPeptideCount" is calculated from the amount of different peptides used
        // for the ratio calculation and can be included.
        QuantColumnsOfInterest      = [|"Quant_Light";"Quant_Heavy"|]
        ProtColumnsOfInterest       = [|"DistinctPeptideCount"; "QValue"|]
        StatisticalMeasurements     = [|("Ratio", StatisticalMeasurement.CV); ("Ratio", StatisticalMeasurement.SEM); ("Ratio", StatisticalMeasurement.StDev)|]
        AggregatorFunction          = AggregationMethod.Median
        AggregatorFunctionIntensity = AggregationMethod.Median
        AggregatorPepToProt         = AggregationMethod.Median
        Tukey                       = [|("Ratio",2., Transform.Log2)|]
    }

let serialized = 
    defaultTableSortParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\TableSortParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\TableSortParams.json")
    |> Json.deserialize<Dto.TableSortParams>
    |> TableSortParams.toDomain
