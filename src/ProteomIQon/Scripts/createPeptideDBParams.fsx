// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\net47\BioFSharp.dll"
#r @"../../../bin\ProteomIQon\net47\BioFSharp.Mz.dll"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
#r @"../../../packages\FSharpAux.IO\lib\net47\FSharpAux.IO.dll"
#r @"../../../packages\FSharpAux.IO\lib\net47\FSharpAux.dll"

open ProteomIQon
open ProteomIQon.Domain
open ProteomIQon.Dto
//open ProteomIQon.Library
open BioFSharp.Mz.SearchDB



let defaultPeptideDBParams: Dto.PeptideDBParams = 
        {
        Name                        = "AraTest"
        FastaPath                   = @"D:\David\N15-test\db\Arabidopsis_thaliana.TAIR10.pep.all.fasta"
        ParseProteinIDRegexPattern  = "id"
        Protease                    = Protease.Trypsin
        MinMissedCleavages          = 0
        MaxMissedCleavages          = 2
        MaxMass                     = 15000.0
        MinPepLength                = 4
        MaxPepLength                = 65
        IsotopicMod                 = [IsotopicMod.N15]
        MassMode                    = MassMode.Monoisotopic
        FixedMods                   = []
        VariableMods                = [Modification.Oxidation'Met';Modification.Acetylation'ProtNTerm';Modification.Carbamidomethyl'Cys';Modification.Pyro_Glu'GlnNterm';Modification.Pyro_Glu'GluNterm']
        VarModThreshold             = 4
        }
  
let serialized = 
    defaultPeptideDBParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\peptideDBParamsGiada.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\peptideDBParams.json")
    |> Json.deserialize<Dto.PeptideDBParams>
    |> PeptideDBParams.toDomain

