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

let defaultPepeParams: Dto.PeptideDBParams = 
        {
        Name                        = "ChlamyTest"
        FastaPath                   = @"C:\Users\david\Source\Repos\netCoreRepos\ProteomiconTest\PeptideDB\Chlamy_JGI5_5(Cp_Mp).fasta"
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
        VariableMods                = [Modification.Oxidation'Met']
        VarModThreshold             = 3
        }
  
let serialized = 
    defaultPepeParams
    |> Json.serialize

System.IO.File.WriteAllText( @"C:\Users\david\Source\Repos\netCoreRepos\ProteomIQon\src\ProteomIQon\defaultParams\peptideDBParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(@"C:\Users\david\Source\Repos\netCoreRepos\ProteomIQon\src\ProteomIQon\defaultParams\peptideDBParams.json")
    |> Json.deserialize<Dto.PeptideDBParams>
    |> PeptideDBParams.toDomain

