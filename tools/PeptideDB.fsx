(**
// can't yet format YamlFrontmatter (["title: PeptideDB"; "category: Tools"; "categoryindex: 1"; "index: 1"], Some { StartLine = 2 StartColumn = 0 EndLine = 6 EndColumn = 8 }) to pynb markdown

[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/PeptideDB.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/PeptideDB.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/PeptideDB.ipynb)

# PeptideDB

MS-based shotgun proteomics estimates protein abundances using a proxy: peptides. An established method to identify acquired MS/MS spectra is the comparison of each spectrum with peptides in a reference database. 
The PeptideDB tool helps to create peptide databases by in silico digestion given proteome information in the FASTA format and a set of parameters that allow the user to mimic conditions of their specific experiment. 
The created database stores peptide protein relationships in a SQLite database which can then be supplied to other ProteomIQon tools.

## Parameters
The following table gives an overview of the parameter set:

| **Parameter**              | **Default Value**                                                 | **Description**                                                     |
|----------------------------|-------------------------------------------------------------------|---------------------------------------------------------------------|
| Name                       | "YourNameHere"                                                    | Name of the database                                                |
| ParseProteinIDRegexPattern | id                                                                | Regex pattern for parsing of the protein IDs in the database        |
| Protease                   | Protease.Trypsin                                                  | Protease used for the digestion of the proteins                     |
| MinMissedCleavages         | 0                                                                 | Minimal amount of missed cleavages a peptide can have               |
| MaxMissedCleavages         | 2                                                                 | Maximal amount of missed cleavages a peptide can have               |
| MaxMass                    | 15000.                                                            | Maximal mass of a peptide in Da                                     |
| MinPepLength               | 4                                                                 | Minimal length of a peptide                                         |
| MaxPepLength               | 65                                                                | Maximal length of a peptide                                         |
| IsotopicMod                | [IsotopicMod.N15]                                                 | List of isotopic modifications in the experiment                    |
| MassMode                   | MassMode.Monoisotopi                                              | Method for mass calcluation (possibilities: Monoisotopic & Average) |
| FixedMods                  | []                                                                | Fixed modifications of the proteins                                 |
| VariableMods               | [Modification.Oxidation'Met';Modification.Acetylation'ProtNTerm'] | Variable Modifications of the proteins                              |
| VarModThreshold            | 4                                                                 | Threshold for variable modifications                                |

## Parameter Generation

Parameters are handed to the cli tool as .json files. You can download an example file [here](https://github.com/CSBiology/ProteomIQon/blob/master/src/ProteomIQon/defaultParams/peptideDBParams.json), 
or use an F# script, which can be downloaded or run in Binder at the top of the page, to write your own parameter file:

*)
#r "nuget: BioFSharp.Mz, 0.1.5-beta"
#r "nuget: Newtonsoft.Json, 12.0.3"
#r "nuget: ProteomIQon, 0.0.1"

open BioFSharp.Mz.SearchDB
open Newtonsoft.Json
open ProteomIQon

let peptideDBParams: Dto.PeptideDBParams = 
    {
    Name                        = "Test"
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
    VariableMods                = [Modification.Oxidation'Met';Modification.Acetylation'ProtNTerm']
    VarModThreshold             = 4
    }

let serialized = 
    peptideDBParams
    |> JsonConvert.SerializeObject

System.IO.File.WriteAllText("path/to/your/params.json",serialized)
(**
## Executing the Tool
To create a peptide data base just call the tool:

	proteomiqon-peptidedb -i "path/to/your/proteom.fasta" -o "path/to/your/outDirectory" -p "path/to/your/params.json"

A detailed description of the CLI arguments the tool expects can be obtained by calling the tool:

	proteomiqon-peptidedb --help

*)