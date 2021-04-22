(**
---
title: PeptideDB
category: Tools
categoryindex: 1
index: 1
---
*)
(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# PeptideDB

This tool takes proteome information in the form of a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file. The proteome information is then processed according 
to the parameters given and stored in a database. The database contains information about the proteins, their identifier and their digested peptides.

## Parameters

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

Parameters are handed to the cli tool as a .json file. you can download the default file [here](https://github.com/CSBiology/ProteomIQon/blob/master/src/ProteomIQon/defaultParams/peptideDBParams.json), 
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

System.IO.File.WriteAllText("YourPathHere",serialized)

(***condition:ipynb***)
#if IPYNB
(**
If you are running this tool in Binder, you can copy the output of the following codeblock and save it in a JSON file.
*)
serialized
#endif // IPYNB