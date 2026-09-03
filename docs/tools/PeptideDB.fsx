(**
---
title: PeptideDB
category: Tools
categoryindex: 1
index: 3
---
*)

(*** hide ***)

(*** condition: prepare ***)
#r "nuget: BioFSharp.Mz, 0.2.1"
#r "nuget: Newtonsoft.Json, 13.0.4"
#r "../../src/ProteomIQon/bin/Release/net10.0/ProteomIQon.dll"

(*** condition: ipynb ***)
#if IPYNB
#r "nuget: ProteomIQon, {{fsdocs-package-version}}"
#endif // IPYNB

(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# PeptideDB

PeptideDB digests the proteins of a FASTA file in silico and stores the resulting peptides, with their masses and modifications, in a SQLite database.
The search tools compare measured spectra against these peptides, so the digestion settings define the search space of the whole analysis. How the protease, the allowed missed cleavages, the length limits and the modifications shape that space is explained in the [BioFSharp.Mz documentation](https://www.biofsharp.com/BioFSharp.Mz/02_02_search_databases.html#Describing-the-search-space), and the [database layout](https://www.biofsharp.com/BioFSharp.Mz/02_02_search_databases.html#Building-and-connecting-to-the-database) is described on the same page.

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | one FASTA file (no directories) | your proteome |
| `-o` | the output directory, created when missing | |
| `-p` | the parameter file in JSON | this page |

The tool writes `<Name>.db` into the output directory, where Name comes from the parameter file. The database is read by [PeptideSpectrumMatching]({{root}}tools/PeptideSpectrumMatching.html), [PSMStatistics]({{root}}tools/PSMStatistics.html), [MsFraggerToPSM]({{root}}tools/MsFraggerToPSM.html), [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html), [PSMBasedQuantificationTIMs]({{root}}tools/PSMBasedQuantificationTIMs.html), [AlignmentBasedQuantification]({{root}}tools/AlignmentBasedQuantification.html) and [ProteinInference]({{root}}tools/ProteinInference.html).
When a `<Name>.db` already exists in the output directory and the parameters stored inside it equal the parameter file, the tool reuses that database and only makes sure its index is present.
The output directory also receives `PeptideDB_log.txt` and `PeptideDB_<Name>_log.txt`.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| Name | "AraTest" | Name of the database and stem of the output file. |
| ParseProteinIDRegexPattern | "id" | How protein identifiers are taken from the FASTA headers. "id" (in any letter case, or an empty string) keeps the whole header line as the identifier. Any other value is a regular expression, and the matched part of the header becomes the identifier, for example "^\S+" for the first word. |
| Protease | Protease.Trypsin | Enzyme used for the digestion. Other values are Protease.Trypsin_P, Protease.LysC, Protease.LysC_P, Protease.Chymotrypsin and Protease.PepsinA. |
| MinMissedCleavages | 0 | Lowest number of missed cleavage sites a peptide may have. |
| MaxMissedCleavages | 2 | Highest number of missed cleavage sites a peptide may have. |
| MaxMass | 15000.0 | Peptides above this mass in Da are dropped. |
| MinPepLength | 4 | Shortest peptide length kept, in residues. |
| MaxPepLength | 65 | Longest peptide length kept, in residues. |
| IsotopicMod | [IsotopicMod.N15] | Isotopic labels of the experiment. Each label adds a labeled variant of every peptide next to the unlabeled one. Use [] for an unlabeled experiment. Other values are IsotopicMod.C13, IsotopicMod.O17, IsotopicMod.O18 and IsotopicMod.D. |
| MassMode | MassMode.Monoisotopic | Mass calculation, MassMode.Monoisotopic or MassMode.Average. |
| FixedMods | [] | Modifications applied to every matching residue. |
| VariableMods | [Modification.Oxidation'Met'; Modification.Acetylation'ProtNTerm'] | Modifications that may or may not be present. Every allowed combination becomes its own database entry. |
| VarModThreshold | 4 | Maximum number of variable modifications on one peptide. |

The available modifications are the cases of ProteomIQon.Modification, for example Modification.Carbamidomethyl'Cys' or Modification.Phosphorylation'Ser'Thr'Tyr'.
The default file is [peptideDBParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/peptideDBParams.json). A variant for Thermo data with carbamidomethylation as variable modification is [peptideDBParamsThermo.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/peptideDBParamsThermo.json).

## Writing a parameter file
*)

open BioFSharp.Mz.SearchDB
open ProteomIQon

let peptideDBParams : Dto.PeptideDBParams =
    {
        Name                       = "AraTest"
        ParseProteinIDRegexPattern = "id"
        Protease                   = Protease.Trypsin
        MinMissedCleavages         = 0
        MaxMissedCleavages         = 2
        MaxMass                    = 15000.0
        MinPepLength               = 4
        MaxPepLength               = 65
        IsotopicMod                = [IsotopicMod.N15]
        MassMode                   = MassMode.Monoisotopic
        FixedMods                  = []
        VariableMods               = [Modification.Oxidation'Met'; Modification.Acetylation'ProtNTerm']
        VarModThreshold            = 4
    }

// Replace the temp folder with your project folder.
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "peptideDBParams.json")

Json.serializeAndWrite outputPath peptideDBParams

(**
## Running the tool

Install the tool with `dotnet tool install --global ProteomIQon.PeptideDB`, then build the database:

```text
proteomiqon-peptidedb -i path/to/proteome.fasta -o path/to/output -p path/to/peptideDBParams.json
```

All flags:

```text
proteomiqon-peptidedb --help
```
*)
