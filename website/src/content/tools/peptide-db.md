---
name: "PeptideDB"
phase: preprocessing
subgroup: Peptide database
order: 3
summary: "PeptideDB digests the proteins of a FASTA file in silico and stores the resulting peptides, with their masses and modifications, in a SQLite database."
inputsFrom: []
---

## Tool Description

PeptideDB digests the proteins of a FASTA file in silico and stores the resulting peptides, with their masses and modifications, in a SQLite database.

## How does PeptideDB work?

1. **Search space**

   The search tools compare measured spectra against these peptides, so the digestion settings define the search space of the whole analysis.

2. **Algorithm**

   The protease, missed-cleavage limits, peptide lengths and modifications define the search space. BioFSharp.Mz describes the [search-space settings](https://www.biofsharp.com/BioFSharp.Mz/02_02_search_databases.html#Describing-the-search-space) and the [database layout](https://www.biofsharp.com/BioFSharp.Mz/02_02_search_databases.html#Building-and-connecting-to-the-database).

## Inputs

These inputs are needed to execute the tool:

| Flag | Meaning | Comes from |
| --- | --- | --- |
| `-i` | one FASTA file (no directories) | choosen Proteome |
| `-o` | the output directory, created when missing | - |
| `-p` | the parameter file in JSON | [Github ProteomIQon Core](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/peptideDBParams.json) |

- The tool writes `<Name>.db` into the output directory, using `Name` from the parameter file. [PeptideSpectrumMatchingTIMs](./peptide-spectrum-matching-tims.md) requires a database with monoisotopic masses.

- The database is read by

  - [PeptideSpectrumMatching](./peptide-spectrum-matching.md)
  - [PeptideSpectrumMatchingTIMs](./peptide-spectrum-matching-tims.md)
  - [PSMStatistics](./psm-statistics.md)
  - [MSFraggerToPSM](./ms-fragger-to-psm.md)
  - [PSMBasedQuantification](./psm-based-quantification.md)
  - [PSMBasedQuantificationTIMs](./psm-based-quantification-tims.md)
  - [AlignmentBasedQuantification](./alignment-based-quantification.md)
  - [Protein Inference](./protein-inference.md)

- When a `<Name>.db` already exists in the output directory and the parameters stored inside it equal the parameter file, the tool reuses that database and only makes sure its index is present.

- The output directory also receives `PeptideDB_log.txt` and `PeptideDB_<Name>_log.txt`.

## Parameter file

| Parameter | Default | Meaning |
|---|---|---|
| Name | "AraTest" | Name of the database and stem of the output file. |
| ParseProteinIDRegexPattern | "id" | How protein identifiers are taken from the FASTA headers. "id" (in any letter case, or an empty string) keeps the whole header line as the identifier. Any other value is a regular expression, and the matched part of the header becomes the identifier, for example "^\S+" for the first word. |
| Protease | Protease.Trypsin | Enzyme used for the digestion. Other values are Protease.Trypsin_P, Protease.LysC, Protease.LysC_P, Protease.Chymotrypsin and Protease.PepsinA. |
| MinMissedCleavages | 0 | Lowest number of missed cleavage sites a peptide may have. |
| MaxMissedCleavages | 2 | Highest number of missed cleavage sites a peptide may have. |
| MaxMass | 15000.0 | Peptides above this mass in Da are dropped. |
| MinPepLength | 6 | Shortest peptide length kept, in residues. |
| MaxPepLength | 65 | Longest peptide length kept, in residues. |
| IsotopicMod | [IsotopicMod.N15] | Isotopic labels of the experiment. Each label adds a labeled variant of every peptide next to the unlabeled one. Use [] for an unlabeled experiment. Other values are IsotopicMod.C13, IsotopicMod.O17, IsotopicMod.O18 and IsotopicMod.D. |
| MassMode | MassMode.Monoisotopic | Mass calculation, MassMode.Monoisotopic or MassMode.Average. |
| FixedMods | [] | Modifications applied to every matching residue. |
| VariableMods | [Modification.Oxidation'Met'; Modification.Acetylation'ProtNTerm'] | Modifications that may or may not be present. Every allowed combination becomes its own database entry. |
| VarModThreshold | 4 | Maximum number of variable modifications on one peptide. |

- The available modifications are the cases of `ProteomIQon.Modification`, for example `Modification.Carbamidomethyl'Cys'` or `Modification.Phosphorylation'Ser'Thr'Tyr'`.
- The default file is [peptideDBParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/peptideDBParams.json). The [Thermo preset](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/peptideDBParamsThermo.json) includes carbamidomethylation as a variable modification.

## Parameter generation

For an F# example that writes a parameter JSON file, see the [parameter-generation example](https://github.com/CSBiology/ProteomIQon/blob/dev/docs/tools/PeptideDB.fsx#L65).

## Installing and executing the tool

1. **Installing**

   Install the Tool with the following command inside the choosen Command-Line

   ```sh
   dotnet tool install --global ProteomIQon.PeptideDB
   ```

2. **Execute**

     ```sh
     proteomiqon-peptidedb -i path/to/proteome.fasta -o path/to/output -p path/to/peptideDBParams.json
     ```

3. **All Flags**

   ```sh
   proteomiqon-peptidedb --help
   ```
