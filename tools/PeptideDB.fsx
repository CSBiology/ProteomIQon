(**
// can't yet format YamlFrontmatter (["title: PeptideDB"; "category: Tools"; "categoryindex: 1"; "index: 1"], Some { StartLine = 2 StartColumn = 0 EndLine = 6 EndColumn = 8 }) to pynb markdown

# PeptideDB

This tool takes proteome information in the form of a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file. The proteome information is then processed according 
to the parameters given and stored in a database. The database contains information about the proteins, their identifier and their digested peptides.

## Parameters

| Parameter                  | Default Value                                                     | Description                                                         |
|----------------------------|-------------------------------------------------------------------|---------------------------------------------------------------------|
| Name                       | "YourNameHere"                                                    | Name of the database                                                |
| ParseProteinIDRegexPattern | id                                                                | Regex pattern for parsing of the protein IDs in the database        |
| Protease                   | Protease.Trypsin                                                  | Protease used for the digestion of the proteins                     |
| MinMissedCleavages         | 0                                                                 | Minimal amount of missed cleavages a peptide can have               |
| MaxMissedCleavages         | 2                                                                 | Maximal amount of missed cleavages a peptide can have               |
| MaxMass                    | 15000                                                             | Maximal mass of a peptide in Da                                     |
| MinPepLength               | 4                                                                 | Minimal length of a peptide                                         |
| MaxPepLength               | 65                                                                | Maximal length of a peptide                                         |
| IsotopicMod                | [IsotopicMod.N15]                                                 | List of isotopic modifications in the experiment                    |
| MassMode                   | MassMode.Monoisotopi                                              | Method for mass calcluation (possibilities: Monoisotopic & Average) |
| FixedMods                  | []                                                                | Fixed modifications of the proteins                                 |
| VariableMods               | [Modification.Oxidation'Met';Modification.Acetylation'ProtNTerm'] | Variable Modifications of the proteins                              |
| VarModThreshold            | 4                                                                 | Threshold for variable modifications                                |


*)