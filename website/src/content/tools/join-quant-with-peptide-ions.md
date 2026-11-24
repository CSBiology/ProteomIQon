---
name: "JoinQuantPepIonsWithProteins"
phase: protein-quantification
subgroup: Join quantification
order: 1
summary: "JoinQuantPepIonsWithProteins joins the peptide ion intensities of a .quant file with the protein groups of a .prot file into one .quantAndProt file per run."
inputsFrom: [psm-based-quantification, psm-based-quantification-tims, alignment-based-quantification, alignment-based-quant-statistics, protein-inference, add-deduced-peptides]
---

## JoinQuantPepIonsWithProteins

A `.quant` file contains the intensities of peptide ions. A `.prot` file provides their protein-group assignments and the associated confidence. JoinQuantPepIonsWithProteins combines these files by expanding the semicolon-separated peptide sequences in the protein table and matching them to quantified peptide ions. Matching ions receive their protein group and q-value.

One `.quantAndProt` file is produced for each run. These files are inputs for [LabelFreeProteinQuantification](./label-free-protein-quantification.md) and [LabeledProteinQuantification](./labeled-protein-quantification.md).

## Inputs

| Flag | Meaning | Comes from |
| --- | --- | --- |
| `-i` | One or more `.quant` files, or a directory containing `*.quant` files. | [PSMBasedQuantification](./psm-based-quantification.md), [PSMBasedQuantificationTIMs](./psm-based-quantification-tims.md) or [AlignmentBasedQuantification](./alignment-based-quantification.md), optionally scored by [AlignmentBasedQuantStatistics](./alignment-based-quant-statistics.md) |
| `-ii` | One or more `.prot` files, or a directory containing `*.prot` files. | [ProteinInference](./protein-inference.md) or [AddDeducedPeptides](./add-deduced-peptides.md) |
| `-o` | Output directory, which is created if it does not already exist. |  |
| `-mf` | Switch to match input files by their base names rather than their positions. |  |
| `-c` | Number of files processed concurrently; default is `1`. |  |

### Quantified ions

Quantified ions with no fitted scan time are excluded. From the protein table the tool reads `ProteinGroup`, `PeptideSequence` and `QValue`. The protein q-value is written to the `ProteinGroup_QValue` output column.

### For every .quant file

Each `.quant` input produces a tab-delimited `<name>.quantAndProt` file. The columns from the original `.quant` rows are retained, with `FileName`, `ProteinGroup` and `ProteinGroup_QValue` placed before them. Peptides assigned to multiple protein groups are represented by one output row for each group.

### Without -mf / With -mf

Without `-mf`, the lists of `.quant` and `.prot` files are matched in their given order. Both lists must have the same number of files; otherwise nothing is written.

With `-mf`, matching is based on the files' base names. A `.quant` file is skipped when no corresponding `.prot` file exists. Use `-mf` when supplying directories.

### Logs

Existing `<name>.quantAndProt` files must be deleted before another run because new rows are appended. The output directory also contains `JoinQuantPepIonsWithProteins_log.txt`.

## Running the tool

The tool installs with:

```text
dotnet tool install --global ProteomIQon.JoinQuantPepIonsWithProteins
```

### A single run

```text
proteomiqon-joinquantpepionswithproteins -i path/to/run.quant -ii path/to/run.prot -o path/to/output
```

### Several runs, paired by position, four at a time

```text
proteomiqon-joinquantpepionswithproteins -i path/to/run1.quant path/to/run2.quant -ii path/to/run1.prot path/to/run2.prot -o path/to/output -c 4
```

### Two directories, paired by file name

```text
proteomiqon-joinquantpepionswithproteins -i path/to/quant -ii path/to/prot -o path/to/output -mf -c 4
```

### All flags

```text
proteomiqon-joinquantpepionswithproteins --help
```
