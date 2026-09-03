(**
---
title: JoinQuantPepIonsWithProteins
category: Tools
categoryindex: 1
index: 14
---
*)

(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# JoinQuantPepIonsWithProteins

A .quant file knows the intensity of every peptide ion, and a .prot file knows which protein group each peptide sequence belongs to and how confident that assignment is. JoinQuantPepIonsWithProteins puts the two together. It reads the protein table, splits its semicolon separated peptide list into one row per sequence, and attaches the protein group and its q-value to every quantified ion with the same sequence. The result is one .quantAndProt file per run, which is what [LabelFreeProteinQuantification]({{root}}tools/LabelFreeProteinQuantification.html) and [LabeledProteinQuantification]({{root}}tools/LabeledProteinQuantification.html) read.

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | one or more `.quant` files, or a directory that is searched for `*.quant` | [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html), [PSMBasedQuantificationTIMs]({{root}}tools/PSMBasedQuantificationTIMs.html) or [AlignmentBasedQuantification]({{root}}tools/AlignmentBasedQuantification.html), optionally scored by [AlignmentBasedQuantStatistics]({{root}}tools/AlignmentBasedQuantStatistics.html) |
| `-ii` | one or more `.prot` files, or a directory that is searched for `*.prot` | [ProteinInference]({{root}}tools/ProteinInference.html) or [AddDeducedPeptides]({{root}}tools/AddDeducedPeptides.html) |
| `-o` | the output directory, created when missing | |
| `-mf` | switch: pair files by base name instead of by position | |
| `-c` | number of files processed at the same time, default 1 | |

Quantified ions without a fitted scan time are dropped. The tool uses the ProteinGroup, PeptideSequence and QValue columns of the protein table. The q-value ends up in the ProteinGroup_QValue column of the output.

For every .quant file the tool writes `<name>.quantAndProt`, a tab separated table with all columns of the .quant row plus FileName, ProteinGroup and ProteinGroup_QValue in front. A peptide that maps to several protein groups gives one row per group.

Without `-mf` the two lists are paired by position, so the first .quant file goes with the first .prot file, and both lists need the same length, otherwise the tool stops without output. With `-mf` a .quant file is paired with the .prot file of the same base name, and .quant files without a partner are skipped. Use `-mf` whenever you pass directories.

Delete an existing `<name>.quantAndProt` in the output directory before a rerun, the tool appends to it. Logs go to `JoinQuantPepIonsWithProteins_log.txt` in the output directory.

## Running the tool

The tool installs with `dotnet tool install --global ProteomIQon.JoinQuantPepIonsWithProteins`. A single run:

```text
proteomiqon-joinquantpepionswithproteins -i path/to/run.quant -ii path/to/run.prot -o path/to/output
```

Several runs, paired by position, four at a time:

```text
proteomiqon-joinquantpepionswithproteins -i path/to/run1.quant path/to/run2.quant -ii path/to/run1.prot path/to/run2.prot -o path/to/output -c 4
```

Two directories, paired by file name:

```text
proteomiqon-joinquantpepionswithproteins -i path/to/quant -ii path/to/prot -o path/to/output -mf -c 4
```

All flags:

```text
proteomiqon-joinquantpepionswithproteins --help
```
*)
