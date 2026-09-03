(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/JoinQuantPepIonsWithProteins.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/JoinQuantPepIonsWithProteins.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/JoinQuantPepIonsWithProteins.ipynb)

# JoinQuantPepIonsWithProteins

A .quant file knows the intensity of every peptide ion, and a .prot file knows which protein group each peptide sequence belongs to and how confident that assignment is. JoinQuantPepIonsWithProteins puts the two together. It reads the protein table, splits its semicolon separated peptide list into one row per sequence, and attaches the protein group and its q-value to every quantified ion with the same sequence. The result is one .quantAndProt file per run, which is what [LabelFreeProteinQuantification](https://csbiology.github.io/ProteomIQon/tools/LabelFreeProteinQuantification.html) and [LabeledProteinQuantification](https://csbiology.github.io/ProteomIQon/tools/LabeledProteinQuantification.html) read.

## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | one or more `.quant` files, or a directory that is searched for `*.quant` | [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html), [PSMBasedQuantificationTIMs](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantificationTIMs.html) or [AlignmentBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html), optionally scored by [AlignmentBasedQuantStatistics](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantStatistics.html)
`-ii` | one or more `.prot` files, or a directory that is searched for `*.prot` | [ProteinInference](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.html) or [AddDeducedPeptides](https://csbiology.github.io/ProteomIQon/tools/AddDeducedPeptides.html)
`-o` | the output directory, created when missing | &#32;
`-mf` | switch: pair files by base name instead of by position | &#32;
`-c` | number of files processed at the same time, default 1 | &#32;


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