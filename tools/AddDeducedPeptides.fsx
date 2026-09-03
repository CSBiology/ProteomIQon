(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/AddDeducedPeptides.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/AddDeducedPeptides.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/AddDeducedPeptides.ipynb)

# AddDeducedPeptides

After [alignment](https://csbiology.github.io/ProteomIQon/tools/QuantBasedAlignment.html) and
[AlignmentBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html) a run's .quant file holds peptides that were never identified
in that run, so the run's .prot file from [ProteinInference](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.html) does not list them. AddDeducedPeptides
takes the protein groups of the combined inference and assigns them to every peptide present in each aligned .quant file, so that
[JoinQuantPepIonsWithProteins](https://csbiology.github.io/ProteomIQon/tools/JoinQuantPepIonsWithProteins.html) can join the transferred quantifications to proteins.

## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | one or more aligned `.quant` files, or a directory that is searched for `*.quant` | [AlignmentBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html)
`-ii` | one or more `.prot` files, or a directory that is searched for `*.prot` | [ProteinInference](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.html)
`-o` | the output directory, created when missing | &#32;


The tool reads all .prot files as one combined inference and does not pair them with the .quant files by name. This only
works when ProteinInference ran over all runs together (GroupFiles set to true), because then every run's .prot file carries the same scores and
the same q-value for a protein group. If a protein group turns up with different scores, the tool stops with the message
"If you are running this tool please infer proteins over all files".

For every input .quant file the tool writes `<quant>.prot` into the output directory, a tab separated table with the columns of a
ProteinInference result. Each protein group keeps its scores and q-value, and its PeptideSequence column lists only the peptides that occur
in that .quant file. Protein groups without any peptide in the file are dropped. JoinQuantPepIonsWithProteins reads the new .prot file
together with the scored .quant file of the same run. There is no parameter file. Delete an existing `<quant>.prot` in the output
directory before a rerun, the tool appends to it. Logs go to `AddDeducedPeptides_log.txt` there.

## Running the tool

Install the tool with `dotnet tool install --global ProteomIQon.AddDeducedPeptides`. One run.

```text
proteomiqon-adddeducedpeptides -i path/to/run.quant -ii path/to/run.prot -o path/to/output
```

All runs of an experiment. Every .prot file in the directory is read, and every .quant file gets a new .prot file.

```text
proteomiqon-adddeducedpeptides -i path/to/quant -ii path/to/prot -o path/to/output
```

Print the description of every argument.

```text
proteomiqon-adddeducedpeptides --help
```

*)