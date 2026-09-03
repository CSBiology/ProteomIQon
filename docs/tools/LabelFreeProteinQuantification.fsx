(**
---
title: LabelFreeProteinQuantification
category: Tools
categoryindex: 1
index: 15
---
*)

(*** hide ***)

(*** condition: prepare ***)
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

# LabelFreeProteinQuantification

After [JoinQuantPepIonsWithProteins]({{root}}tools/JoinQuantPepIonsWithProteins.html) every quantified peptide ion carries its protein group. LabelFreeProteinQuantification turns these ion intensities into one intensity per protein group and run. It works in up to three steps. The first merges the charge states of a peptide, the second merges the modified forms of a peptide, and the third merges all peptides of a protein group. The first two steps are optional. Each step has the same shape: transform the values, drop single values or outliers within a group, then aggregate what is left. Only light rows (GlobalMod 0) go in, so this tool is for unlabeled samples. For labeled samples use [LabeledProteinQuantification]({{root}}tools/LabeledProteinQuantification.html).

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | one or more `.quantAndProt` files, or a directory that is searched for `*.quantAndProt` | [JoinQuantPepIonsWithProteins]({{root}}tools/JoinQuantPepIonsWithProteins.html) |
| `-o` | the output directory, created when missing | |
| `-p` | the parameter file in JSON | this page |

Logs go to `LabeledProteinQuantification_log.txt` and `LabeledQuantification_log.txt` in the output directory.

The tool writes up to four tab separated tables into the output directory. ChargeAggregation.txt appears only when the charge step runs and ModificationAggregation.txt only when the modification step runs. ProteinAggregation.txt has one row per protein group and run with the aggregated intensity, the number of values that went into it, and its CV, standard deviation and SEM. LabelFreeQuant.txt is the same table pivoted, one row per protein group and one block of columns per input file, named `<run>.Quant_Light`, `<run>.ItemsUsedForQuant_Light` and so on. Protein groups missing in a run have empty cells there.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| Alignment_QValue | None | Rows that came from alignment and have an alignment q-value at or above this value are dropped, None keeps everything. |
| ModificationFilter | UseModifiedPeptides.All | Which peptides enter the aggregation. All keeps every peptide, No drops every peptide with a modification, UseOnly mods keeps unmodified peptides and peptides whose modifications are all in the list. |
| AggregatePeptideChargeStatesParams | None | Some params runs the charge step with the given AggregationParams. None skips it. |
| AggregateModifiedPeptidesParams | None | Some params runs the modification step. None skips it. |
| AggregateToProteinGroupsParams | { Transform = None; SingleFilters = None; GroupFilters = None; Aggregation = { Light = Mean } } | The protein step, always run. |

An AggregationParams record describes one step.

| Parameter | Default | Meaning |
|---|---|---|
| Transform | None | Some { Light = Some Log2 } applies a NumericTransform (Log2, Add, Substract, MultiplyBy, DivideBy) to every intensity before filtering and aggregation. |
| SingleFilters | None | Some { Light = Some (seq [IsBiggerThan 4.0]) } keeps only intensities that pass every NumericFilter (IsBiggerThan, IsSmallerThan). |
| GroupFilters | None | Some { Light = Some (seq [Tukey 1.5]) } removes outliers within each group before aggregation. GroupFilter is Tukey factor, Stdev factor or TopX count. |
| Aggregation | { Light = Mean } | How the remaining values of a group are combined. NumericAggregation is Mean, Median or Sum. |

Four default files exist. [LabelFreeQuantificationParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams.json) is the table above and only runs the protein step. [LabelFreeQuantificationParams_ChargeAgg.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_ChargeAgg.json) adds the charge step with Mean. [LabelFreeQuantificationParams_ChargeAgg_ModAgg.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_ChargeAgg_ModAgg.json) adds the charge and the modification step, both with Mean. [LabelFreeQuantificationParams_Transform_Filter_Sum.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_Transform_Filter_Sum.json) runs only the protein step, but log2 transforms the intensities, keeps those above 4.0 and sums them.

## Writing a parameter file
*)

open System.IO
open ProteomIQon

// The protein step of the default file: no transform, no filters, mean of the peptide intensities.
let proteinStep : Common.LabelFreeQuantification.AggregationParams =
    {
        Transform     = None
        SingleFilters = None
        GroupFilters  = None
        Aggregation   = { Light = NumericAggregation.Mean }
    }

let labelFreeQuantificationParams : Dto.LabelFreeQuantificationParams =
    {
        Alignment_QValue                   = None
        ModificationFilter                 = UseModifiedPeptides.All
        AggregatePeptideChargeStatesParams = None
        AggregateModifiedPeptidesParams    = None
        AggregateToProteinGroupsParams     = proteinStep
    }

// Replace the temp folder with your project folder.
let outputPath = Path.Combine(Path.GetTempPath(), "LabelFreeQuantificationParams.json")

Json.serializeAndWrite outputPath labelFreeQuantificationParams

(**
## Running the tool

The tool installs with `dotnet tool install --global ProteomIQon.LabelFreeProteinQuantification`. A single run:

```text
proteomiqon-labelfreeproteinquantification -i path/to/run.quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
```

Several runs in one table, from a list or from a directory:

```text
proteomiqon-labelfreeproteinquantification -i path/to/run1.quantAndProt path/to/run2.quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
proteomiqon-labelfreeproteinquantification -i path/to/quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
```

All flags:

```text
proteomiqon-labelfreeproteinquantification --help
```
*)
