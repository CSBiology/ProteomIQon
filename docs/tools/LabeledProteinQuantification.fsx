(**
---
title: LabeledProteinQuantification
category: Tools
categoryindex: 1
index: 16
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

# LabeledProteinQuantification

In a 15N metabolic labeling experiment one sample grows on 15N, is mixed with the 14N sample, and both versions of every peptide are measured in the same run, so the ratio of the light and the heavy peak is a direct comparison of the two samples. How the label enters the search is described in [BioFSharp.Mz](https://www.biofsharp.com/BioFSharp.Mz/02_02_search_databases.html#Describing-the-search-space), and how the ratio is formed from the two peaks in the [quantification chapter](https://www.biofsharp.com/BioFSharp.Mz/03_01_quantification.html#From-one-peak-to-a-pipeline). LabeledProteinQuantification takes the .quantAndProt rows from [JoinQuantPepIonsWithProteins]({{root}}tools/JoinQuantPepIonsWithProteins.html) and aggregates them to one light intensity, one heavy intensity and one ratio per protein group and run. It works in up to four steps. The first step is mandatory: every row already holds a light and a heavy intensity, the step divides them into Ratio_LightByHeavy and merges the rows identified from the light and from the heavy form of the same peptide ion. The second step merges charge states and the third merges modified forms of a peptide, both optional. The fourth merges all peptides of a protein group. Every step transforms, filters and aggregates light, heavy and ratio separately.

<img src="{{root}}img/LabeledQuant.png" width="1000" height="750" />

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | one or more `.quantAndProt` files, or a directory that is searched for `*.quantAndProt` | [JoinQuantPepIonsWithProteins]({{root}}tools/JoinQuantPepIonsWithProteins.html) |
| `-o` | the output directory, created when missing | |
| `-p` | the parameter file in JSON | this page |

Logs go to `LabeledProteinQuantification_log.txt` and `LabeledQuantification_log.txt` in the output directory.

The tool writes up to five tab separated tables into the output directory. GlobModAggregation.txt is the result of the first step, one row per peptide ion and run with Quant_Light, Quant_Heavy and Ratio_LightByHeavy. ChargeAggregation.txt appears only when the charge step runs and ModificationAggregation.txt only when the modification step runs. ProteinAggregation.txt has one row per protein group and run with the aggregated light, heavy and ratio values, and for each the number of values used, CV, standard deviation and SEM. LabeledQuant.txt is the same table pivoted, one row per protein group and one block of columns per input file, named `<run>.Ratio_LightByHeavy`, `<run>.Quant_Heavy` and so on. [RatioLFQ]({{root}}tools/RatioLFQ.html) reads this file to turn the ratios into comparable intensities per run.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| Correlation_Light_Heavy_Threshold | None | Some 0.9 keeps only rows whose light and heavy elution profiles correlate above 0.9. None keeps everything. The filter is applied in the first step. |
| Alignment_QValue | None | Rows that came from alignment and have an alignment q-value at or above this value are dropped, None keeps everything. |
| ModificationFilter | UseModifiedPeptides.All | Which peptides enter the aggregation. All keeps every peptide, No drops every peptide with a modification, UseOnly mods keeps unmodified peptides and peptides whose modifications are all in the list. |
| AggregateGlobalModificationsParams | { LabeledTransform = None; LabeledSingleFilters = None; LabeledGroupFilters = None; LabeledAggregation = { Light = Mean; Heavy = Mean; Ratio = Mean } } | The first step, always run. |
| AggregatePeptideChargeStatesParams | None | Some params runs the charge step with the given AggregationParams. None skips it. |
| AggregateModifiedPeptidesParams | None | Some params runs the modification step. None skips it. |
| AggregateToProteinGroupsParams | { LabeledTransform = None; LabeledSingleFilters = None; LabeledGroupFilters = None; LabeledAggregation = { Light = Mean; Heavy = Mean; Ratio = Mean } } | The protein step, always run. |

An AggregationParams record describes one step. Each of its fields has a Light, a Heavy and a Ratio entry, so the three value series can be treated differently.

| Parameter | Default | Meaning |
|---|---|---|
| LabeledTransform | None | Some { Light = Some Log2; Heavy = Some Log2; Ratio = Some Log2 } applies a NumericTransform (Log2, Add, Substract, MultiplyBy, DivideBy) before filtering and aggregation. None inside means no transform for that series. |
| LabeledSingleFilters | None | Some { Light = None; Heavy = None; Ratio = Some (seq [IsSmallerThan 0.0]) } keeps only values that pass every NumericFilter (IsBiggerThan, IsSmallerThan) of their series. |
| LabeledGroupFilters | None | Some { Light = Some (seq [Tukey 1.5]); Heavy = None; Ratio = None } removes outliers within each group before aggregation. GroupFilter is Tukey factor, Stdev factor or TopX count. |
| LabeledAggregation | { Light = Mean; Heavy = Mean; Ratio = Mean } | How the remaining values of a group are combined per series. NumericAggregation is Mean, Median or Sum. |

Five default files exist. [LabeledQuantificationParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabeledQuantificationParams.json) is the table above. [LabeledQuantificationParams_withCorrFilter.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabeledQuantificationParams_withCorrFilter.json) sets Correlation_Light_Heavy_Threshold to Some 0.9. [LabeledQuantificationParams_CorrFilter_ChargeAgg.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabeledQuantificationParams_CorrFilter_ChargeAgg.json) sets the threshold to Some 0.0 and adds the charge step with Mean. [LabeledQuantificationParams_CorrFilter_ChargeAgg_ModAgg.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabeledQuantificationParams_CorrFilter_ChargeAgg_ModAgg.json) sets the threshold to Some 0.0 and adds the charge and the modification step, both with Mean. [LabeledQuantificationParams_Transform_Filter_Sum.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabeledQuantificationParams_Transform_Filter_Sum.json) log2 transforms light, heavy and ratio in the first step, keeps only log2 ratios below 0.0 there, and sums light and heavy in the protein step while the ratio is still averaged.

## Writing a parameter file
*)

open System.IO
open ProteomIQon

// Used for the first step and the protein step of the default file: no transform, no filters, mean of light, heavy and ratio.
let meanOfAll : Common.LabeledProteinQuantification.AggregationParams =
    {
        LabeledTransform     = None
        LabeledSingleFilters = None
        LabeledGroupFilters  = None
        LabeledAggregation   =
            {
                Light = NumericAggregation.Mean
                Heavy = NumericAggregation.Mean
                Ratio = NumericAggregation.Mean
            }
    }

let labeledQuantificationParams : Dto.LabeledQuantificationParams =
    {
        Correlation_Light_Heavy_Threshold  = None
        Alignment_QValue                   = None
        ModificationFilter                 = UseModifiedPeptides.All
        AggregateGlobalModificationsParams = meanOfAll
        AggregatePeptideChargeStatesParams = None
        AggregateModifiedPeptidesParams    = None
        AggregateToProteinGroupsParams     = meanOfAll
    }

// Replace the temp folder with your project folder.
let outputPath = Path.Combine(Path.GetTempPath(), "LabeledQuantificationParams.json")

Json.serializeAndWrite outputPath labeledQuantificationParams

(**
## Running the tool

The tool installs with `dotnet tool install --global ProteomIQon.LabeledProteinQuantification`. A single run:

```text
proteomiqon-labeledproteinquantification -i path/to/run.quantAndProt -o path/to/output -p path/to/LabeledQuantificationParams.json
```

Several runs in one table, from a list or from a directory:

```text
proteomiqon-labeledproteinquantification -i path/to/run1.quantAndProt path/to/run2.quantAndProt -o path/to/output -p path/to/LabeledQuantificationParams.json
proteomiqon-labeledproteinquantification -i path/to/quantAndProt -o path/to/output -p path/to/LabeledQuantificationParams.json
```

All flags:

```text
proteomiqon-labeledproteinquantification --help
```
*)
