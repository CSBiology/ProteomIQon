---
name: "LabelFreeProteinQuantification"
phase: protein-quantification
subgroup: Label-free and labeled quantification
order: 2
summary: "LabelFreeProteinQuantification turns peptide ion intensities into one intensity per protein group and run."
inputsFrom: [join-quant-with-peptide-ions]
---

## Tool Description

After [JoinQuantPepIonsWithProteins](./join-quant-with-peptide-ions.md) every quantified peptide ion carries its protein group. LabelFreeProteinQuantification turns these ion intensities into one intensity per protein group and run.

## How does Label-Free Protein Quantification work?

1. **Aggregate peptide charge states (optional)**

   The first step merges the charge states of a peptide

2. **Aggregate modified peptides (optional)**

   The second step merges the modified forms of a peptide

3. **Aggregate to protein groups**

   The third step merges all peptides of a protein group

The first two steps are optional. Each step has the same shape: transform the values, drop single values or outliers within a group, then aggregate what is left. Only light rows (GlobalMod 0) go in, so this tool is for unlabeled samples. For labeled samples use [LabeledProteinQuantification](./labeled-protein-quantification.md).

### What happens when stages are skipped?

If charge-state aggregation is disabled, differently charged versions of the same peptide remain separate for the next stage. If modification aggregation is disabled, differently modified versions remain separate. If both are disabled, protein-group aggregation works directly on the unaggregated peptide ions, so a peptide sequence may still occur with multiple charge states and modifications.

## Inputs

| Flag | Meaning | Comes from |
| --- | --- | --- |
| `-i` | One or more `.quantAndProt` files, or a directory containing those files. Each file contains peptide quantification and protein-group assignments. | [JoinQuantPepIonsWithProteins](./join-quant-with-peptide-ions.md) |
| `-o` | Directory in which the aggregation results are written. It is created when missing. | Chosen output location |
| `-p` | the parameter file in JSON | [Default parameter file](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams.json) |

**Important:** The input must already combine peptide-level quantification with protein inference. The tool does not perform those upstream steps itself.

## Parameter file

### General settings

| Parameter | Default | Description |
| --- | --- | --- |
| Alignment_QValue | `Some 0.05`| Rows that came from alignment and have an alignment Q-value at or above this value are dropped. [PSM based quantifications](./psm-based-quantification.md) carry the Q-value 0, and rows without an alignment Q-value are kept. None keeps every row. |
| ModificationFilter | `UseModifiedPeptides.All` | Which peptides enter the aggregation. All keeps every peptide, No drops every peptide with a modification, UseOnly mods keeps unmodified peptides and peptides whose modifications are all in the list. | 
| AggregatePeptideChargeStatesParams | `None` | Some params runs the charge step with the given `AggregationParams`. `None` skips it. | 
| AggregateModifiedPeptidesParams | `None` | Some params runs the modification step. None skips it. | 
| AggregateToProteinGroupsParams | `{ Transform = None; SingleFilters = None; GroupFilters = None; Aggregation = { Light = Mean } }` | The protein step, always run.| 

An `AggregationParams` record describes one step.

|Parameter |Default | Description | 
| --- | --- | --- |
| Transform  | `None` | Some { Light = Some Log2 } applies a NumericTransform (Log2, Add, Substract, MultiplyBy, DivideBy) to every intensity before filtering and aggregation. |
| SingleFilters |  `None`  | Some { Light = Some (seq [IsBiggerThan 4.0]) } keeps only intensities that pass every NumericFilter (IsBiggerThan, IsSmallerThan). | 
| GroupFilters | `None` | Some { Light = Some (seq [Tukey 1.5]) } removes outliers within each group before aggregation. GroupFilter is Tukey factor, Stdev factor or TopX count. | 
|Aggregation | `{ Light = Mean }`| How the remaining values of a group are combined. NumericAggregation is Mean, Median or Sum.


All four presets use `Alignment_QValue = Some 0.05`.

1. [LabelFreeQuantificationParams.json](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams.json) only runs the protein step.

2. [LabelFreeQuantificationParams_ChargeAgg.json](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_ChargeAgg.json) adds the charge step with Mean.

3. [LabelFreeQuantificationParams_ChargeAgg_ModAgg.json](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_ChargeAgg_ModAgg.json) adds the charge and the modification step, both with Mean. 

4. [LabelFreeQuantificationParams_Transform_Filter_Sum.json](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_Transform_Filter_Sum.json) runs only the protein step, but log2 transforms the intensities, keeps those above 4.0 and sums them. 

## Outputs

The tool writes tab-separated tables.

| File | Contents | When produced |
| --- | --- | --- |
| `ChargeAggregation.txt` | Intensities combined across charge states of each modified peptide | When charge-state aggregation is enabled |
| `ModificationAggregation.txt` | Intensities combined across modified forms of each peptide | When modification aggregation is enabled |
| `ProteinAggregation.txt` | One row per protein group and run with aggregated intensity, number of values used, coefficient of variation, standard deviation and standard error of the mean | Always |
| `LabelFreeQuant.txt` | The same results in wide format, with one row per protein group and a set of columns per run, such as `<run>.Quant_Light` and `<run>.ItemsUsedForQuant_Light`. Missing protein groups have empty cells for the corresponding runs. | Always |

Logs go to `LabeledProteinQuantification_log.txt` and `LabeledQuantification_log.txt` in the output directory.

## Parameter generation

From the ProteomIQon repository root, build the library with `dotnet build src/ProteomIQon/ProteomIQon.fsproj -c Release`. Save the following as `generate-params.fsx` in that same directory and run `dotnet fsi generate-params.fsx`. The script writes the JSON file to the system temporary directory.

```fsharp
#r "nuget: Newtonsoft.Json, 13.0.4"
#r "./src/ProteomIQon/bin/Release/net10.0/ProteomIQon.dll"

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
        Alignment_QValue                   = Some 0.05
        ModificationFilter                 = UseModifiedPeptides.All
        AggregatePeptideChargeStatesParams = None
        AggregateModifiedPeptidesParams    = None
        AggregateToProteinGroupsParams     = proteinStep
    }

// Replace the temp folder with your project folder.
let outputPath = Path.Combine(Path.GetTempPath(), "LabelFreeQuantificationParams.json")

Json.serializeAndWrite outputPath labelFreeQuantificationParams
```

## Running the tool
1. **Installing**

   Install the tool globally with the .NET CLI.

   ```sh
   dotnet tool install --global ProteomIQon.LabelFreeProteinQuantification
   ```

2. **Quantify one run**

   ```sh
   proteomiqon-labelfreeproteinquantification -i path/to/run.quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
   ```

3. **Several runs in one table, from a list or from a directory**


   ```sh
   proteomiqon-labelfreeproteinquantification -i path/to/run1.quantAndProt path/to/run2.quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
   proteomiqon-labelfreeproteinquantification -i path/to/quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
   ```

4. **All Flags**

   ```sh
   proteomiqon-labelfreeproteinquantification --help
   ```
