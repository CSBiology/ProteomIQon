---
name: "Label-Free Protein Quantification"
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
| `-i` | One or more `.quantAndProt` input files containing peptide quantification and protein-group assignments. | [JoinQuantPepIonsWithProteins](./join-quant-with-peptide-ions.md) |
| `-o` | Directory in which the aggregation results are written. | Chosen output location |
| `-p` | Path to a JSON file specifying modification filtering and aggregation settings. | [Default parameter file](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams.json) |

**Important:** The input must already combine peptide-level quantification with protein inference. The tool does not perform those upstream steps itself.

## Parameter file

### General settings

| Parameter | Default | Description |
| --- | --- | --- |
| Alignment_QValue | `Some 0.05`| Rows that came from `alignment` and have an alignment Q-value at or above this value are dropped. [PSM based quantifications](./psm-based-quantification.md) carry the Q-value 0, and rows without an alignment Q-value are kept. None keeps every row. | 
| ModificationFilter | `UseModifiedPeptides.All` | Which peptides enter the aggregation. All keeps every peptide, No drops every peptide with a modification, UseOnly mods keeps unmodified peptides and peptides whose modifications are all in the list. | 
| AggregatePeptideChargeStatesParams | `None` | Some params runs the charge step with the given `AggregationParams`. `None` skips it. | 
| AggregateModifiedPeptidesParams | `None` | Some params runs the modification step. None skips it. | 
| AggregateToProteinGroupsParams | `{ Transform = None; SingleFilters = None; GroupFilters = None; Aggregation = { Light = Mean } }` | The protein step, always run.| 

An AggregationParams record describes one step.
|Parameter |Default | Description | 
| ... | ... | ... | 
| Transform  | `None` | Some { Light = Some Log2 } applies a NumericTransform (Log2, Add, Substract, MultiplyBy, DivideBy) to every intensity before filtering and aggregation. |
| SingleFilters |  `None`  | Some { Light = Some (seq [IsBiggerThan 4.0]) } keeps only intensities that pass every NumericFilter (IsBiggerThan, IsSmallerThan). | 
| GroupFilters | `None` | Some { Light = Some (seq [Tukey 1.5]) } removes outliers within each group before aggregation. GroupFilter is Tukey factor, Stdev factor or TopX count. | 
|Aggregation | `{ Light = Mean }`| How the remaining values of a group are combined. NumericAggregation is Mean, Median or Sum.


**Four default files exists** 
1. [LablefreeQuantification Params](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams.json) only runs the protein step.

2. [LabelFreeQuantificationParams_ChargeAgg.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_ChargeAgg.json) adds the chrage step with Mean. 

3. [LabelFreeQuantificationParams_ChargeAgg_ModAgg.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_ChargeAgg_ModAgg.json) adds the charge and the modification step, both with Mean. 

4. [LabelFreeQuantificationParams_Transform_Filter_Sum.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_Transform_Filter_Sum.json) runs only the protein step, but log2 transforms the intensities, keeps those above 4.0 and sums them. 

## Running the tool 
1. **Installing**

   Install the Tool with the following command inside the choosen Command-Line

   ```sh
   dotnet tool install --global ProteomIQon.LabelFreeProteinQuantification
   ```

2. **Several runs in one table, from a list or from a directory**


     ```sh
      proteomiqon-labelfreeproteinquantification -i path/to/run1.quantAndProt path/to/run2.quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
      proteomiqon-labelfreeproteinquantification -i path/to/quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
     ```

3. **All Flags**

   ```sh
    proteomiqon-labelfreeproteinquantification --help
     ```












**Modification filter options:** `UseModifiedPeptides.All` (use all), `UseModifiedPeptides.No` (exclude modified peptides), or `UseModifiedPeptides.UseOnly` (allow specified modifications).

**About the defaults:** The documentation lists the optional stage parameters as `Some AggregationParams`, but the repository also supplies a *minimal* JSON preset in which both optional fields are `null`. Which stages run therefore depends on the chosen preset.

### Charge-state aggregation

When enabled, this stage combines versions of the same modified peptide that differ only in precursor charge. Set `AggregatePeptideChargeStatesParams` to a `Some` record containing the settings listed under *Shared aggregation settings*. Set it to `None` / JSON `null` to skip the stage.

In the supplied charge-aggregation JSON preset, filtering and transformation are disabled and light intensities are aggregated using the mean.

### Modification aggregation

When enabled, this stage combines peptide versions with the same amino-acid sequence but different modifications. Set `AggregateModifiedPeptidesParams` to a `Some` record using the shared aggregation settings; `None` / JSON `null` skips the stage.

In the supplied charge-and-modification JSON preset, light intensities are aggregated using the mean and the three optional transformation/filter fields are disabled.

### Protein-group aggregation

The final aggregation combines the peptide intensities assigned to each protein group. Unlike the previous two stages, `AggregateToProteinGroupsParams` is not optional and always contains an aggregation-parameter record.

The official JSON presets use the mean of light intensities without additional transformations or filtering for this step.

### Shared aggregation settings

Each enabled stage uses the same four-field `LabelFreeQuantification.AggregationParams` record.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `Transform` | `None` | Optional transformation applied to intensity values. |
| `SingleFilters` | `None` | Optional lower and/or upper limits for individual intensities. |
| `GroupFilters` | `None` | Optional filtering of the group of intensities being aggregated. |
| `Aggregation` | `{Light = NumericAggregation.Mean}` | Specifies how the intensity values in a group are combined. |

**Numeric aggregation choices:** `NumericAggregation.Mean`, `NumericAggregation.Median` and `NumericAggregation.Sum`.

In JSON, the `None` values are represented by `null`, and the mean aggregator is encoded as `"Aggregation": {"Light": {"Case": "Mean"}}`.

**Default parameter files:** Choose the preset matching the aggregation stages you want to run.

- [Minimal: protein-group aggregation only ↗](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams.json)
- [Charge-state + protein-group aggregation ↗](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_ChargeAgg.json)
- [Charge-state + modification + protein-group aggregation ↗](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/LabelFreeQuantificationParams_ChargeAgg_ModAgg.json)

## Outputs

| Stage | Output content | When produced |
| --- | --- | --- |
| Charge-state aggregation | Intensities combined across charge states of the same modified peptide. | Only if `AggregatePeptideChargeStatesParams` is enabled. |
| Modification aggregation | Intensities combined across differently modified versions of the same peptide sequence. | Only if `AggregateModifiedPeptidesParams` is enabled. |
| Protein-group aggregation | Peptide intensities aggregated to protein-group level. | Always performed. |

**Output format:** tab-separated files. The referenced documentation does not specify exact output filenames or a column-by-column output schema, so no unverified filenames or columns are assumed here.

## Installing and executing the tool

1. **Installing**

   Install the Tool with the following command inside the choosen Command-Line

   ```sh
   dotnet tool install --global ProteomIQon.LabelFreeProteinQuantification
   ```

2. **Execute**

   - Running the tool

     ```sh
     proteomiqon-labelfreeproteinquantification -i path/to/run.quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
     ```

3. **Several runs in one table, from a list or from a directory**

   ```sh
   proteomiqon-labelfreeproteinquantification -i path/to/run1.quantAndProt path/to/run2.quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
   proteomiqon-labelfreeproteinquantification -i path/to/quantAndProt -o path/to/output -p path/to/LabelFreeQuantificationParams.json
   ```

4. **Two directories, paired by file name**

   ```sh
   proteomiqon-psmbasedquantification -i path/to/mzliteFolder -ii path/to/qpsmFolder -d path/to/database.db -o path/to/output -p path/to/QuantificationParams.json -c 3 -mf
   ```

5. **All Flags**

   ```sh
   proteomiqon-labelfreeproteinquantification --help
   ```
