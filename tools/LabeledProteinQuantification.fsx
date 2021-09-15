(**
// can't yet format YamlFrontmatter (["title: LabeledProteinQuantification"; "category: Tools"; "categoryindex: 1"; "index: 10"], Some { StartLine = 2 StartColumn = 0 EndLine = 6 EndColumn = 9 }) to pynb markdown

[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/LabeledProteinQuantification.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/LabeledProteinQuantification.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/LabeledProteinQuantification.ipynb)

# Labeled Protein Quantification
**Disclaimer** this tool needs a [quantAndProt](https://csbiology.github.io/ProteomIQon/tools/JoinQuantPepIonsWithProteins.html) file, which combines the results from [ProteinInference](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.html) 
and [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html).

After quantification and protein inference are performed, it is known which peptide originated from which protein, as well as the intensity of each peptide. The information available for each peptide now needs to be aggragated 
for their proteins. 

This tool performs the aggregation from the peptides to the protein in several steps. The first step for the labeled protein quantification is the aggregation of the differently labeled peptides. Peptides with the same sequence, modifications and 
charge are aggregated and the ratio between the intensity from the light and heavy version is calculated. The next two aggregation steps are optional. One of them is the aggregation based on charge state. Similarily to the first step, peptides with 
the same sequence and modifications, but different charge states are being aggregated. The next optional step does the same for peptides with the same sequence, but different modification. Those steps build upon each other. The last step is the aggregation of 
all peptides of a protein. The result of each aggregation step is given as a tab separated file. The aggregation is performed according to the given parameters for each step. If an optional aggregation is not performed, the next step takes the result from the prior aggregation. For example, if aggregation by charge and 
modification are skipped, the protein aggregation gets a collection of peptides, where a peptidesequence can occur with different charge states and midifications.

## Parameters

### LabeledQuantificationParams

The following table gives an overview of the parameter set:

| **Parameter**                       | **Default Value**                                   | **Description**                                                                             |
|-------------------------------------|-----------------------------------------------------|---------------------------------------------------------------------------------------------|
| Correlation\_Light\_Heavy_Threshold | Some 0.0                                            | Optional minimum value for the correlation between light and heavy peptide                  |
| ModificationFilter                  | UseModifiedPeptides.All                             | Specifies which modifications are used during the aggregation steps                         |
| AggregateGlobalModificationsParams  | LabeledProteinQuantification.AggregationParams      | Specifies how the differently labeled versions of a peptide are aggregated                  |
| AggregatePeptideChargeStatesParams  | Some LabeledProteinQuantification.AggregationParams | Specifies how the differently charged versions of a peptide are aggregated (optional step)  |
| AggregateModifiedPeptidesParams     | Some LabeledProteinQuantification.AggregationParams | Specifies how the differently modified versions of a peptide are aggregated (optional step) |
| AggregateToProteinGroupsParams      | LabeledProteinQuantification.AggregationParams      | Specifies how the peptides of a protein are aggregated                                      |

### LabeledProteinQuantification.AggregationParams

The following table gives an overview of the parameter set:

| **Parameter**        | **Default Value**                                                                              | **Description**                                                                               |
|----------------------|------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| LabeledTransform     | None                                                                                           | Possibility to add different transformations to the light or heavy intensities or their ratio |
| LabeledSingleFilters | None                                                                                           | Possibility to add minimal and/or maximal values for the individual intensity or ratio        |
| LabeledGroupFilters  | None                                                                                           | Possibility to apply a filtering step to the intensities or ratios that are being aggregated  |
| LabeledAggregation   | {Light= NumericAggregation.Mean; Heavy=NumericAggregation.Mean; Ratio=NumericAggregation.Mean} | Specifies how the light and heavy intensities and the ratio are aggregated                    |


*)