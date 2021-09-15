(**
---
title: LabeledProteinQuantification
category: Tools
categoryindex: 1
index: 10
---
*)

(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# Labeled Protein Quantification

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