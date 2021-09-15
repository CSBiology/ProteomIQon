(**
---
title: LabelFreeProteinQuantification
category: Tools
categoryindex: 1
index: 9
---
*)

(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# Labelfree Protein Quantification

**Disclaimer** this tool needs a [quantAndProt]({{root}}tools/JoinQuantPepIonsWithProteins.html) file, which combines the results from [ProteinInference]({{root}}tools/ProteinInference.html) 
and [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html).

After quantification and protein inference are performed, it is known which peptide originated from which protein, as well as the intensity of each peptide. The information available for each peptide now needs to be aggragated 
for their proteins. 

This tool performs the aggregation from the peptides to the protein in several steps. The first two aggregation steps are optional. One of them is the aggregation based on charge state. Peptides with 
the same sequence and modifications, but different charge states are being aggregated. The next optional step does the same for peptides with the same sequence, but different modifications. Those steps build upon each other. The last step is the aggregation of 
all peptides of a protein. The result of each aggregation step is given as a tab separated file. The aggregation is performed according to the given parameters for each step. If an optional aggregation is not performed, the next step takes the result from the prior aggregation. For example, if aggregation by charge and 
modification are skipped, the protein aggregation is performed on previously unaggregated peptides, where a peptidesequence can occur with different charge states and modifications.

## Parameters

### LabelFreeQuantificationParams

The following table gives an overview of the parameter set:

| **Parameter**                       | **Default Value**                              | **Description**                                                                             |
|-------------------------------------|------------------------------------------------|---------------------------------------------------------------------------------------------|
| ModificationFilter                  | UseModifiedPeptides.All                        | Specifies which modifications are used during the aggregation steps                         |
| AggregatePeptideChargeStatesParams  | Some LabelFreeQuantification.AggregationParams | Specifies how the differently charged versions of a peptide are aggregated (optional step)  |
| AggregateModifiedPeptidesParams     | Some LabelFreeQuantification.AggregationParams | Specifies how the differently modified versions of a peptide are aggregated (optional step) |
| AggregateToProteinGroupsParams      | LabelFreeQuantification.AggregationParams      | Specifies how the peptides of a protein are aggregated                                      |

### LabelFreeQuantification.AggregationParams

The following table gives an overview of the parameter set:

| **Parameter** | **Default Value**                | **Description**                                                                    |
|---------------|----------------------------------|------------------------------------------------------------------------------------|
| Transform     | None                             | Possibility to add different transformations the intensities                       |
| SingleFilters | None                             | Possibility to add minimal and/or maximal values for the individual intensity      |
| GroupFilters  | None                             | Possibility to apply a filtering step to the intensities that are being aggregated |
| Aggregation   | {Light= NumericAggregation.Mean} | Specifies how the intensities are aggregated                                       |

*)