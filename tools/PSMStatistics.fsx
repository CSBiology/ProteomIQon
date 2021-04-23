(**
// can't yet format YamlFrontmatter (["title: PSMStatistics"; "category: Tools"; "categoryindex: 1"; "index: 3"], Some { StartLine = 2 StartColumn = 0 EndLine = 6 EndColumn = 8 }) to pynb markdown

[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/PSMStatistics.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.ipynb)

# Peptide Spectrum Matching Statistics
**Disclaimer** this tool relies on the output of the tools [PeptideDB](https://csbiology.github.io/ProteomIQon/tools/peptideDb.html) and [PeptideSpectrumMatching](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.html).

An established method to identify acquired MS/MS spectra is the [comparison](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.html) of each spectrum with peptides in a [reference database](https://csbiology.github.io/ProteomIQon/tools/peptideDb.html). 

To measure the similarity of in silico generated spectra and measured MS/MS scans we use our own implementations of three established search enginge scores: SEQUEST, Andromeda and XTandem. 
Additionally, we also record quality control parameters such as the mass difference between the precursor ion and the theoretically calulated mass or the uniquness of each score in comparison to 'competing'
peptides within the search space. The PSMStatistics tool utilizes semi supervised machine learning techniques to integrate search engine scores as well as the mentioned quality scores into one single consensus score. 

<img src="https://csbiology.github.io/ProteomIQon/img/SemiSupervisedScoring.png" width="1000" height="750" />
<img src="https://csbiology.github.io/ProteomIQon/img/SemiSupervisedScoring.png" width="1000" height="750" />

Since the search space is extended by so called decoys - reversed counterparts of peptides within the search space - we can estimate the distribution of 'true negatives' and calculate local (PEP values) and global (Q values)
false discovery rates at each consensus score. 
The reported peptides at user defined local and global FDR cutoffs can then be used as inputs for any downstream analysis be it [ProteinInference](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.html) or [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html) 

## Parameters
The following table gives an overview of the parameter set:

| **Parameter**                  | **Default Value**                                                                                                                         | **Description**                                                    |
|--------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------|
| Threshold                      | Threshold.Estimate {QValueThreshold = 0.01; PepValueThreshold = 0.05; MaxIterations=15; MinimumIncreaseBetweenIterations=0.005 }          | Parameters used for FDR based filtering of scored peptides         |
| ParseProteinIDRegexPattern     | id                                                                                                                                        | Regex pattern for parsing of the protein IDs in the database       |
| KeepTemporaryFiles             | false                                                                                                                                     | Indicates if temporary files should be kept or discarded           |

## Parameter Generation

Parameters are handed to the cli tool as a .json file. you can download the default file [here](https://github.com/CSBiology/ProteomIQon/blob/master/src/ProteomIQon/defaultParams/psmstatisticsparams.json), 
or use an F# script, which can be downloaded or run in Binder at the top of the page, to write your own parameter file:

*)
#r "nuget: BioFSharp.Mz, 0.1.5-beta"
#r "nuget: ProteomIQon, 0.0.1"

open ProteomIQon
open ProteomIQon.Domain
open BioFSharp.Mz


let defaultPSMStatistics : Dto.PSMStatisticsParams = 
    {
        Threshold = Threshold.Estimate {QValueThreshold = 0.01; PepValueThreshold = 0.05;MaxIterations=15;MinimumIncreaseBetweenIterations=0.005}
        ParseProteinIDRegexPattern  = "id"
        KeepTemporaryFiles          = true
    }

let serialized = 
    defaultPSMStatistics
    |> Json.serialize
(**
## Executing the Tool
**Disclaimer** this tool relies on the output of the tools [PeptideDB](https://csbiology.github.io/ProteomIQon/tools/peptideDb.html) and [PeptideSpectrumMatching](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.html).

To rescore all MS/MS to identify 'true' psms call: 


	proteomiqon-psmstatistics -i "path/to/your/run.psm" -d "path/to/your/database.sqlite" -o "path/to/your/outDirectory" -p "path/to/your/params.json"

It is also possible to call the tool on a list of .psm files. If you have a mulitcore cpu it is possible to rescore multiple runs in parallel using the -c flag:

	proteomiqon-psmstatistics -i "path/to/your/run1.psm" "path/to/your/run2.psm" "path/to/your/run3.psm" -d "path/to/your/database.sqlite" -o "path/to/your/outDirectory" -p "path/to/your/params.json" -c 3

A detailed description of the CLI arguments the tool expects can be obtained by calling the tool:

	proteomiqon-psmstatistics --help

*)