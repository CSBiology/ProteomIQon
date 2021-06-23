(**
---
title: PSMStatistics
category: Tools
categoryindex: 1
index: 5
---
*)
(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# Peptide Spectrum Matching Statistics
**Disclaimer** this tool relies on the output of the tools [PeptideDB]({{root}}tools/peptideDb.html) and [PeptideSpectrumMatching]({{root}}tools/PeptideSpectrumMatching.html).

An established method to identify acquired MS/MS spectra is the [comparison]({{root}}tools/PeptideSpectrumMatching.html) of each spectrum with peptides in a [reference database]({{root}}tools/peptideDb.html). 

To measure the similarity of in silico generated spectra and measured MS/MS scans we use our own implementations of three established search enginge scores: SEQUEST, Andromeda and XTandem. 
Additionally, we also record quality control parameters such as the mass difference between the precursor ion and the theoretically calulated mass or the uniquness of each score in comparison to 'competing'
peptides within the search space. The PSMStatistics tool utilizes semi supervised machine learning techniques to integrate search engine scores as well as the mentioned quality scores into one single consensus score. 

<img src="{{root}}img/SemiSupervisedScoring.png" width="1000" height="750" />

Since the search space is extended by so called decoys - reversed counterparts of peptides within the search space - we can estimate the distribution of 'true negatives' and calculate local (PEP values) and global (Q values)
false discovery rates at each consensus score. 
The reported peptides at user defined local and global FDR cutoffs can then be used as inputs for any downstream analysis be it [ProteinInference]({{root}}tools/ProteinInference.html) or [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html) 

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

(***condition:ipynb***)
#if IPYNB
(**
If you are running this tool in Binder, you can copy the output of the following codeblock and save it in a JSON file.
*)
serialized
#endif // IPYNB

(**
## Executing the Tool
**Disclaimer** this tool relies on the output of the tools [PeptideDB]({{root}}tools/peptideDb.html) and [PeptideSpectrumMatching]({{root}}tools/PeptideSpectrumMatching.html).

To rescore all MS/MS to identify 'true' psms call: 

*)

(**
	proteomiqon-psmstatistics -i "path/to/your/run.psm" -d "path/to/your/database.sqlite" -o "path/to/your/outDirectory" -p "path/to/your/params.json"
*)

(**
It is also possible to call the tool on a list of .psm files. If you have a mulitcore cpu it is possible to rescore multiple runs in parallel using the -c flag:
*)

(**
	proteomiqon-psmstatistics -i "path/to/your/run1.psm" "path/to/your/run2.psm" "path/to/your/run3.psm" -d "path/to/your/database.sqlite" -o "path/to/your/outDirectory" -p "path/to/your/params.json" -c 3
*)

(**
To create diagnostic plots which show the performance of the psm scorer after a iteration, you can specify the -dc flag:
*)

(**
	proteomiqon-psmstatistics -i "path/to/your/run1.psm" "path/to/your/run2.psm" "path/to/your/run3.psm" -d "path/to/your/database.sqlite" -o "path/to/your/outDirectory" -p "path/to/your/params.json" -c 3 -dc
*)

(**
A detailed description of the CLI arguments the tool expects can be obtained by calling the tool:
*)

(**
	proteomiqon-psmstatistics --help
*)

