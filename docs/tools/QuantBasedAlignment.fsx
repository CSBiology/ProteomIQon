(**
---
title: QuantBasedAlignment
category: Tools
categoryindex: 1
index: 7
---
*)
(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# Quantification based Alignment
**Disclaimer** this tool needs [aligned quant]({{root}}tools/PSMBasedQuantification.html) files.

One of the drawbacks of data-dependent acquisition is the stochastic nature of peptide ion selection for MSMS fragmentation as a prerequisite for
peptide identification and quantification. A way to overcome this drawback is the transfer of identified ions from one run to another using the assumption that the run is merely
lacking a successful MSMS scan, but still containing the peptide itself. 
While this ion transfer is commonly referred to as 'run alignment', strictly speaking, the alignment itself is only the process of finding a function f 
that can map from one point of the rt/mz/intensity plane of a source run (parameter -ii) to a target run (parameter -i). 
In this tool, we use smoothing splines as a nonparametric regression technique to estimate the function f for all given source runs (parameter -ii). 
With this estimators at hand, we are able to predict the scan time of individual ions in the target runs. This calculation can then be used as a starting point for
further refinements (e.g. local alignments) or quantification procedures.    

## Executing the Tool
*)

(**
	proteomiqon-quantbasedalignment -i "path/to/your/targetRun1.quant" -ii "path/to/your/sourceRun1.quant" -o "path/to/your/outDirectory" 
*)

(**
It is also possible to call the tool on a lists of target and source files. If you have a mulitcore cpu it is possible to align multiple runs in parallel using the -c flag:
*)

(**
	proteomiqon-quantbasedalignment -i "path/to/your/targetRun1.quant" "path/to/your/targetRun2.quant" -ii "path/to/your/sourceRun1.quant" "path/to/your/sourceRun2.quant" -o "path/to/your/outDirectory" -c
*)

(**
Alternatively, it is possible to call the tool on folders of target and source files.
*)

(**
	proteomiqon-quantbasedalignment -i "path/to/your/targets" -ii "path/to/your/source" -o "path/to/your/outDirectory" 
*)

(**
A detailed description of the CLI arguments the tool expects can be obtained by calling the tool:
*)

(**
	proteomiqon-quantbasedalignment --help
*)
