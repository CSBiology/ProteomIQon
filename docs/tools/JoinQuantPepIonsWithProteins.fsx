(**
---
title: JoinQuantPepIonsWithProteins
category: Tools
categoryindex: 1
index: 9
---
*)

(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# Join Quant Peptide-Ions With Proteins

**Disclaimer** this tool needs the results from [ProteinInference]({{root}}tools/ProteinInference.html) and [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html).

Results from [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html) contain detailed information about the quantification of every peptide, but only basic informations about the protein they originated from. 
This information is obtained during the [ProteinInference]({{root}}tools/ProteinInference.html). This tool takes the results from both tools and combines the qauntification information of each peptide with the more accurate information about 
the corresponding protein including its q-value.

*)
(**
## Executing the Tool

**Disclaimer** this tool needs the results from [ProteinInference]({{root}}tools/ProteinInference.html) and [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html).

To combine all quantified peptides with their infered protein simply call:

*)

(**
	proteomiqon-joinquantpepionswithproteins -i "path/to/your/run/quant" -ii "path/to/your/run/prot" -o "path/to/your/outDirectory"
*)

(**
It is also possible to call the tool on a list of qaunt and prot files:
*)

(**
	proteomiqon-joinquantpepionswithproteins -i "path/to/your/run1.quant" "path/to/your/run2.quant" "path/to/your/run3.quant" -ii "path/to/your/run1.prot" "path/to/your/run2.prot" "path/to/your/run3.prot" -o "path/to/your/outDirectory"
*)

(**
A detailed description of the CLI arguments the tool expects can be obtained by calling the tool:
*)

(**
	proteomiqon-joinquantpepionswithproteins --help
*)