(**
---
title: JoinQuantPepIonsWithProteins
category: Tools
categoryindex: 1
index: 8
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