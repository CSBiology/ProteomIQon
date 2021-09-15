(**
// can't yet format YamlFrontmatter (["title: JoinQuantPepIonsWithProteins"; "category: Tools"; "categoryindex: 1"; "index: 8"], Some { StartLine = 2 StartColumn = 0 EndLine = 6 EndColumn = 8 }) to pynb markdown

[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/JoinQuantPepIonsWithProteins.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/JoinQuantPepIonsWithProteins.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/JoinQuantPepIonsWithProteins.ipynb)

# Join Quant Peptide-Ions With Proteins
**Disclaimer** this tool needs the results from [ProteinInference](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.html) and [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html).

Results from [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html) contain detailed information about the quantification of every peptide, but only basic informations about the protein they originated from. 
This information is obtained during the [ProteinInference](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.html). This tool takes the results from both tools and combines the qauntification information of each peptide with the more accurate information about 
the corresponding protein including its q-value.


*)