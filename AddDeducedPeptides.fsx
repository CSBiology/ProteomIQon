(**
---
title: AddDeducedPeptides
category: Tools
categoryindex: 1
index: 7
---
*)
(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# Add Deduced Peptides
**Disclaimer** this tool needs [protein inference]({{root}}tools/proteinInference.html) and [aligned quant]({{root}}tools/psmBasedQuantification.html)
files.

After the alignment based quantification, the protein inference file does no longer match the peptides present in the quantification file. The quantification file now contains additional quantified peptides based on 
information from other quantification files in the same run. If the protein inference was performed on the combined files, then all inference information is still contained in the combined protein inference files. This tool 
takes this combined inference information and assigns it to the peptides present in the new quantification file, thereby creating a new protein inference file.

## Executing the Tool
*)

(**
	proteomiqon-adddeducedpeptides -i "path/to/your/run.quant" -ii "path/to/your/prot" -o "path/to/your/outDirectory"
*)

(**
It is also possible to call the tool on a lists of quant and prot files.
*)

(**
	proteomiqon-adddeducedpeptides -i "path/to/your/run1.quant" "path/to/your/run2.quant" "path/to/your/run3.quant" -ii "path/to/your/run1.prot" "path/to/your/run2.prot" "path/to/your/run3.prot" -o "path/to/your/outDirectory"
*)

(**
Alternatively, it is possible to call the tool on folders of quant and prot files.
*)

(**
	proteomiqon-adddeducedpeptides -i "path/to/your/quant" -ii "path/to/your/prot" -o "path/to/your/outDirectory"
*)

(**
A detailed description of the CLI arguments the tool expects can be obtained by calling the tool:
*)

(**
	proteomiqon-psmbasedquantification --help
*)
