(**
// can't yet format YamlFrontmatter (["title: AddDeducedPeptides"; "category: Tools"; "categoryindex: 1"; "index: 10"], Some { StartLine = 2 StartColumn = 0 EndLine = 6 EndColumn = 9 }) to pynb markdown

[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/AddDeducedPeptides.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/AddDeducedPeptides.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/AddDeducedPeptides.ipynb)

# Add Deduced Peptides
**Disclaimer** this tool needs [protein inference](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.html) and [aligned quant](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html)
files.

After the alignment based quantification, the protein inference files do no longer match the peptides present in the quantification files. The quantification file now contains additional quantified peptides based on 
information from other quantification files in the same run. If the protein inference was performed on the combined files earlier, then all inference information is contained in the combined protein inference files. This tool 
takes this combined inference information and assigns it to the peptides present in the new quantification file, thereby creating a new protein inference result for each quantification.

## Executing the Tool

	proteomiqon-adddeducedpeptides -i "path/to/your/run.quant" -ii "path/to/your/prot" -o "path/to/your/outDirectory"

It is also possible to call the tool on a lists of quant and prot files.

	proteomiqon-adddeducedpeptides -i "path/to/your/run1.quant" "path/to/your/run2.quant" "path/to/your/run3.quant" -ii "path/to/your/run1.prot" "path/to/your/run2.prot" "path/to/your/run3.prot" -o "path/to/your/outDirectory"

Alternatively, it is possible to call the tool on folders of quant and prot files.

	proteomiqon-adddeducedpeptides -i "path/to/your/quant" -ii "path/to/your/prot" -o "path/to/your/outDirectory"

A detailed description of the CLI arguments the tool expects can be obtained by calling the tool:

	proteomiqon-adddeducedpeptides --help

*)