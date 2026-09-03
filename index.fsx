(**
[![Made with F#](https://img.shields.io/badge/Made%20with-FSharp-rgb(184,69,252).svg)](https://fsharp.org/)
![GitHub contributors](https://img.shields.io/github/contributors/CSBiology/ProteomIQon)

## The ProteomIQon

ProteomIQon is a set of command line tools for the analysis of mass spectrometry based proteomics data, written in F#.
Each tool does one step and writes plain files that the next tool reads, so a pipeline is a sequence of tool calls
that can be rerun from any point. Most tools append to an output file that already exists, so delete the old
file before rerunning a step into the same directory. The chain covers signal detection, peptide identification,
quantification, retention time alignment between runs and protein inference.

<img src="https://csbiology.github.io/ProteomIQon/img/PillarsOfCompProt.png" width="1000" height="750" />
The proteomics theory behind the steps is documented in [BioFSharp.Mz](https://www.biofsharp.com/BioFSharp.Mz/), the
library most tools build on. The tool pages here link to it where a concept needs more than a sentence.

## Installing a tool

Every tool is a .NET tool on [nuget.org](https://www.nuget.org/profiles/CSBiology), package `ProteomIQon.<ToolName>`,
command `proteomiqon-<toolname>`. The one exception is QuantBasedAlignment, which ships as
`ProteomIQon.QuantBasedAlignment_win_x64` and `ProteomIQon.QuantBasedAlignment_linux_x64`. The package
`ProteomIQon.QuantBasedAlignment` is the library both wrap and installs no command.

```text
dotnet tool install --global ProteomIQon.PeptideDB
proteomiqon-peptidedb --help
```

Tools that take parameters read them from a JSON file. Each tool page shows the parameters with a script that writes
the file, and the command line calls.

## How the tools chain

A run of a data dependent acquisition experiment goes through these steps. Each tool writes the files in its last
column, and the next tool reads them.

Step | Tool | Reads | Writes
--- | --- | --- | ---
Convert raw data | [MzMLToMzLite](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLite.html) | .mzML (from msconvert) | .mzlite
Build the search space | [PeptideDB](https://csbiology.github.io/ProteomIQon/tools/PeptideDB.html) | FASTA | .db (SQLite)
Identify spectra | [PeptideSpectrumMatching](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.html) | .mzlite, .db | .psm
Control the FDR | [PSMStatistics](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html) | .psm, .db | .qpsm
Quantify identified peptide ions | [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html) | .mzlite, .qpsm, .db | .quant
Infer proteins | [ProteinInference](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.html) | .qpsm, .db | .prot
Join peptides and proteins | [JoinQuantPepIonsWithProteins](https://csbiology.github.io/ProteomIQon/tools/JoinQuantPepIonsWithProteins.html) | .quant, .prot | .quantAndProt
Aggregate to proteins | [LabelFreeProteinQuantification](https://csbiology.github.io/ProteomIQon/tools/LabelFreeProteinQuantification.html) or [LabeledProteinQuantification](https://csbiology.github.io/ProteomIQon/tools/LabeledProteinQuantification.html) | .quantAndProt | LabelFreeQuant.txt or LabeledQuant.txt


With several runs, the alignment tools transfer identifications between runs before the join. They sit between
PSMBasedQuantification and JoinQuantPepIonsWithProteins:

Step | Tool | Reads | Writes
--- | --- | --- | ---
Map scan times between runs | [QuantBasedAlignment](https://csbiology.github.io/ProteomIQon/tools/QuantBasedAlignment.html) | .quant of all runs | .align, .alignmetric
Quantify transferred peptide ions | [AlignmentBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html) | .mzlite, .align, .alignmetric, .quant, .db | .quant
Score the transfers | [AlignmentBasedQuantStatistics](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantStatistics.html) | .quant, .align | .quant
Reassign protein groups | [AddDeducedPeptides](https://csbiology.github.io/ProteomIQon/tools/AddDeducedPeptides.html) | .quant, .prot | .prot


Four more tools cover other inputs. For 15N labeled samples, [RatioLFQ](https://csbiology.github.io/ProteomIQon/tools/RatioLFQ.html) turns the
light to heavy ratios in LabeledQuant.txt into one intensity per run and protein. For timsTOF data with ion mobility,
[MzMLToMzLiteIonMobility](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLiteIonMobility.html) converts the mzML,
[MsFraggerToPSM](https://csbiology.github.io/ProteomIQon/tools/MsFraggerToPSM.html) imports an MSFragger search, and
[PSMBasedQuantificationTIMs](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantificationTIMs.html) quantifies in retention time and ion
mobility.

## The core project

The ProteomIQon core is referenced by all tools. It contains mainly serializable data transfer objects such as tool results and tool parameters, as well as their mapping to domain specific types.
This is also the place for any kind of code reusable across tools such as thin wrappers around data readers, logging or CLI formatting.

## Documentation

The documentation is generated with [fsdocs](https://fsprojects.github.io/FSharp.Formatting/) from the .fsx files in the docs folder.
If you find a typo, please submit a pull request. [How to document your work](https://csbiology.github.io/ProteomIQon/develop/How_to_document_your_work.html) explains the setup.

## Contributing

Please refer to the CSB [Contribution guidelines](https://github.com/CSBiology/BioFSharp/blob/developer/.github/CONTRIBUTING.md)

## Community/Social

Want to get in touch with us? We recently joined the twitter crowd:

[![Twitter Follow](https://img.shields.io/twitter/follow/cs_biology.svg?style=social)](https://twitter.com/cs_biology)

## Citation

When using ProteomIQon in scientific or commercial releases, please cite us using the DOI `10.5281/zenodo.6335068`.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6335068.svg)](https://doi.org/10.5281/zenodo.6335068)

*)