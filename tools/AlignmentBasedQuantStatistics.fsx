(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/AlignmentBasedQuantStatistics.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantStatistics.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantStatistics.ipynb)

# AlignmentBasedQuantStatistics

A peptide ion transferred by [alignment](https://csbiology.github.io/ProteomIQon/tools/QuantBasedAlignment.html) and quantified by
[AlignmentBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html) has no identification behind it, so nothing says yet whether
the peak found at the predicted scan time is the peptide. AlignmentBasedQuantStatistics gives every transferred quantification an alignment
score and an alignment q-value, so downstream tools can filter transferred quantifications the way identifications are filtered by q-value.
It learns from peptide ions where the truth is known: the .align file also predicts ions the run had identified itself, and for those the
transferred quantification can be compared with the quantification of the identified ion. Ions where m/z and intensity agree closely are
positive training examples, ions where they disagree are negative ones. A gradient boosted tree classifier trained on them scores every
transferred ion, and the q-value follows from the score distribution of the negative examples.

## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | .quant of the run, before alignment | [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html)
`-ii` | .align of the run | [QuantBasedAlignment](https://csbiology.github.io/ProteomIQon/tools/QuantBasedAlignment.html)
`-iii` | .quant of the run, after alignment | [AlignmentBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html)
`-l` | .quant files of the runs used for training | [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html)
`-ll` | .align files of the runs used for training | [QuantBasedAlignment](https://csbiology.github.io/ProteomIQon/tools/QuantBasedAlignment.html)
`-lll` | .quant files after alignment for the runs used for training | [AlignmentBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html)
`-o` | output directory, created when missing | &#32;
`-p` | parameter file (JSON) | this page
`-c` | number of runs scored in parallel, default 1 | &#32;
`-ts` | switch: write the training examples selected from each run | &#32;
`-dc` | switch: write diagnostic charts | &#32;
`-mf` | switch: pair input files by name instead of position | &#32;


`-i` to `-iii` name the runs to score, `-l` to `-lll` the runs to learn from. The runs to learn from can be the same runs. Every one
of the six flags takes one file, several files, or a directory that is searched for *.quant or *.align. Without `-mf` the quant, align and
aligned quant files are paired by position. With `-mf` they are paired by file name without extension. Every scored run trains a classifier on
the pooled examples of all training runs.

The tool writes a file with the same name as the `-i` .quant file into the output directory. It holds the rows of that file with
AlignmentScore and AlignmentQValue filled for the rows whose QuantificationSource is Alignment.
[JoinQuantPepIonsWithProteins](https://csbiology.github.io/ProteomIQon/tools/JoinQuantPepIonsWithProteins.html) reads it, and
[LabelFreeProteinQuantification](https://csbiology.github.io/ProteomIQon/tools/LabelFreeProteinQuantification.html) and
[LabeledProteinQuantification](https://csbiology.github.io/ProteomIQon/tools/LabeledProteinQuantification.html) filter on the q-value column through their Alignment_QValue
parameter. With `-ts` the tool also writes `<quant file name>testset`, the training examples taken from the run of the same name with
their aligned quantification columns, AlignmentScore, AlignmentQValue and the label column WasPositiveSet. With `-dc` it writes ProbabilityHistogram.html and QValueDistribution.html into a directory
named after the run. Logs go to `AlignmentBasedQuantStatistics_log.txt` and one `<run>_log.txt` per scored run in the output directory.

## Parameters

The four cutoffs decide which training ions count as positive and which as negative. Each compares a transferred value with the value of the
identified ion (light or heavy, matching the ion's GlobalMod) and expresses the allowed difference as a fraction of the identified value.

Parameter | Default | Meaning
--- | --- | ---
PositiveQuantMzCutoff | 0.01 | Positive when the m/z difference is below this fraction of the identified m/z and the intensity condition holds too.
NegativeQuantMzCutoff | 0.99 | Negative when the m/z difference is above (1 - NegativeQuantMzCutoff) times the identified m/z, so above 1 percent with the default.
PositiveQuantCutoff | 0.1 | Positive when the intensity difference is below this fraction of the identified intensity and the m/z condition holds too.
NegativeQuantCutoff | 0.9 | Negative when the intensity difference is above (1 - NegativeQuantCutoff) times the identified intensity, so above 10 percent with the default.


The positive test runs first. An ion that fails it is negative, whether or not it passes one of the negative conditions. The default file is
[AlignmentBasedQuantStatisticsParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/AlignmentBasedQuantStatisticsParams.json).

## Writing a parameter file

*)
open ProteomIQon

let alignmentBasedQuantStatisticsParams : Dto.AlignmentBasedQuantStatisticsParams =
    {
        PositiveQuantMzCutoff = 0.01
        NegativeQuantMzCutoff = 0.99
        PositiveQuantCutoff   = 0.1
        NegativeQuantCutoff   = 0.9
    }

// Replace the temp folder with your project folder.
let outputPath =
    System.IO.Path.Combine(System.IO.Path.GetTempPath(), "AlignmentBasedQuantStatisticsParams.json")

Json.serializeAndWrite outputPath alignmentBasedQuantStatisticsParams
(**
## Running the tool

Install the tool with `dotnet tool install --global ProteomIQon.AlignmentBasedQuantStatistics`. Score one run and learn from the same run.

```text
proteomiqon-alignmentbasedquantstatistics -i path/to/run.quant -ii path/to/run.align -iii path/to/aligned/run.quant -l path/to/run.quant -ll path/to/run.align -lll path/to/aligned/run.quant -o path/to/output -p path/to/AlignmentBasedQuantStatisticsParams.json
```

Score every run in the directories, learn from all of them, pair the files by name, and work on four runs at a time.

```text
proteomiqon-alignmentbasedquantstatistics -i path/to/quant -ii path/to/align -iii path/to/aligned -l path/to/quant -ll path/to/align -lll path/to/aligned -o path/to/output -p path/to/AlignmentBasedQuantStatisticsParams.json -mf -c 4
```

Add `-ts` to keep the training examples and `-dc` to get the charts.

```text
proteomiqon-alignmentbasedquantstatistics -i path/to/quant -ii path/to/align -iii path/to/aligned -l path/to/quant -ll path/to/align -lll path/to/aligned -o path/to/output -p path/to/AlignmentBasedQuantStatisticsParams.json -mf -ts -dc
```

Print the description of every argument.

```text
proteomiqon-alignmentbasedquantstatistics --help
```

*)