(**
---
title: AlignmentBasedQuantification
category: Tools
categoryindex: 1
index: 11
---
*)

(*** hide ***)

(*** condition: prepare ***)
#r "nuget: FSharp.Stats, 0.6.0"
#r "nuget: Newtonsoft.Json, 13.0.4"
#r "../../src/ProteomIQon/bin/Release/net10.0/ProteomIQon.dll"

(*** condition: ipynb ***)
#if IPYNB
#r "nuget: ProteomIQon, {{fsdocs-package-version}}"
#endif // IPYNB

(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# AlignmentBasedQuantification

[QuantBasedAlignment]({{root}}tools/QuantBasedAlignment.html) predicts where a peptide ion identified in another run should elute in this run.
AlignmentBasedQuantification goes to that place. For every peptide ion in the .align file it extracts an ion chromatogram from the .mzlite file
around the predicted scan time, and, with PerformLocalWarp set, moves the scan time estimate by dynamic time warping of the extracted chromatogram
against the chromatogram the ion had in the source run. Baseline correction, peak detection, peak fitting and the labeled counterpart work as in
[PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html). The result is a .quant file that holds the quantifications of this run's own
identifications together with the transferred ones.

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | .mzlite of the run | [MzMLToMzLite]({{root}}tools/MzMLToMzLite.html) |
| `-ii` | .align of the run | [QuantBasedAlignment]({{root}}tools/QuantBasedAlignment.html) |
| `-iii` | .alignmetric of the run | [QuantBasedAlignment]({{root}}tools/QuantBasedAlignment.html) |
| `-iv` | .quant of the same run | [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html) |
| `-d` | peptide database (SQLite) | [PeptideDB]({{root}}tools/PeptideDB.html) |
| `-o` | output directory, created when missing | |
| `-p` | parameter file (JSON) | this page |
| `-mf` | switch: pair the four file lists by base name instead of by position | |
| `-dc` | switch: write chromatogram and m/z correction charts | |
| `-c` | number of runs processed in parallel, default 1, in the directory case | |

Each of the four input flags takes exactly one path. When `-i` is an existing file, the tool processes that single run and expects single files for
`-ii` through -iv as well. When `-i` and `-ii` are directories, the tool searches `-i` for *.mzlite, `-ii` for *.align, `-iii` for *.alignmetric and -iv for
*.quant and processes every run it can pair. Without `-mf` the four file lists are paired by position. With `-mf` they are paired by file name
without extension, and a run that lacks one of the four files is skipped with a note in the log.

The tool writes `<run>.quant` into the output directory, a tab separated table with the rows of the input .quant file followed by one row per
transferred peptide ion. The column QuantificationSource tells them apart (PSM or Alignment). AlignmentScore and AlignmentQValue hold NaN
until [AlignmentBasedQuantStatistics]({{root}}tools/AlignmentBasedQuantStatistics.html) fills them. The file is also read by
[AddDeducedPeptides]({{root}}tools/AddDeducedPeptides.html) and
[JoinQuantPepIonsWithProteins]({{root}}tools/JoinQuantPepIonsWithProteins.html). Next to it the tool writes `<run>.alignquantMetrics`, a
table of the alignment metrics in quantification format.
With `-dc` the chromatograms and the m/z correction go as HTML charts into `<run>_plots`. Delete an existing `<run>.quant` or `<run>.alignquantMetrics` in the
output directory before a rerun, the tool appends to them. Logs go to `AlignmentBasedQuantification_log.txt` and one `<run>_log.txt` per run.

## Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| PerformLabeledQuantification | true | Also extract and quantify the differentially labeled counterpart of every transferred ion. |
| PerformLocalWarp | true | Refine the predicted scan time by dynamic time warping against the source chromatogram before peak detection. |
| XicExtraction | see below | How the ion chromatogram is extracted and processed. |
| BaseLineCorrection | Some { MaxIterations = 10; Lambda = 6; P = 0.05 } | Asymmetric least squares baseline correction of the chromatogram. None switches it off. |

XicExtraction is the same record PSMBasedQuantification uses.

| Field | Default | Meaning |
|-------|---------|---------|
| ScanTimeWindow | 2.0 | Half width in minutes of the scan time range extracted around the predicted scan time. |
| MzWindow_Da | Window.Estimate | Keep Window.Estimate. The tool estimates the m/z tolerance from the run itself. |
| XicProcessing | XicProcessing.Wavelet waveletParams | Peak detection method. The wavelet fields are described on the [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html) page. |
| TopKPSMs | None | Keep None. Transferred ions have no identifications to rank. |

The default file is [AlignmentBasedQuantificationParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/AlignmentBasedQuantificationParams.json).

## Writing a parameter file
*)

open ProteomIQon
open ProteomIQon.Domain
open FSharp.Stats.Signal

let waveletParams : WaveletParameters =
    {
        Borderpadding          = None
        BorderPadMethod        = Padding.BorderPaddingMethod.Random
        InternalPaddingMethod  = Padding.InternalPaddingMethod.LinearInterpolation
        HugeGapPaddingMethod   = Padding.HugeGapPaddingMethod.Zero
        HugeGapPaddingDistance = 100.
        MinPeakDistance        = None
        MinPeakLength          = Some 0.1
        MaxPeakLength          = 1.5
        NoiseQuantile          = 0.01
        MinSNR                 = 0.01
    }

let alignmentBasedQuantificationParams : Dto.AlignmentBasedQuantificationParams =
    {
        PerformLabeledQuantification = true
        PerformLocalWarp             = true
        XicExtraction =
            {
                ScanTimeWindow = 2.
                MzWindow_Da    = Window.Estimate
                XicProcessing  = XicProcessing.Wavelet waveletParams
                TopKPSMs       = None
            }
        BaseLineCorrection = Some { MaxIterations = 10; Lambda = 6; P = 0.05 }
    }

// Replace the temp folder with your project folder.
let outputPath =
    System.IO.Path.Combine(System.IO.Path.GetTempPath(), "AlignmentBasedQuantificationParams.json")

Json.serializeAndWrite outputPath alignmentBasedQuantificationParams

(**
## Running the tool

Install the tool with `dotnet tool install --global ProteomIQon.AlignmentBasedQuantification`. The four run specific inputs share the run name, so the same base name appears four times.

```text
proteomiqon-alignmentbasedquantification -i path/to/run.mzlite -ii path/to/run.align -iii path/to/run.alignmetric -iv path/to/run.quant -d path/to/peptides.db -o path/to/output -p path/to/AlignmentBasedQuantificationParams.json
```

Process every run in a set of directories, paired by position, three at a time.

```text
proteomiqon-alignmentbasedquantification -i path/to/mzlite -ii path/to/align -iii path/to/alignmetric -iv path/to/quant -d path/to/peptides.db -o path/to/output -p path/to/AlignmentBasedQuantificationParams.json -c 3
```

The same, paired by file name.

```text
proteomiqon-alignmentbasedquantification -i path/to/mzlite -ii path/to/align -iii path/to/alignmetric -iv path/to/quant -d path/to/peptides.db -o path/to/output -p path/to/AlignmentBasedQuantificationParams.json -c 3 -mf
```

Print the description of every argument.

```text
proteomiqon-alignmentbasedquantification --help
```
*)
