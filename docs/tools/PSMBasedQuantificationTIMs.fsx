(**
---
title: PSMBasedQuantificationTIMs
category: Tools
categoryindex: 1
index: 9
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

# PSMBasedQuantificationTIMs

PSMBasedQuantificationTIMs is the [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html) variant for timsTOF runs. On these instruments a peptide ion is a peak in two dimensions, retention time and ion mobility, so the tool extracts the signal around every identification in both dimensions.

An XIC point is the sum of the intensity within 0.05 of the identification's ion mobility (`IonMobilityWindow`).

The points form a retention time [ion chromatogram](https://www.biofsharp.com/BioFSharp.Mz/03_01_quantification.html#Building-an-XIC-trace). As PSMBasedQuantification does, the tool detects the peaks on the trace and [fits and integrates](https://www.biofsharp.com/BioFSharp.Mz/03_01_quantification.html#Fitting-and-integrating-the-main-peak) the one that contains the identification's retention time. The output has the layout of PSMBasedQuantification.

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | one or more `.mzlite` files with ion mobility, or a directory | [MzMLToMzLiteIonMobility]({{root}}tools/MzMLToMzLiteIonMobility.html) |
| `-ii` | one or more `.qpsm` files, or a directory that is searched for `*.qpsm` | [PSMStatistics]({{root}}tools/PSMStatistics.html) or [MsFraggerToPSM]({{root}}tools/MsFraggerToPSM.html) |
| `-d` | the SQLite peptide database | [PeptideDB]({{root}}tools/PeptideDB.html) |
| `-o` | the output directory, created when missing | |
| `-p` | the parameter file in JSON | this page |
| `-mf` | switch: pair files by base name instead of by position | |
| `-dc` | switch: write diagnostic charts | |
| `-z` | switch: pack the chart directory into a zip file and delete it | |
| `-c` | number of runs quantified at the same time, default 1 | |

The tool refuses any instrument file whose name does not end in `.mzlite`. It reads the qpsm files in the column layout of PSMStatistics, which MsFraggerToPSM writes too, and stops with the missing columns when a file has another layout. The qpsm carries no ion mobility, so the tool reads the ion mobility of every identification from its MS2 spectrum in the mzlite, the spectrum whose id is the `PSMId`. The tool warns about identifications whose spectrum is missing or has no ion mobility. These identifications do not locate the trace. An ion without any identification with ion mobility gets no row, and the tool stops when no identification of the file has an ion mobility.
With one file per flag the pair is used as given. With lists, the n-th mzlite file is paired with the n-th qpsm file unless you pass `-mf`, which pairs files by their base name and skips mzlite files without a partner.

For every run the tool writes `run.quant`, tab separated with a header and the columns described on the [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html) page.
[QuantBasedAlignment]({{root}}tools/QuantBasedAlignment.html), [AlignmentBasedQuantification]({{root}}tools/AlignmentBasedQuantification.html), [AlignmentBasedQuantStatistics]({{root}}tools/AlignmentBasedQuantStatistics.html) and [JoinQuantPepIonsWithProteins]({{root}}tools/JoinQuantPepIonsWithProteins.html) read `.quant` files.

With `-dc` the tool writes charts into `run_plots`, and `-z` packs that directory into `run_plots.zip` and deletes it. Delete an existing `run.quant` in the output directory before a rerun, the tool appends to it. Logs go to `PSMBasedQuantification_log.txt` and `run_log.txt` in the output directory.

## Parameters

The parameter record is the one of PSMBasedQuantification. The defaults below are the timsTOF defaults.

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `PerformLabeledQuantification` | `Labeling.Unlabeled` | Label free quantification of the identified ion. The three labeled cases work as in PSMBasedQuantification and quantify the partner ion too. |
| `XicExtraction.ScanTimeWindow` | `2.0` | Half width of the retention time range around the identification that is extracted, in the scan time unit of the run. |
| `XicExtraction.MzWindow_Da` | `Window.Estimate` | Half width of the m/z window, in Da. The m/z correction and `Window.Estimate` work as in PSMBasedQuantification. `Window.Fixed 0.07` sets 0.07 Da. |
| `XicExtraction.XicProcessing` | `XicProcessing.Wavelet waveletParams` | Peak detection on the retention time trace, summed over the ion mobility window. The wavelet parameters are listed below. `SecondDerivative` is the other supported case. `Gabor3D` is in development and does not work currently. |
| `XicExtraction.TopKPSMs` | `None` | `Some k` keeps only the k identifications with the highest `ModelScore` per peptide ion when placing the peak. |
| `XicExtraction.IonMobilityWindow` | `Window.Fixed 0.05` | Half width of the ion mobility window around the weighted ion mobility of the identifications that centre the trace (`CentrePSMs`), in 1/K0. `Window.Estimate` is not available here and stops the tool. |
| `BaseLineCorrection` | `Some { MaxIterations = 10; Lambda = 6; P = 0.05 }` | Asymmetric least squares baseline subtraction before peak detection. `None` skips it. |
| `QualityFilterLog2Limit` | `2.0` | The quality filter keeps rows with -limit < log2 ((area / median area) / (apex / median apex)) < limit. Rows with a non-finite or non-positive area, apex or ratio are dropped. |
| `CentrePSMs` | `3` | Number of identifications per peptide ion with an ion mobility, ranked by the summed intensity of their matched fragment ions, whose intensity weighted scan time and ion mobility centre the trace. It is applied after `TopKPSMs`, which ranks by `ModelScore`. At least 1. |

The wavelet parameters:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `Borderpadding` | `None` | Number of points padded at both ends of the trace. `None` derives it from the data. |
| `BorderPadMethod` | `Padding.BorderPaddingMethod.Random` | How the border points are filled (`Random` or `Zero`). |
| `InternalPaddingMethod` | `Padding.InternalPaddingMethod.LinearInterpolation` | How small gaps inside the trace are filled. |
| `HugeGapPaddingMethod` | `Padding.HugeGapPaddingMethod.Zero` | How gaps longer than `HugeGapPaddingDistance` are filled. |
| `HugeGapPaddingDistance` | `100.` | Gap length from which a gap counts as huge. |
| `MinPeakDistance` | `Some 0.1` | Peaks closer than this (in minutes) are merged. The distance is at least 1.5 times the median point spacing of the trace, which is also what `None` uses. |
| `MinPeakLength` | `Some 0.1` | Smallest peak width the wavelet scales start at. |
| `MaxPeakLength` | `1.5` | Largest peak width the wavelet scales go up to. |
| `NoiseQuantile` | `0.01` | Quantile of the wavelet correlation used as noise level. |
| `MinSNR` | `0.01` | A peak needs a correlation above `MinSNR` times the noise level. |

The default file is [QuantificationTIMsParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/QuantificationTIMsParams.json).

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
        MinPeakDistance        = Some 0.1
        MinPeakLength          = Some 0.1
        MaxPeakLength          = 1.5
        NoiseQuantile          = 0.01
        MinSNR                 = 0.01
    }

let quantificationTIMsParams : Dto.QuantificationParams =
    {
        PerformLabeledQuantification = Labeling.Unlabeled
        XicExtraction =
            {
                ScanTimeWindow = 2.
                MzWindow_Da    = Window.Estimate
                XicProcessing  = XicProcessing.Wavelet waveletParams
                TopKPSMs       = None
                IonMobilityWindow = Window.Fixed 0.05
            }
        BaseLineCorrection = Some { MaxIterations = 10; Lambda = 6; P = 0.05 }
        QualityFilterLog2Limit = 2.
        CentrePSMs = 3
    }

// Replace the temp folder with your project folder.
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "QuantificationTIMsParams.json")

Json.serializeAndWrite outputPath quantificationTIMsParams

(**
## Running the tool

Install with `dotnet tool install --global ProteomIQon.PSMBasedQuantificationTIMs`, then quantify one run:

```text
proteomiqon-psmbasedquantificationtims -i path/to/run.mzlite -ii path/to/run.qpsm -d path/to/database.db -o path/to/output -p path/to/QuantificationTIMsParams.json
```

Several runs, paired by position in the two lists, three in parallel:

```text
proteomiqon-psmbasedquantificationtims -i path/to/run1.mzlite path/to/run2.mzlite path/to/run3.mzlite -ii path/to/run1.qpsm path/to/run2.qpsm path/to/run3.qpsm -d path/to/database.db -o path/to/output -p path/to/QuantificationTIMsParams.json -c 3
```

Two directories, paired by file name:

```text
proteomiqon-psmbasedquantificationtims -i path/to/mzliteFolder -ii path/to/qpsmFolder -d path/to/database.db -o path/to/output -p path/to/QuantificationTIMsParams.json -c 3 -mf
```

All flags:

```text
proteomiqon-psmbasedquantificationtims --help
```
*)
