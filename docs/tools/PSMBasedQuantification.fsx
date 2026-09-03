(**
---
title: PSMBasedQuantification
category: Tools
categoryindex: 1
index: 7
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

# PSMBasedQuantification

PSMBasedQuantification turns the identifications from [PSMStatistics]({{root}}tools/PSMStatistics.html) into peptide ion intensities.
For every identified peptide ion (sequence, modification state and charge) it takes the scan times of its MS/MS identifications, extracts an [ion chromatogram](https://www.biofsharp.com/BioFSharp.Mz/03_01_quantification.html#Building-an-XIC-trace) of the monoisotopic m/z in a window around them, detects the peaks in that trace and [fits the peak](https://www.biofsharp.com/BioFSharp.Mz/03_01_quantification.html#Fitting-and-integrating-the-main-peak) closest to the identification. The fitted area is the reported quantity.
For metabolically labeled samples the tool computes the m/z of the partner ion (the heavy form of a light identification, or the other way round) and quantifies that trace as well, so every row can carry a light and a heavy value.

<img src="{{root}}img/LabeledQuant.png" width="1000" height="750" />

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | one or more `.mzlite` files, or a directory | [MzMLToMzLite]({{root}}tools/MzMLToMzLite.html) |
| `-ii` | one or more `.qpsm` files, or a directory that is searched for `*.qpsm` | [PSMStatistics]({{root}}tools/PSMStatistics.html) |
| `-d` | the SQLite peptide database | [PeptideDB]({{root}}tools/PeptideDB.html) |
| `-o` | the output directory, created when missing | |
| `-p` | the parameter file in JSON | this page |
| `-mf` | switch: pair files by base name instead of by position | |
| `-dc` | switch: write diagnostic charts | |
| `-z` | switch: pack the chart directory into a zip file and delete it | |
| `-c` | number of runs quantified at the same time, default 1 | |

The input must be `.mzlite`. A directory given to `-i` is also scanned for `.mzML`, but the tool opens every file with the mzlite reader and fails on anything else, so convert first.
With one file per flag the pair is used as given. With lists, the n-th mzlite file is paired with the n-th qpsm file unless you pass `-mf`, which pairs files by their base name (`run1.mzlite` with `run1.qpsm`) and skips mzlite files without a partner.

For every run the tool writes `run.quant`, tab separated with a header and one row per quantified peptide ion. Next to the identifiers and the best q-value and PEP value of the ion it reports the fitted quantity (`Quant_Light`, `Quant_Heavy`), the apex intensity, the fit parameters, the difference between identification and fitted scan time, the extracted traces and the isotope pattern for both channels.
[QuantBasedAlignment]({{root}}tools/QuantBasedAlignment.html), [AlignmentBasedQuantification]({{root}}tools/AlignmentBasedQuantification.html), [AlignmentBasedQuantStatistics]({{root}}tools/AlignmentBasedQuantStatistics.html) and [JoinQuantPepIonsWithProteins]({{root}}tools/JoinQuantPepIonsWithProteins.html) read `.quant` files.

With `-dc` the tool writes one chart per peptide ion into `run_plots`, plus `mzErrorAndCorrection.html` and `precMzCorrected.html` for the m/z calibration it derives from the identifications. `-z` packs that directory into `run_plots.zip` and deletes it. Delete an existing `run.quant` in the output directory before a rerun, the tool appends to it. Logs go to `PSMBasedQuantification_log.txt` and `run_log.txt` in the output directory.

## Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `PerformLabeledQuantification` | `Labeling.N15Labeling` | `Unlabeled` quantifies only the identified ion. `N15Labeling` quantifies the partner ion too and keeps rows where both fits pass the quality filter. `N15LabelingOnly` applies the filter to the heavy channel alone. `Labelshift` quantifies both and skips the filter. |
| `FragPipe` | `false` | Marks qpsm input in the FragPipe layout, which [PSMBasedQuantificationTIMs]({{root}}tools/PSMBasedQuantificationTIMs.html) reads. Keep it at `false` for this tool. |
| `XicExtraction.ScanTimeWindow` | `2.0` | Half width of the scan time window around the identification that is extracted and searched for the peak, in the scan time unit of the run. |
| `XicExtraction.MzWindow_Da` | `Window.Estimate` | Half width of the m/z window of the trace, in Da. `Estimate` uses four times the MS1 mass error the tool measures on the precursor m/z of the identifications. `Window.Fixed 0.07` sets it to 0.07 Da. |
| `XicExtraction.XicProcessing` | `XicProcessing.Wavelet waveletParams` | Peak detection on the trace. The wavelet parameters are listed below. `SecondDerivative` is the other supported case. `Gabor3D` belongs to the TIMs tool and makes this tool fail. |
| `XicExtraction.TopKPSMs` | `None` | `Some k` keeps only the k identifications with the highest `SequestScore` per peptide ion when placing the scan time. |
| `BaseLineCorrection` | `Some { MaxIterations = 10; Lambda = 6; P = 0.05 }` | Asymmetric least squares baseline subtraction on the trace before peak detection. `None` skips it. |

The wavelet parameters:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `Borderpadding` | `None` | Number of points padded at both ends of the trace. `None` derives it from the data. |
| `BorderPadMethod` | `Padding.BorderPaddingMethod.Random` | How the border points are filled (`Random` or `Zero`). |
| `InternalPaddingMethod` | `Padding.InternalPaddingMethod.LinearInterpolation` | How small gaps inside the trace are filled. |
| `HugeGapPaddingMethod` | `Padding.HugeGapPaddingMethod.Zero` | How gaps longer than `HugeGapPaddingDistance` are filled. |
| `HugeGapPaddingDistance` | `100.` | Gap length from which a gap counts as huge. |
| `MinPeakDistance` | `None` | Peaks closer than this are merged. `None` uses the point spacing of the trace. |
| `MinPeakLength` | `Some 0.1` | Smallest peak width the wavelet scales start at. |
| `MaxPeakLength` | `1.5` | Largest peak width the wavelet scales go up to. |
| `NoiseQuantile` | `0.01` | Quantile of the wavelet correlation used as noise level. |
| `MinSNR` | `0.01` | A peak needs a correlation above `MinSNR` times the noise level. |

The default file is [QuantificationParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/QuantificationParams.json). To run without baseline correction set `BaseLineCorrection = None` in the script below.

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

let quantificationParams : Dto.QuantificationParams =
    {
        PerformLabeledQuantification = Labeling.N15Labeling
        FragPipe                     = false
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
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "QuantificationParams.json")

Json.serializeAndWrite outputPath quantificationParams

(**
## Running the tool

Install with `dotnet tool install --global ProteomIQon.PSMBasedQuantification`, then quantify one run:

```text
proteomiqon-psmbasedquantification -i path/to/run.mzlite -ii path/to/run.qpsm -d path/to/database.db -o path/to/output -p path/to/QuantificationParams.json
```

Several runs, paired by position in the two lists, three in parallel:

```text
proteomiqon-psmbasedquantification -i path/to/run1.mzlite path/to/run2.mzlite path/to/run3.mzlite -ii path/to/run1.qpsm path/to/run2.qpsm path/to/run3.qpsm -d path/to/database.db -o path/to/output -p path/to/QuantificationParams.json -c 3
```

Two directories, paired by file name:

```text
proteomiqon-psmbasedquantification -i path/to/mzliteFolder -ii path/to/qpsmFolder -d path/to/database.db -o path/to/output -p path/to/QuantificationParams.json -c 3 -mf
```

All flags:

```text
proteomiqon-psmbasedquantification --help
```
*)
