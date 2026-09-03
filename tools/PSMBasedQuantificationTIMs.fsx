(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/PSMBasedQuantificationTIMs.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantificationTIMs.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantificationTIMs.ipynb)

# PSMBasedQuantificationTIMs

PSMBasedQuantificationTIMs is the [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html) variant for timsTOF runs. On these instruments a peptide ion is a peak in two dimensions, retention time and ion mobility, so the tool extracts an intensity grid over that plane around every identification instead of a one dimensional [ion chromatogram](https://www.biofsharp.com/BioFSharp.Mz/03_01_quantification.html#Building-an-XIC-trace).
By default it locates the peak on the grid with a Gabor filter (a Gaussian envelope multiplied with a cosine wave, tuned in size, width, frequency and angle), then [fits and integrates](https://www.biofsharp.com/BioFSharp.Mz/03_01_quantification.html#Fitting-and-integrating-the-main-peak) it. The output has the same layout as the output of PSMBasedQuantification, so the downstream tools do not care which of the two produced it.

## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | one or more `.mzlite` files with ion mobility, or a directory | [MzMLToMzLiteIonMobility](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLiteIonMobility.html)
`-ii` | one or more `.qpsm` files in the FragPipe column layout, or a directory that is searched for `*.qpsm` | [MsFraggerToPSM](https://csbiology.github.io/ProteomIQon/tools/MsFraggerToPSM.html)
`-d` | the SQLite peptide database | [PeptideDB](https://csbiology.github.io/ProteomIQon/tools/PeptideDB.html)
`-o` | the output directory, created when missing | &#32;
`-p` | the parameter file in JSON | this page
`-mf` | switch: pair files by base name instead of by position | &#32;
`-dc` | switch: write diagnostic charts | &#32;
`-z` | switch: pack the chart directory into a zip file and delete it | &#32;
`-c` | number of runs quantified at the same time, default 1 | &#32;


The tool refuses any instrument file whose name does not end in `.mzlite`. It always reads the qpsm files as the FragPipe layout written by MsFraggerToPSM (with `Hyperscore` and `IonMobility` columns), so the files from [PSMStatistics](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html) do not fit here.
With one file per flag the pair is used as given. With lists, the n-th mzlite file is paired with the n-th qpsm file unless you pass `-mf`, which pairs files by their base name and skips mzlite files without a partner.

For every run the tool writes `run.quant`, tab separated with a header and the columns described on the [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html) page.
[QuantBasedAlignment](https://csbiology.github.io/ProteomIQon/tools/QuantBasedAlignment.html), [AlignmentBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html), [AlignmentBasedQuantStatistics](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantStatistics.html) and [JoinQuantPepIonsWithProteins](https://csbiology.github.io/ProteomIQon/tools/JoinQuantPepIonsWithProteins.html) read `.quant` files.

With `-dc` the tool writes charts into `run_plots`, and `-z` packs that directory into `run_plots.zip` and deletes it. Delete an existing `run.quant` in the output directory before a rerun, the tool appends to it. Logs go to `PSMBasedQuantification_log.txt` and `run_log.txt` in the output directory.

## Parameters

The parameter record is the one of PSMBasedQuantification. The defaults below are the timsTOF defaults.

Parameter | Default | Meaning
--- | --- | ---
`PerformLabeledQuantification` | `Labeling.Unlabeled` | Label free quantification of the identified ion. The three labeled cases work as in PSMBasedQuantification and quantify the partner ion too.
`FragPipe` | `true` | Marks that the qpsm input has the FragPipe layout written by MsFraggerToPSM. Keep it at `true` for this tool.
`XicExtraction.ScanTimeWindow` | `2.0` | Half width of the retention time range around the identification that is extracted, in the scan time unit of the run.
`XicExtraction.MzWindow_Da` | `Window.Estimate` | Half width of the m/z window, in Da. `Estimate` uses four times the MS1 mass error measured on the identifications. `Window.Fixed 0.07` sets it to 0.07 Da.
`XicExtraction.XicProcessing` | `XicProcessing.Gabor3D gaborParams` | Peak detection on the retention time and mobility grid. The Gabor parameters are listed below. `Wavelet` and `SecondDerivative` fall back to the one dimensional processing of PSMBasedQuantification.
`XicExtraction.TopKPSMs` | `None` | `Some k` keeps only the k identifications with the highest `Hyperscore` per peptide ion when placing the peak.
`BaseLineCorrection` | `None` | No baseline subtraction. `Some { MaxIterations = 10; Lambda = 6; P = 0.05 }` enables the asymmetric least squares baseline.


The Gabor parameters. The tool rejects a kernel size at or below zero and any width, frequency or angle that is not finite (widths must also be positive).

Parameter | Default | Meaning
--- | --- | ---
`sizeX` | `41` | Kernel size along retention time, in grid bins. Half of it is the kernel radius, which widens the search window around the expected position.
`sizeY` | `41` | Kernel size along ion mobility, in grid bins.
`sigmaX` | `8.0` | Width of the Gaussian envelope along the first kernel axis, in bins.
`sigmaY` | `8.0` | Width of the Gaussian envelope along the second kernel axis, in bins.
`frequency` | `2.0` | Frequency of the cosine wave inside the kernel.
`theta` | `0.06` | Rotation of the kernel in the plane, in radians.


The default file is [QuantificationTIMsParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/QuantificationTIMsParams.json).

## Writing a parameter file

*)
open ProteomIQon
open ProteomIQon.Domain

let gaborParams : Gabor3DParams =
    {
        sizeX     = 41
        sizeY     = 41
        sigmaX    = 8.0
        sigmaY    = 8.0
        frequency = 2.0
        theta     = 0.06
    }

let quantificationTIMsParams : Dto.QuantificationParams =
    {
        PerformLabeledQuantification = Labeling.Unlabeled
        FragPipe                     = true
        XicExtraction =
            {
                ScanTimeWindow = 2.
                MzWindow_Da    = Window.Estimate
                XicProcessing  = XicProcessing.Gabor3D gaborParams
                TopKPSMs       = None
            }
        BaseLineCorrection = None
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