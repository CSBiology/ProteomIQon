(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/MzMLToMzLiteIonMobility.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLiteIonMobility.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLiteIonMobility.ipynb)

# MzMLToMzLiteIonMobility

MzMLToMzLiteIonMobility converts mzML files from ion mobility instruments such as the Bruker timsTOF into mzlite files. In these files every peak carries an ion mobility value next to its m/z and intensity, and the tool keeps that value in the mzlite.
It is the ion mobility counterpart of [MzMLToMzLite](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLite.html), which also explains how to get mzML from raw files and what mzlite is.
The tool reads the same parameter record as MzMLToMzLite. The spectra are written with their m/z, intensity and ion mobility arrays as they are in the mzML, and Compress decides how those arrays are stored.

## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | one or more mzML files with ion mobility data per peak, or a directory that is searched for `*.mzML` and `*.mzml` (not recursive) | msconvert
`-o` | the output directory, created when missing | &#32;
`-p` | the parameter file in JSON, by default TIMsMzMLtoMzLiteParams.json | this page
`-c` | number of files converted at the same time, default 1 | &#32;
`-f` | switch: rewrite every input file in place and remove `&quot` | &#32;


The reader fails on mzML without ion mobility data. `-f` rewrites every input file in place and removes the string `&quot`, for mzML exports where that string breaks the XML reader.

The tool writes one `<run>.mzlite` per input into the output directory. Spectra with an MS level other than 1 or 2 stop the conversion with an error. Scan times are stored in minutes.
The mzlite is read by [MsFraggerToPSM](https://csbiology.github.io/ProteomIQon/tools/MsFraggerToPSM.html) and [PSMBasedQuantificationTIMs](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantificationTIMs.html).
The output directory also receives `MzMLToMzLite_log.txt` and one `<run>_log.txt` per input.

## Parameters

Parameter | Default | Meaning
--- | --- | ---
Compress | Compression.NoCompression | Compression of the peak arrays in the mzlite. Compression.ZLib and Compression.NumPress are the other options, Compression.NumPressZLib combines both.
StartRetentionTime | None | Same field as in MzMLToMzLite. Keep None, the spectra are copied as they are.
EndRetentionTime | None | Same field as in MzMLToMzLite. Keep None.
MS1PeakPicking | PeakPicking.Centroid (CentroidizationMode.Wavelet ms1Wavelet) | Same field as in MzMLToMzLite. The default file carries the wavelet record shown in the script below, and the MS1 peaks are copied as they are.
MS2PeakPicking | PeakPicking.Centroid (CentroidizationMode.Wavelet ms2Wavelet) | Same field as in MzMLToMzLite, with padding settings in the default file. The MS2 peaks are copied as they are.


The wavelet parameters are explained on the [MzMLToMzLite](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLite.html) page. The default file is [TIMsMzMLtoMzLiteParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/TIMsMzMLtoMzLiteParams.json).

## Writing a parameter file

*)
open ProteomIQon
open ProteomIQon.Domain

let ms1Wavelet : WaveletPeakPickingParams =
    {
        NumberOfScales       = 3
        YThreshold           = YThreshold.Fixed 1.0
        Centroid_MzTolerance = 0.1
        SNRS_Percentile      = 95.0
        MinSNR               = 1.0
        RefineMZ             = false
        SumIntensities       = false
        PaddingParams        = None
    }

let ms2Padding : PaddingParams =
    {
        MaximumPaddingPoints = Some 7
        Padding_MzTolerance  = 0.05
        WindowSize           = 150
        SpacingPerc          = 95.0
    }

let ms2Wavelet : WaveletPeakPickingParams =
    {
        NumberOfScales       = 10
        YThreshold           = YThreshold.MinSpectrumIntensity
        Centroid_MzTolerance = 0.1
        SNRS_Percentile      = 95.0
        MinSNR               = 1.0
        RefineMZ             = false
        SumIntensities       = false
        PaddingParams        = Some ms2Padding
    }

let timsMzMLToMzLiteParams : Dto.MzMLtoMzLiteParams =
    {
        Compress           = Compression.NoCompression
        StartRetentionTime = None
        EndRetentionTime   = None
        MS1PeakPicking     = PeakPicking.Centroid (CentroidizationMode.Wavelet ms1Wavelet)
        MS2PeakPicking     = PeakPicking.Centroid (CentroidizationMode.Wavelet ms2Wavelet)
    }

// Replace the temp folder with your project folder.
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "TIMsMzMLtoMzLiteParams.json")

Json.serializeAndWrite outputPath timsMzMLToMzLiteParams
(**
## Running the tool

Install the tool with `dotnet tool install --global ProteomIQon.MzMLToMzLiteIonMobility`, then convert one run:

```text
proteomiqon-mzmltomzliteionmobility -i path/to/run.mzML -o path/to/output -p path/to/TIMsMzMLtoMzLiteParams.json
```

Several runs, three of them converted at the same time:

```text
proteomiqon-mzmltomzliteionmobility -i path/to/run1.mzML path/to/run2.mzML path/to/run3.mzML -o path/to/output -p path/to/TIMsMzMLtoMzLiteParams.json -c 3
```

A file that fails to parse because of a stray `&quot`:

```text
proteomiqon-mzmltomzliteionmobility -i path/to/run.mzML -o path/to/output -p path/to/TIMsMzMLtoMzLiteParams.json -f
```

All flags:

```text
proteomiqon-mzmltomzliteionmobility --help
```

*)