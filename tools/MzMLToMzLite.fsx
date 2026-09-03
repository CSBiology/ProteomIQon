(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/MzMLToMzLite.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLite.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLite.ipynb)

# MzMLToMzLite

MzMLToMzLite converts an mzML file into an mzlite file. It is the first tool of the chain, every later tool reads the mzlite it writes.
mzlite is the SQLite based storage format of [MzIO](https://github.com/CSBiology/MzIO). It holds the same spectra and metadata as the mzML, in a form that supports random access to single spectra.
To get mzML from vendor raw files, use [msconvert](https://www.nature.com/articles/nbt.2377).
Galaxy Europe runs msconvert in the browser without a local install, see [their guide](https://galaxyproject.eu/posts/2019/03/24/msconvert/).
During the conversion the tool can restrict the spectra to a retention time window and can centroid profile spectra with a wavelet peak picker. Centroiding replaces the many points that describe one peak by a single m/z and intensity pair. The algorithm is described in the [BioFSharp.Mz documentation](https://www.biofsharp.com/BioFSharp.Mz/01_02_signal_detection.html#Centroiding-with-the-wavelet-peak-picker).

## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | one or more mzML files, or a directory that is searched for `*.mzML` and `*.mzml` (not recursive) | msconvert
`-o` | the output directory, created when missing | &#32;
`-p` | the parameter file in JSON | this page
`-c` | number of files converted at the same time, default 1, no effect on a single file | &#32;
`-f` | switch: rewrite every input file in place and remove `&quot` | &#32;


`-f` rewrites every input file in place and removes the string `&quot` from it. Some mzML exports contain that string in places where the XML reader cannot handle it. Leave the flag off unless the conversion fails on such a file.

The tool writes one `<run>.mzlite` per input `<run>.mzML` into the output directory. Only spectra with MS level 1 or 2 are copied, and spectra whose peak array is empty after peak picking are skipped. Scan times in the mzlite are stored in minutes, whatever unit the mzML used.
The mzlite is read by [PeptideSpectrumMatching](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.html), [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html) and [AlignmentBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html).
The output directory also receives `MzMLToMzLite_log.txt` and one `<run>_log.txt` per input.

## Parameters

Parameter | Default | Meaning
--- | --- | ---
Compress | Compression.NoCompression | Compression of the peak arrays in the mzlite. Compression.ZLib and Compression.NumPress are the other options, Compression.NumPressZLib combines both.
StartRetentionTime | None | With Some t, only spectra with a scan time above t are copied. The comparison uses the time unit of the mzML file.
EndRetentionTime | None | With Some t, only spectra with a scan time below t are copied, same unit as above.
MS1PeakPicking | PeakPicking.ProfilePeaks | How MS1 peaks are extracted. PeakPicking.ProfilePeaks copies the peak arrays as they are. PeakPicking.Centroid (CentroidizationMode.Wavelet w) runs the wavelet peak picker with the parameters w.
MS2PeakPicking | PeakPicking.ProfilePeaks | The same choice for MS2 spectra.


CentroidizationMode.Manufacturer is for vendor raw files. For mzML input use ProfilePeaks or the wavelet centroiding.

When a spectrum type is set to PeakPicking.Centroid (CentroidizationMode.Wavelet w), w is a WaveletPeakPickingParams record. The values below are the MS1 settings from the ion mobility default file and serve as a starting point.

Parameter | Example | Meaning
--- | --- | ---
NumberOfScales | 3 | Number of wavelet widths the transform is evaluated at.
YThreshold | YThreshold.Fixed 1.0 | Intensity floor for the transform. YThreshold.MinSpectrumIntensity uses the lowest intensity of each spectrum instead of a fixed value.
Centroid_MzTolerance | 0.1 | m/z tolerance within which neighbouring maxima are merged into one centroid.
SNRS_Percentile | 95.0 | Percentile of the signal used to estimate the noise level.
MinSNR | 1.0 | Minimum signal to noise ratio a centroid has to reach.
RefineMZ | false | Whether the centroid m/z is refined from the profile points around the maximum.
SumIntensities | false | Whether the centroid intensity is the sum of the profile intensities under the peak instead of the maximum.
PaddingParams | None | Some p pads gaps in sparse profile spectra before the transform. p is a PaddingParams record with MaximumPaddingPoints (int option), Padding_MzTolerance (float), WindowSize (int) and SpacingPerc (float). The MS2 settings of the ion mobility default file use Some { MaximumPaddingPoints = Some 7; Padding_MzTolerance = 0.05; WindowSize = 150; SpacingPerc = 95.0 }.


The default file is [mzMLToMzLiteParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/mzMLToMzLiteParams.json). A file with wavelet centroiding on both levels is [TIMsMzMLtoMzLiteParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/TIMsMzMLtoMzLiteParams.json).

## Writing a parameter file

*)
open ProteomIQon
open ProteomIQon.Domain

let mzMLToMzLiteParams : Dto.MzMLtoMzLiteParams =
    {
        Compress           = Compression.NoCompression
        StartRetentionTime = None
        EndRetentionTime   = None
        MS1PeakPicking     = PeakPicking.ProfilePeaks
        MS2PeakPicking     = PeakPicking.ProfilePeaks
    }

// Replace the temp folder with your project folder.
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "mzMLToMzLiteParams.json")

Json.serializeAndWrite outputPath mzMLToMzLiteParams
(**
## Running the tool

Install the tool with `dotnet tool install --global ProteomIQon.MzMLToMzLite`, then convert one run:

```text
proteomiqon-mzmltomzlite -i path/to/run.mzML -o path/to/output -p path/to/mzMLToMzLiteParams.json
```

Several runs, three of them converted at the same time:

```text
proteomiqon-mzmltomzlite -i path/to/run1.mzML path/to/run2.mzML path/to/run3.mzML -o path/to/output -p path/to/mzMLToMzLiteParams.json -c 3
```

A file that fails to parse because of a stray `&quot`:

```text
proteomiqon-mzmltomzlite -i path/to/run.mzML -o path/to/output -p path/to/mzMLToMzLiteParams.json -f
```

All flags:

```text
proteomiqon-mzmltomzlite --help
```

*)