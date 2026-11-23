---
name: "MzMLToMzLite"
phase: preprocessing
subgroup: File conversion
order: 1
summary: "MzMLToMzLite converts an mzML file into an mzlite file."
inputsFrom: []
---

## Tool Description

MzMLToMzLite converts an mzML file into an mzlite file. It is the first tool of the chain, every later tool reads the mzlite it writes. Mzlite is the SQLite based storage format of [MzIO](https://github.com/CSBiology/MzIO).

## How does MzMLToMzLite work?

1. **MzML to Mzlite**

   The mzlite holds the same spectra and metadata as the mzML, in a form that supports random access to single spectra.

2. **How to get MzML**

   To get mzML from vendor raw files, use [msconvert](https://www.nature.com/articles/nbt.2377). Galaxy Europe runs msconvert in the browser without a local install, see [their guide](https://galaxyproject.eu/posts/2019/03/24/msconvert/).

3. **Converting**

   During the conversion the tool can restrict the spectra to a retention time window and can centroid profile spectra with a wavelet peak picker. Centroiding replaces the many points that describe one peak by a single m/z and intensity pair. The algorithm is described in the [BioFSharp.Mz documentation](https://www.biofsharp.com/BioFSharp.Mz/01_02_signal_detection.html#Centroiding-with-the-wavelet-peak-picker).

## Inputs

| Flag | Description | Source |
| --- | --- | --- |
| `-i` | One or more mzML files, or a directory searched for `*.mzML` and `*.mzml` (not recursively). | MSConvert |
| `-o` | Output directory. Created if it does not exist. | — |
| `-p` | JSON parameter file. | [ProteomIQon default parameters](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/mzMLToMzLiteParams.json) |
| `-c` | Number of files converted in parallel. Default: `1`; no effect for a single file. | — |
| `-f` | switch: rewrite every input file in place and remove `&quot` | — |

## Outputs

| Output | Quantity | Description |
| --- | --- | --- |
| `run.mzlite` | One per input mzML | MzIO SQLite file containing spectra and metadata. Copies MS1 and MS2 spectra only; skips spectra with an empty peak array after peak picking. Scan times are stored in minutes, regardless of the mzML unit. |
| `MzMLToMzLite_log.txt` | One per conversion | Overall conversion log in the output directory. |
| `run_log.txt` | One per input file | Per-run log in the output directory. |

The generated mzlite is read by [PeptideSpectrumMatching](./peptide-spectrum-matching.md), [PSMBasedQuantification](./psm-based-quantification.md) and [AlignmentBasedQuantification](./alignment-based-quantification.md).

## Parameter file

| Parameter | Default | Meaning |
| --- | --- | --- |
| `Compress` | Compression.NoCompression | Compression of the peak arrays in the mzlite. Compression.ZLib and Compression.NumPress are the other options, Compression.NumPressZLib combines both. |
| `StartRetentionTime` | None | With Some t, only spectra with a scan time above t are copied. The comparison uses the time unit of the mzML file. |
| `EndRetentionTime` | None | With Some t, only spectra with a scan time below t are copied, same unit as above. |
| `MS1PeakPicking` | PeakPicking.ProfilePeaks | How MS1 peaks are extracted. PeakPicking.ProfilePeaks copies the peak arrays as they are. PeakPicking.Centroid (CentroidizationMode.Wavelet w) runs the wavelet peak picker with the parameters w. |
| `MS2PeakPicking` | PeakPicking.ProfilePeaks | The same choice for MS2 spectra. |

- CentroidizationMode.Manufacturer is for vendor raw files. For mzML input use ProfilePeaks or the wavelet centroiding.

- When a spectrum type is set to PeakPicking.Centroid (CentroidizationMode.Wavelet w), w is a WaveletPeakPickingParams record. The values below are the MS1 settings from the ion mobility default file and serve as a starting point.

  | Parameter | Example | Meaning |
  | --- | --- | --- |
  | `NumberOfScales` | `3` | Number of wavelet widths the transform is evaluated at. |
  | `YThreshold` | `YThreshold.Fixed 1.0` | Intensity floor for the transform. `YThreshold.MinSpectrumIntensity` uses the lowest intensity in each spectrum instead of a fixed value. |
  | `Centroid_MzTolerance` | `0.1` | m/z tolerance within which neighbouring maxima are merged into one centroid. |
  | `SNRS_Percentile` | `95.0` | Percentile of the signal used to estimate the noise level. |
  | `MinSNR` | `1.0` | Minimum signal-to-noise ratio a centroid must reach. |
  | `RefineMZ` | `false` | Whether the centroid m/z is refined from the profile points around the maximum. |
  | `SumIntensities` | `false` | Whether centroid intensity is the sum of the profile intensities under the peak instead of the maximum. |
  | `PaddingParams` | `None` | Optional padding settings for gaps in sparse profile spectra. |

- `Some p` pads gaps in sparse profile spectra before the transform. The `PaddingParams` record contains `MaximumPaddingPoints` (`int option`), `Padding_MzTolerance` (`float`), `WindowSize` (`int`) and `SpacingPerc` (`float`).
- The MS2 settings in the ion mobility default file use `Some { MaximumPaddingPoints = Some 7; Padding_MzTolerance = 0.05; WindowSize = 150; SpacingPerc = 95.0 }` as an example.

- The default file is [mzMLToMzLiteParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/mzMLToMzLiteParams.json). A file with wavelet centroiding on both levels is [TIMsMzMLtoMzLiteParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/TIMsMzMLtoMzLiteParams.json).

## Parameter generation

From the ProteomIQon repository root, build the library with `dotnet build src/ProteomIQon/ProteomIQon.fsproj -c Release`. Save the following as `generate-params.fsx` in that same directory and run `dotnet fsi generate-params.fsx`. The script writes the JSON file to the system temporary directory.

```fsharp
#r "nuget: Newtonsoft.Json, 13.0.4"
#r "./src/ProteomIQon/bin/Release/net10.0/ProteomIQon.dll"

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
```

## Installing and executing the tool

1. **Installing**

   Install the tool globally with the .NET CLI.

   ```sh
   dotnet tool install --global ProteomIQon.MzMLToMzLite
   ```

2. **Execute**

   - Convert one run

     ```sh
     proteomiqon-mzmltomzlite -i path/to/run.mzML -o path/to/output -p path/to/mzMLToMzLiteParams.json
     ```

   - Convert multiple runs

     ```sh
     proteomiqon-mzmltomzlite -i path/to/run1.mzML path/to/run2.mzML path/to/run3.mzML -o path/to/output -p path/to/mzMLToMzLiteParams.json -c 3
     ```

3. **Help**

   If the XML reader fails because the input contains a stray `&quot`, use `-f` to remove that string from the input file in place.

   ```sh
   proteomiqon-mzmltomzlite -i path/to/run.mzML -o path/to/output -p path/to/mzMLToMzLiteParams.json -f
   ```

4. **All Flags**

   ```sh
   proteomiqon-mzmltomzlite --help
   ```
