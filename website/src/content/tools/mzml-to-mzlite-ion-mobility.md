---
name: "MzMLToMzLiteIonMobility"
phase: preprocessing
subgroup: File conversion
order: 2
summary: "MzMLToMzLiteIonMobility converts mzML files from ion mobility instruments such as the Bruker timsTOF into mzlite files."
inputsFrom: []
---

## Tool Description

MzMLToMzLiteIonMobility converts mzML files from ion mobility instruments such as the Bruker timsTOF into mzlite files.

## How does MzMLToMzLiteIonMobility work?

1. **MzML to Mzlite**

   In these files every peak carries an ion mobility value next to its m/z and intensity, and the tool keeps that value in the mzlite. It reads the same parameter record as [MzMLToMzLite](./mzml-to-mzlite.md). That page also explains how to obtain mzML from vendor files.

2. **Preserve peak arrays**

   The spectra are written with their m/z, intensity and ion mobility arrays as they are in the mzML, and Compress decides how those arrays are stored.

## Inputs

These inputs are needed to execute the tool:

| Flag | Meaning | Comes from |
| --- | --- | --- |
| `-i` | one or more mzML files with ion mobility data per peak, or a directory that is searched for `*.mzML` and `*.mzml` (not recursive) | MSConvert |
| `-o` | the output directory, created when missing | - |
| `-p` | the parameter file in JSON | [Github ProteomIQon Core](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/TIMsMzMLtoMzLiteParams.json) |
| `-c` | number of files converted at the same time, default 1, no effect on a single file | - |
| `-f` | switch: rewrite every input file in place and remove `&quot` | - |

- The reader fails on mzML without ion mobility data.
- The tool writes one `run.mzlite` per input into the output directory.
- Spectra with an MS level other than 1 or 2 stop the conversion with an error. Scan times are stored in minutes.
- The mzlite is read by [PeptideSpectrumMatchingTIMs](./peptide-spectrum-matching-tims.md), [PSMBasedQuantificationTIMs](./psm-based-quantification-tims.md) and [AlignmentBasedQuantification](./alignment-based-quantification.md). [MsFraggerToPSM](./ms-fragger-to-psm.md) reads it when identification comes from FragPipe. The output directory also receives `MzMLToMzLite_log.txt` and one `run_log.txt` per input.

## Parameter file

| Parameter | Default | Meaning |
| --- | --- | --- |
| `Compress` | Compression.NoCompression | Compression of the peak arrays in the mzlite. Compression.ZLib and Compression.NumPress are the other options, Compression.NumPressZLib combines both. |
| `StartRetentionTime` | None | Unused by this converter. Leave at `None`. |
| `EndRetentionTime` | None | Unused by this converter. Leave at `None`. |
| `MS1PeakPicking` | PeakPicking.Centroid (CentroidizationMode.Wavelet ms1Wavelet) | Same field as in [MzMLToMzLite](./mzml-to-mzlite.md). The default file contains a wavelet record for MS1. The TIMs converter copies the MS1 peak arrays as they are. |
| `MS2PeakPicking` | PeakPicking.Centroid (CentroidizationMode.Wavelet ms2Wavelet) | Same field as in [MzMLToMzLite](./mzml-to-mzlite.md), with padding settings in the default file. The MS2 peaks are copied as they are. |

- The wavelet parameters are explained on the [MzMLToMzLite page](./mzml-to-mzlite.md). The default file is [TIMsMzMLtoMzLiteParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/TIMsMzMLtoMzLiteParams.json).

## Parameter generation

From the ProteomIQon repository root, build the library with `dotnet build src/ProteomIQon/ProteomIQon.fsproj -c Release`. Save the following as `generate-params.fsx` in that same directory and run `dotnet fsi generate-params.fsx`. The script writes the JSON file to the system temporary directory.

```fsharp
#r "nuget: Newtonsoft.Json, 13.0.4"
#r "./src/ProteomIQon/bin/Release/net10.0/ProteomIQon.dll"

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
```

## Installing and executing the tool

1. **Installing**

   Install the tool globally with the .NET CLI.

   ```sh
   dotnet tool install --global ProteomIQon.MzMLToMzLiteIonMobility
   ```

2. **Execute**

   - Convert one run

     ```sh
     proteomiqon-mzmltomzliteionmobility -i path/to/run.mzML -o path/to/output -p path/to/TIMsMzMLtoMzLiteParams.json
     ```

   - Convert multiple runs

     ```sh
     proteomiqon-mzmltomzliteionmobility -i path/to/run1.mzML path/to/run2.mzML path/to/run3.mzML -o path/to/output -p path/to/TIMsMzMLtoMzLiteParams.json -c 3
     ```

3. **Help**

   If the XML reader fails because the input contains a stray `&quot`, use `-f` to remove that string from the input file in place.

   ```sh
   proteomiqon-mzmltomzliteionmobility -i path/to/run.mzML -o path/to/output -p path/to/TIMsMzMLtoMzLiteParams.json -f
   ```

4. **All Flags**

   ```sh
   proteomiqon-mzmltomzliteionmobility --help
   ```
