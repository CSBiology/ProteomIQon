---
name: "PSMBasedQuantificationTIMs"
phase: peptide-quantification
order: 2
summary: "PSMBasedQuantificationTIMs is the PSMBasedQuantification variant for timsTOF runs."
inputsFrom: [mzml-to-mzlite, peptide-db]
---

## Tool Description

PSMBasedQuantificationTIMs is the PSMBasedQuantification variant for timsTOF runs. On these instruments a peptide ion is a peak in two dimensions, retention time and ion mobility, so the tool extracts the signal around every identification in both dimensions.

## How does PSMBasedQuantificationTIMs work?

1. **XIC Point Calculation and Ion Mobility Window**

   An XIC point is the sum of the intensity of all peaks inside the m/z window (`MzWindow_Da`) and inside the ion mobility window (`IonMobilityWindow`), so both widths are part of the quantity.
   With `IonMobilityWindow` set to `Window.Fixed v`, the ion mobility window is ±v around the ion mobility centre of the ion. The centre is the intensity weighted ion mobility of the signal within v of the weighted ion mobility of the identifications, around their weighted scan time (see `CentrePSMs`). Without such signal it is the weighted ion mobility of the identifications.

2. **Peak Detection and Quantification**
   The points form a retention time [ion chromatogram](https://www.biofsharp.com/BioFSharp.Mz/03_01_quantification.html#Building-an-XIC-trace).
   As PSMBasedQuantification does, the tool detects the peaks on the trace and [fits and integrates](https://www.biofsharp.com/BioFSharp.Mz/03_01_quantification.html#Fitting-and-integrating-the-main-peak) the peak that contains the weighted scan time of the identifications.
   The measured isotope pattern of a row (the `KLDiv` and `IsotopicPattern` columns) is read with the same windows, on the spectrum closest to the fitted scan time.

3. **Output**

   The output has the same layout as the output of PSMBasedQuantification, so the downstream tools do not care which of the two produced it.

## Inputs

These inputs are needed to execute the tool:

| Flag | Meaning | Comes from |
| --- | --- | --- |
| `-i` | one or more `.mzlite` files with ion mobility, or a directory that is searched for `*.mzlite` | [MzMLToMzLiteIonMobility](./mzml-to-mzlite-ion-mobility.md) |
| `-ii` | one or more `.qpsm` files in the PSMStatistics layout, or a directory searched for `*.qpsm` | [PSMStatistics](./psm-statistics.md) or [MsFraggerToPSM](./ms-fragger-to-psm.md) |
| `-o` | the output directory, created when missing | - |
| `-d` | the SQLite peptide database | [PeptideDB](./peptide-db.md) |
| `-p` | the parameter file in JSON | [Github ProteomIQon Core](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/QuantificationTIMsParams.json) |
| `-mf` | switch: pair files by base name instead of by position | - |
| `-dc` | switch: write diagnostic charts | - |
| `-z` | switch: pack the chart directory into a zip file and delete it | - |
| `-c` | number of runs quantified at the same time, default 1 | - |

- The tool accepts instrument files ending in `.mzlite`. It checks the `.qpsm` header against the PSMStatistics layout, which MsFraggerToPSM also writes. The tool uses the `PSMId` to retrieve ion mobility from the corresponding MS2 spectrum in the `.mzlite`. It warns about missing spectra or ion mobility values and skips ions without a usable identification. If none of the identifications in a file has ion mobility, the tool stops.
- With one file per flag the pair is used as given. With lists, the n-th mzlite file is paired with the n-th qpsm file unless you pass -mf, which pairs files by their base name and skips mzlite files without a partner.
- For every run the tool writes `run.quant`, tab separated with a header.
- With `-dc` the tool writes charts into `run_plots`, and `-z` packs that directory into `run_plots.zip` and deletes it.
- Delete an existing `run.quant` in the output directory before a rerun, the tool appends to it. Logs go to `PSMBasedQuantification_log.txt` and `run_log.txt` in the output directory.

## Parameter file

The parameter record is shared with PSMBasedQuantification. The following default values are specific to timsTOF runs.

### General settings

| Parameter | Default | Description |
|-----------|---------|---------|
| `PerformLabeledQuantification` | `Labeling.Unlabeled` | Label free quantification of the identified ion. The three labeled cases work as in PSMBasedQuantification and quantify the partner ion too. |
| `XicExtraction.ScanTimeWindow` | `2.0` | Half width of the retention time range around the identification that is extracted, in the scan time unit of the run. |
| `XicExtraction.MzWindow_Da` | `Window.Fixed 0.008` | Half width of the m/z window, in Da. Every peak inside it counts in an XIC point. The m/z correction and `Window.Estimate` work as in PSMBasedQuantification. |
| `XicExtraction.XicProcessing` | `XicProcessing.Wavelet waveletParams` | Peak detection on the retention time trace, summed over the ion mobility window. The wavelet parameters are listed below. `SecondDerivative` is the other supported case. `Gabor3D` is in development and does not work currently. |
| `XicExtraction.TopKPSMs` | `None` | `Some k` keeps only the k identifications with the highest `ModelScore` per peptide ion when placing the peak. |
| `XicExtraction.IonMobilityWindow` | `Window.Fixed 0.025` | Half width v of the ion mobility window around the ion mobility centre of the ion, in 1/K0. `Window.Estimate` is not available here and stops the tool with an error that asks for `Fixed`. |
| `BaseLineCorrection` | `Some { MaxIterations = 10; Lambda = 6; P = 0.05 }` | Asymmetric least squares baseline subtraction before peak detection. `None` skips it. |
| `QualityFilterLog2Limit` | `2.0` | The quality filter keeps rows with -limit < log2 ((area / median area) / (apex / median apex)) < limit. Rows with a non-finite or non-positive area, apex or ratio are dropped. |
| `CentrePSMs` | `3` | Number of identifications per peptide ion with an ion mobility, ranked by the summed intensity of their matched fragment ions, whose intensity weighted scan time and ion mobility locate the trace. It is applied after `TopKPSMs`, which ranks by `ModelScore`. At least 1. |
| `RequirePeakAtScanTime` | `true` | `true` quantifies only a peak that contains the weighted scan time of the identifications, so an ion without such a peak gets no row. `false` takes the peak with the closest apex when no peak contains it. |

The wavelet parameters:

| Parameter | Default | Description |
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

**Default parameter file:** [Default JSON File](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/QuantificationTIMsParams.json)

For a given dataset, use the same fixed m/z and ion mobility windows in [AlignmentBasedQuantification](./alignment-based-quantification.md). The appropriate windows depend on the instrument and acquisition method.

## Outputs

[QuantBasedAlignment](./quant-based-alignment.md), [AlignmentBasedQuantification](./alignment-based-quantification.md), [AlignmentBasedQuantStatistics](./alignment-based-quant-statistics.md) and [JoinQuantPepIonsWithProteins](./join-quant-with-peptide-ions.md) read these `.quant` files.

The output table contains 47 columns. The columns are grouped by category below.

### Identification and PSM confidence

| Column | Description |
| --- | --- |
| `StringSequence` | Amino Acid sequence (one-letter code) from PeptideDB which matches the psm candidate |
| `GlobalMod` | Indicator, if a peptide ion species is labeled or unlabeled |
| `Charge` | Precursor ion charge state |
| `PepSequenceID` | Unique identifier of the unmodified peptide sequence, which points to PeptideDB |
| `ModSequenceID` | Unique identifier of the modified peptide sequence (including PTMs, e.g., methylation), which points to PeptideDB |
| `PrecursorMZ` | Precursor ion mass-to-charge ratio (m/z) |
| `MeasuredMass` | Measured mass of precursor ions in Dalton |
| `TheoMass` | Theoretical peptide mass in the spectrum (based on amino acid composition) in Dalton |
| `AbsDeltaMass` | Absolute mass deviation between theoretical and measured mass (mass error) |
| `MeanPercolatorScore` | Mean `ModelScore` of the identifications of the peptide ion |
| `QValue` | Q-value (False-Discovery-Rate) based on combined model scores |
| `PEPValue` | Posterior Error Probability |
| `ProteinNames` | Protein Names out of the PeptideDB |

### Light quantification (8)

| Column | Description |
| --- | --- |
| `QuantMz_Light` | Mass-to-charge ratio from the extracted XIC of the unlabeled peptide ion |
| `Quant_Light` | Peak area from the extracted XIC of the unlabeled peptide ion |
| `MeasuredApex_Light` | Measured Apex of unlabeled peak |
| `Seo_Light` | Standard Error of Prediction for quantification |
| `Params_Light` | Describes the shape of the fitted peak for the unlabeled version of a peptide. Contains information about peak height, estimated peak time, and peak width. |
| `Difference_SearchRT_FittedRT_Light` | Difference between the retention time originally determined using PSMs and the retention time calculated from the measured peak |
| `KLDiv_Observed_Theoretical_Light` | Describes how well the measured unlabeled isotope pattern matches the theoretically expected pattern. |
| `KLDiv_CorrectedObserved_Theoretical_Light` | Corrected Kullback–Leibler divergence of isotopic pattern of unlabeled peptide ion versus the measured pattern |

### Heavy quantification (8)

| Column | Description |
| --- | --- |
| `QuantMz_Heavy` | Mass-to-charge ratio from the extracted XIC of the labeled peptide ion |
| `Quant_Heavy` | Peak area from the extracted XIC of the labeled peptide ion |
| `MeasuredApex_Heavy` | Measured Apex of labeled peak |
| `Seo_Heavy` | Standard Error of Prediction of quantification |
| `Params_Heavy` | Describes the shape of the fitted peak for the labeled version of a peptide. Contains information about peak height, estimated peak time, and peak width. |
| `Difference_SearchRT_FittedRT_Heavy` | Difference between the retention time originally determined using PSMs and the retention time calculated from the measured peak |
| `KLDiv_Observed_Theoretical_Heavy` | Describes how well the measured labeled isotope pattern matches the theoretically expected pattern. |
| `KLDiv_CorrectedObserved_Theoretical_Heavy` | Corrected Kullback–Leibler divergence of isotopic pattern of labeled peptide ion versus the measured pattern |

### Light isotope patterns and traces (6)

| Column | Description |
| --- | --- |
| `IsotopicPatternMz_Light` | M/z-distribution of isotope cluster in spectra for unlabeled peaks |
| `IsotopicPatternIntensity_Observed_Light` | Observed intensity distribution of unlabeled isotopic clusters |
| `IsotopicPatternIntensity_Corrected_Light` | Corrected intensity distribution of unlabeled isotopic clusters |
| `RtTrace_Light` | Retention Time Profile of unlabeled peptide ion |
| `IntensityTrace_Observed_Light` | Intensity trace over for an observed unlabeled isotopic peak cluster |
| `IntensityTrace_Corrected_Light` | Corrected intensity trace between an unlabeled peak and a measured isotopic peak cluster |

### Heavy isotope patterns and traces (6)

| Column | Description |
| --- | --- |
| `IsotopicPatternMz_Heavy` | M/z-distribution of isotope cluster in spectra for labeled peaks |
| `IsotopicPatternIntensity_Observed_Heavy` | Observed intensity distribution of labeled isotopic clusters |
| `IsotopicPatternIntensity_Corrected_Heavy` | Corrected intensity distribution of labeled isotopic clusters |
| `RtTrace_Heavy` | Retention Time Profile of labeled peptide ion |
| `IntensityTrace_Observed_Heavy` | Intensity trace for an observed labeled isotopic peak cluster |
| `IntensityTrace_Corrected_Heavy` | Corrected intensity trace between a labeled peak and a measured isotopic peak cluster |

### Light / Heavy comparison and quantification source (2)

| Column | Description |
| --- | --- |
| `Correlation_Light_Heavy` | Correlation calculated based on Pearson between unlabeled and labeled peaks |
| `QuantificationSource` | If spectra are found over Alignments or Peptide Spectrum Matching |

### Alignment quality (2)

Alignment scores and associated q-values

| Column | Description |
| --- | --- |
| `AlignmentScore` | `NaN` for PSM-based rows. For transferred rows, the probability from AlignmentBasedQuantStatistics that the quantified peak belongs to the ion. |
| `AlignmentQValue` | `0` for PSM-based rows. For transferred rows, the q-value from AlignmentBasedQuantStatistics. |

### Ion Mobility 
| Column | Description |
| --- | --- |
|IonMobility | Ion mobility of the identifications, weighted by their matched fragment intensity, NaN without ion mobility |
| IonMobilityCentre | Ion mobility the XIC of the row is summed around, NaN without an ion mobility window |

## Parameter generation

From the ProteomIQon repository root, build the library with `dotnet build src/ProteomIQon/ProteomIQon.fsproj -c Release`. Save the following as `generate-params.fsx` in that same directory and run `dotnet fsi generate-params.fsx`. The script writes the JSON file to the system temporary directory.

```fsharp
#r "nuget: FSharp.Stats, 0.6.0"
#r "nuget: Newtonsoft.Json, 13.0.4"
#r "./src/ProteomIQon/bin/Release/net10.0/ProteomIQon.dll"

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
                MzWindow_Da    = Window.Fixed 0.008
                XicProcessing  = XicProcessing.Wavelet waveletParams
                TopKPSMs       = None
                IonMobilityWindow = Window.Fixed 0.025
            }
        BaseLineCorrection = Some { MaxIterations = 10; Lambda = 6; P = 0.05 }
        QualityFilterLog2Limit = 2.
        CentrePSMs = 3
        RequirePeakAtScanTime = true
    }

// Replace the temp folder with your project folder.
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "QuantificationTIMsParams.json")

Json.serializeAndWrite outputPath quantificationTIMsParams
```

## Installing and executing the tool

1. **Installing**

   Install the tool globally with the .NET CLI.

   ```sh
   dotnet tool install --global ProteomIQon.PSMBasedQuantificationTIMs
   ```

2. **Execute**

     ```sh
     proteomiqon-psmbasedquantificationtims -i path/to/run.mzlite -ii path/to/run.qpsm -d path/to/database.db -o path/to/output -p path/to/QuantificationTIMsParams.json
     ```

3. **Running the tool for multiple runs, for three at the same time**

   ```sh
   proteomiqon-psmbasedquantificationtims -i path/to/run1.mzlite path/to/run2.mzlite path/to/run3.mzlite -ii path/to/run1.qpsm path/to/run2.qpsm path/to/run3.qpsm -d path/to/database.db -o path/to/output -p path/to/QuantificationTIMsParams.json -c 3
   ```

4. **Two directories, paired by file name**

   ```sh
   proteomiqon-psmbasedquantificationtims -i path/to/mzliteFolder -ii path/to/qpsmFolder -d path/to/database.db -o path/to/output -p path/to/QuantificationTIMsParams.json -c 3 -mf
   ```

5. **All Flags**

   ```sh
   proteomiqon-psmbasedquantificationtims --help
   ```
