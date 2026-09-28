(**
---
title: PeptideSpectrumMatchingTIMs
category: Tools
categoryindex: 1
index: 5
---
*)

(*** hide ***)

(*** condition: prepare ***)
#r "nuget: BioFSharp.Mz, 0.2.1"
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

# PeptideSpectrumMatchingTIMs

PeptideSpectrumMatchingTIMs identifies peptides in timsTOF runs with precursor ion mobility. It accepts mzlite input and follows the search pattern used by MSFragger. [PeptideSpectrumMatching]({{root}}tools/PeptideSpectrumMatching.html) performs the corresponding search for mzlite and mzML files.

The fragment index holds the b and y ions of target peptides only. Reversed decoy peptides use the same index through a water shift. A decoy b ion is a target y ion minus water, and a decoy y ion is a target b ion plus water. The index uses 4 bytes per fragment entry.

For each MS2 spectrum, the tool removes the peaks around the precursor m/z and the peaks below `MinimumPeakRatio` of the base peak. It removes isotope partner peaks and keeps the `TopNPeaks` most intense peaks. Every remaining peak is looked up in the index, restricted to the peptides whose precursor mass fits one of the configured isotope errors. The precursor charge and the selected ion m/z come from the mzlite file. When a spectrum carries no charge, the charge states in `FallbackChargeStates` are searched. The tool reports the `ReportedHitsPerLabel` best targets and the same number of best decoys.

The index matches whole bins of `FragmentIndexBinWidth`, so a peak reaches every fragment within its tolerance plus up to one bin width on either side. Peaks whose tolerance windows lie closer than one bin are merged with the larger intensity before the lookup, so every fragment of a candidate is counted at most once in the index pass. The hyperscore is the X!Tandem hyperscore of BioFSharp.Mz, `ln(sum of matched intensities) + ln(Nb!) + ln(Ny!)`, with intensities relative to a base peak of 100 and a floor of 0. For the reported hits the tool computes it again exactly with the fragment tolerance, which gives the `Hyperscore` and `MatchedIons` columns and decides the rank. `MinMatchedFragments` applies to that exact count.

The expectation value follows the model of X!Tandem, which MSFragger also follows. The index pass hyperscores of all candidates of a spectrum with at least `MinFragmentsModelling` matched fragments form a histogram with unit bins on the X!Tandem scale, which is 4 log10 of the product form, so 4 / ln 10 times the natural log form above. The `Expectscore` of a hit is evaluated at its index pass hyperscore, the same quantity the model was built from, so it can order the hits of a spectrum differently from `Hyperscore`. The model removes candidates that sit above a gap in the upper part of the survival function. It fits a line to the log10 survival function from half of the remaining count down to ten survivors. When fewer than 200 candidates remain, the model uses X!Tandem's default line, `10^(3.5 - 0.18 x)`. The expectation value has a floor of `1e-15`.

The reported targets and decoys then receive the SEQUEST-like score used by PeptideSpectrumMatching. The tool also computes its Andromeda-like and X!Tandem-like scores for these hits. The `NormDeltaBestToRest` and `NormDeltaNext` columns are computed over the target and the decoy form of every peptide the index pass selected, as a target or as a decoy, which includes candidates that the exact rescoring drops afterwards. PeptideSpectrumMatching computes these columns over every peptide in the precursor window, so values from the two tools are not comparable.

Parallelization follows the library convention at the file level. `-c` sets the number of runs processed at the same time and defaults to 1, as in PeptideSpectrumMatching. The index build at the start of an invocation uses that many workers as well, after that the tool searches one run on one thread. The peptide database has to use monoisotopic masses. Most of the time per spectrum goes into the classic scoring of the reported hits, so the run time scales with `ReportedHitsPerLabel`.

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | one or more `.mzlite` files, or a directory that is searched for `*.mzlite` | [MzMLToMzLiteIonMobility]({{root}}tools/MzMLToMzLiteIonMobility.html) |
| `-d` | the SQLite peptide database | [PeptideDB]({{root}}tools/PeptideDB.html) |
| `-o` | the output directory, created when missing | |
| `-p` | the parameter file in JSON | this page |
| `-c` | number of runs processed at the same time, default 1 | |

All flags except `-c` are mandatory. The input directory search is not recursive.

The tool accepts `.mzlite` files and directories. An mzML path stops the tool with a message.

The tool writes one `<run>.psm` per input into the output directory. The file is a tab separated table with a header. The tool overwrites it when it already exists. Label is 1 for a target and -1 for a reversed decoy. `MissCleavages` is -1, as in PeptideSpectrumMatching. PSMStatistics recomputes it.

The output record is `Dto.PeptideSpectrumMatchingResult`, the same record written by [PeptideSpectrumMatching]({{root}}tools/PeptideSpectrumMatching.html). It has 28 columns. The record has the 23 classic PeptideSpectrumMatching columns, followed by these five columns.

| Column | Meaning |
|--------|---------|
| `IonMobility` | Inverse reduced ion mobility, 1/K0, of the precursor. It is `NaN` when the run has no ion mobility value. |
| `Hyperscore` | Hyperscore from the matched b and y ions. |
| `Expectscore` | Expectation value from the candidate hyperscore model. |
| `MatchedIons` | Number of matched b and y ions used for the hyperscore. |
| `TotalIons` | Number of theoretical b and y ions of the peptide. |

`ScanNr` is the position of the MS2 spectrum in scan time order. `PSMId` is the spectrum id with spaces replaced by "-", followed by `_ScanNr_charge_rank`. The rank is the position among the reported targets or decoys of that spectrum, with the best hyperscore first. `AbsDeltaMass` is the mass error after the isotope error that fits the candidate best.

PSMStatistics reads it directly and uses `Hyperscore` and `Expectscore` when they are present. Its `.qpsm` output does not carry the ion mobility yet, and [PSMBasedQuantificationTIMs]({{root}}tools/PSMBasedQuantificationTIMs.html) still expects the FragPipe column layout of [MsFraggerToPSM]({{root}}tools/MsFraggerToPSM.html), so the quantification of a run identified with this tool needs those two tools to be adapted first.

The output directory receives `PeptideSpectrumMatchingTIMs_log.txt` and one `<run>_log.txt` per input.

## Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `PrecursorTolerancePPM` | `20.0` | Precursor neutral mass tolerance in ppm. |
| `FragmentTolerancePPM` | `20.0` | Fragment m/z tolerance in ppm. |
| `IsotopeErrors` | `[0, 1, 2]` | Precursor isotope errors searched. `0` is the monoisotopic peak. |
| `MaxFragmentCharge` | `2` | Highest fragment charge, at most the precursor charge minus one and at least 1. |
| `FallbackChargeStates` | `[2, 3]` | Charge states searched when the spectrum has no precursor charge. |
| `TopNPeaks` | `150` | Number of most intense peaks kept per spectrum after filtering. |
| `MinimumPeakRatio` | `0.01` | Peaks below this fraction of the base peak are removed. |
| `RemovePrecursorRange` | `1.5` | Peaks within this m/z distance of the precursor m/z are removed. |
| `Deisotope` | `true` | Removes isotope partner peaks before matching. |
| `MinimumPeaks` | `15` | Spectra with fewer peaks after preprocessing are skipped. |
| `MinMatchedFragments` | `4` | Minimum matched fragments for a candidate to be reported. |
| `MinFragmentsModelling` | `1` | Minimum matched fragments for a candidate to enter the expectation value model. |
| `ReportedHitsPerLabel` | `10` | Number of best targets and best decoys reported per spectrum. The run time scales with this number. |
| `FragmentIndexBinWidth` | `0.01` | Fragment index bin width in Da. |
| `nTerminalSeries` | `B` | N-terminal ion series used by the classic scoring functions. |
| `cTerminalSeries` | `Y` | C-terminal ion series used by the classic scoring functions. |
| `Andromeda` | `{ PMinPMax = 4, 10; MatchingIonTolerancePPM = 100.0 }` | Andromeda-like scoring settings. |

`nTerminalSeries` and `cTerminalSeries` only affect the classic scores of the reported hits. The fragment index always uses b and y ions.

Andromeda is a `ProteomIQon.Domain.AndromedaParams` record.

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `PMinPMax` | `4, 10` | Lowest and highest number of most intense peaks kept per 100 Da window. Every count in this range is tried and the best score is kept. |
| `MatchingIonTolerancePPM` | `100.0` | Tolerance in ppm for matching a theoretical fragment to a measured peak. |

The default file is [peptideSpectrumMatchingTIMsParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/peptideSpectrumMatchingTIMsParams.json).

## Writing a parameter file
*)

open BioFSharp.Mz
open ProteomIQon

let andromedaParams : Domain.AndromedaParams =
    {
        PMinPMax                = 4, 10
        MatchingIonTolerancePPM = 100.0
    }

let peptideSpectrumMatchingTIMsParams : Dto.PeptideSpectrumMatchingTIMsParams =
    {
        PrecursorTolerancePPM  = 20.0
        FragmentTolerancePPM   = 20.0
        IsotopeErrors          = [0; 1; 2]
        MaxFragmentCharge      = 2
        FallbackChargeStates   = [2; 3]
        TopNPeaks               = 150
        MinimumPeakRatio        = 0.01
        RemovePrecursorRange    = 1.5
        Deisotope               = true
        MinimumPeaks            = 15
        MinMatchedFragments     = 4
        MinFragmentsModelling   = 1
        ReportedHitsPerLabel    = 10
        FragmentIndexBinWidth   = 0.01
        nTerminalSeries         = NTerminalSeries.B
        cTerminalSeries         = CTerminalSeries.Y
        Andromeda               = andromedaParams
    }

// Replace the temp folder with your project folder.
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "peptideSpectrumMatchingTIMsParams.json")

Json.serializeAndWrite outputPath peptideSpectrumMatchingTIMsParams

(**
## Running the tool

Install the tool with `dotnet tool install --global ProteomIQon.PeptideSpectrumMatchingTIMs`, then score one run:

```text
proteomiqon-peptidespectrummatchingtims -i path/to/run.mzlite -d path/to/AraTest.db -o path/to/output -p path/to/peptideSpectrumMatchingTIMsParams.json
```

Process three runs with `-c 3`, so three runs are processed at the same time:

```text
proteomiqon-peptidespectrummatchingtims -i path/to/run1.mzlite path/to/run2.mzlite path/to/run3.mzlite -d path/to/AraTest.db -o path/to/output -p path/to/peptideSpectrumMatchingTIMsParams.json -c 3
```

All flags:

```text
proteomiqon-peptidespectrummatchingtims --help
```
*)
