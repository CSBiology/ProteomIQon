---
name: "PeptideSpectrumMatchingTIMs"
phase: peptide-identification
subgroup: Spectrum matching
order: 2
summary: "PeptideSpectrumMatchingTIMs identifies peptides from timsTOF runs using precursor ion mobility data."
inputsFrom: [mzml-to-mzlite-ion-mobility, peptide-db]
---

## Tool Description

PeptideSpectrumMatchingTIMs identifies peptides from timsTOF runs using precursor ion mobility data. It reads `.mzlite` files and applies a fragment-index search approach related to MSFragger. For mzlite or mzML searches without this TIMs-specific workflow, see [PeptideSpectrumMatching ↗](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.html).

The program merges the mobility scans of each MS2 spectrum, selects candidate peptides from a fragment index, ranks their spectrum matches and generates one `.psm` result file for each run.

## How does PeptideSpectrumMatchingTIMs work?

1. **Combine ion mobility scans**

   A timsTOF mzlite MS2 spectrum stores peaks by mobility scan. Since the fragment ions inherit the precursor's mobility, the tool combines these scans before matching. Working from low to high m/z, it groups peaks lying within `--mobility-merge-ppm` of the first peak of the group. Group intensity is summed, while the representative m/z is intensity-weighted.

   The default merging width equals `FragmentTolerancePPM`. Setting `--mobility-merge-ppm 0` combines only exactly equal m/z values; it does *not* disable summation. Spectra without peak-level mobility information pass through unchanged. Separate PASEF measurements remain separate MS2 spectra.

2. **Build the target fragment index**

   The index stores theoretical b and y fragments for target peptides. Reversed decoys can share that index: a decoy b fragment maps to a target y fragment minus water, whereas a decoy y fragment maps to a target b fragment plus water.

   Only peptides in the measured runs' precursor-mass windows are indexed, accounting for isotope errors and precursor tolerance. Fragment masses are capped by the largest applicable precursor window plus 1 Da or the optional `-f` value, and never above the largest peptide mass in the database plus 1 Da. Each fragment entry occupies 4 bytes.

   ### Memory management: index slices

   With `-s`, the peptide table is divided into that many index slices, built and searched sequentially so only one slice is held in memory. Normally the spectra are reread for each slice. Adding `-k` caches preprocessed spectra instead (roughly 2.4 kB per spectrum and charge). Slicing changes memory and I/O requirements, not the results.

3. **Preprocess the spectra and retrieve candidates**

   For each spectrum, remove peaks near the precursor m/z, peaks below the `MinimumPeakRatio` fraction of the base peak, and isotope partners; then retain the `TopNPeaks` strongest peaks. The remaining peaks query the fragment index, while peptide candidates must fit the precursor mass for at least one configured isotope error.

   Precursor m/z and charge normally come from the mzlite input. Missing charge information triggers the configured `FallbackChargeStates`. The search retains the top `ReportedHitsPerLabel` targets and the same number of decoys.

   ### Index binning and duplicate-hit prevention

   Lookup operates on complete `FragmentIndexBinWidth` bins, which can extend the effective candidate-retrieval interval by up to one bin width on either side of the specified tolerance. Nearby peak tolerance ranges closer than a bin are merged using the stronger intensity, preventing repeated fragment counts during index lookup.

4. **Calculate and refine the hyperscore**

   The initial X!Tandem-style hyperscore is calculated using the BioFSharp.Mz implementation. With base-peak-relative intensities scaled to 100 and a minimum intensity floor of zero, the natural-log version is:

   `ln(sum of matched intensities) + ln(Nb!) + ln(Ny!)`

   `Nb` and `Ny` denote the counts of matched b and y fragments. For reported candidates, the tool repeats matching at the exact fragment tolerance. This second calculation determines `Hyperscore`, `MatchedIons` and hyperscore rank. `MinMatchedFragments` is checked against that exact match count.

5. **Model expectation values**

   The expectation-value model follows X!Tandem's approach, also used by MSFragger. For each spectrum, index-pass candidate scores with at least `MinFragmentsModelling` matched fragments populate a histogram with bins of width one on the X!Tandem score scale.

   The histogram score scale is `4 × log10(product-form hyperscore)`, equivalent to `4 / ln(10)` times the natural-log expression above. The fit discards candidates above a gap in the upper survival-function tail and fits a line to the log10 survival function, from half of the remaining candidate count down to ten survivors.

   When fewer than 200 candidates remain, the X!Tandem fallback is `10^(3.5 - 0.18 × x)`. Expectation values cannot go below `1e-15`. Importantly, `Expectscore` is evaluated using the *index-pass* hyperscore, whereas `Hyperscore` is the exact rescored value. Consequently, the two columns need not rank candidates identically.

6. **Compute classic scores and score differences**

   The reported targets and decoys also receive the SEQUEST-like score from PeptideSpectrumMatching and Andromeda-like and X!Tandem-like scores. The `NormDeltaBestToRest` and `NormDeltaNext` metrics account for both target and decoy forms of all peptides selected during index retrieval, including candidates subsequently removed by exact rescoring.

   In ordinary PeptideSpectrumMatching, these differences use all peptides in the precursor window. Therefore, those normalized delta values are **not directly comparable** between the two tools.

7. **Process multiple runs**

   With `-c`, the tool handles groups of runs concurrently, using one thread per run and one shared index per group. Building the index uses one worker per run. A new group starts only after all runs in the current group finish; the default is `-c 1`.

   The peptide database must contain **monoisotopic masses**. Classic rescoring accounts for much of the processing time, so increasing `ReportedHitsPerLabel` also increases runtime.

## Inputs

| Flag | Meaning | Comes from |
| --- | --- | --- |
| `-i` | One or more `.mzlite` files, or a directory containing `*.mzlite` files | [MzMLToMzLiteIonMobility ↗](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLiteIonMobility.html) |
| `-d` | SQLite peptide database | [PeptideDB ↗](https://csbiology.github.io/ProteomIQon/tools/PeptideDB.html) |
| `-o` | Output directory, automatically created if absent | — |
| `-p` | JSON file containing search parameters | Parameters section below |
| `-c` | Concurrent runs (one thread per run) | Default: `1` |
| `-s` | Number of sequential fragment-index slices | Default: `1` |
| `-k` | Cache preprocessed spectra between index slices | Optional switch |
| `-f` | Maximum indexed fragment mass, in Da | Default: largest precursor window + 1 Da |
| `--mobility-merge-ppm` | m/z grouping width when combining mobility scans, in ppm | Default: `FragmentTolerancePPM` |

- `-i`, `-d`, `-o` and `-p` are mandatory. The input directory search is not recursive.

**Input restrictions:** Directory scanning is not recursive. An mzML input path produces an error; only mzlite inputs are supported.

## Parameters

The search parameters are provided through the JSON file passed with `-p`. The table lists every parameter and its documented default.

| Parameter | Default | Description |
| --- | --- | --- |
| `PrecursorTolerancePPM` | 20.0 | Tolerance for neutral precursor mass, in ppm |
| `FragmentTolerancePPM` | 20.0 | Tolerance for fragment m/z matching, in ppm |
| `IsotopeErrors` | `[0, 1, 2]` | Allowed precursor isotope offsets; zero denotes the monoisotopic peak |
| `MaxFragmentCharge` | 2 | Upper fragment charge, restricted to the precursor charge minus one and at least one |
| `FallbackChargeStates` | `[2, 3]` | Charge states used when precursor charge is missing; must contain at least one |
| `TopNPeaks` | 150 | Maximum strongest peaks retained after preprocessing |
| `MinimumPeakRatio` | 0.01 | Discard peaks below this base-peak intensity fraction |
| `RemovePrecursorRange` | 1.5 | Excluded m/z distance around the precursor |
| `Deisotope` | `true` | Remove isotope-related peaks before searching |
| `MinimumPeaks` | 15 | Skip spectra with fewer peaks after filtering |
| `MinMatchedFragments` | 4 | Minimum exact fragment matches required to report a candidate |
| `MinFragmentsModelling` | 1 | Minimum index-pass fragment matches for expectation-model inclusion |
| `ReportedHitsPerLabel` | 10 | Top hits to report separately for targets and decoys; affects runtime |
| `FragmentIndexBinWidth` | 0.01 | Index bin width, in Da |
| `nTerminalSeries` | `B` | N-terminal ion series used for classic rescoring |
| `cTerminalSeries` | `Y` | C-terminal ion series used for classic rescoring |
| `Andromeda` | `{ PMinPMax = 4, 10; MatchingIonTolerancePPM = 100.0 }` | Nested Andromeda-like scoring configuration |

**Note:** `nTerminalSeries` and `cTerminalSeries` only change the classic rescoring stage; the fragment index always contains b and y ions.

### Andromeda subparameters

`Andromeda` is a `ProteomIQon.Domain.AndromedaParams` record.

| Parameter | Default | Description |
| --- | --- | --- |
| `PMinPMax` | `4, 10` | Minimum and maximum count of strongest peaks per 100 Da window. Every count in the range is tested; the best result is chosen. |
| `MatchingIonTolerancePPM` | 100.0 | Theoretical-to-observed fragment matching tolerance, in ppm |

[Find the reference parameter JSON on GitHub ↗](https://github.com/CSBiology/ProteomIQon/search?q=peptideSpectrumMatchingTIMsParams.json)

## Outputs

Every input run generates `<run>.psm` in the output directory, as a header-containing, tab-delimited file. Existing PSM files are overwritten. Its data record is `Dto.PeptideSpectrumMatchingResult`, shared with PeptideSpectrumMatching.

The table contains 28 columns: 23 classic PSM columns and 5 TIMs-specific additions. `Label` is `1` for targets and `-1` for reversed decoys; `MissCleavages` is set to `-1` and recalculated by PSMStatistics.

### Classic PSM columns (23)

| Column | Description |
| --- | --- |
| PSMId | MS/MS spectrum match identifier |
| GlobalMod | Whether the peptide ion species is labeled or unlabeled |
| PepSequenceID | Unmodified peptide sequence identifier in PeptideDB |
| ModSequenceID | Modified peptide sequence identifier in PeptideDB, including PTMs |
| Label | Target (+1) or decoy (-1) |
| ScanNr | MS2 position in scan-time order |
| ScanTime | MS2 retention time in minutes |
| Charge | Precursor charge |
| PrecursorMZ | Precursor mass-to-charge ratio |
| TheoMass | Theoretical peptide mass in Da |
| AbsDeltaMass | Mass error using the best-matching isotope offset |
| PeptideLength | Number of peptide amino acids |
| MissCleavages | Value -1 at this stage; recomputed by PSMStatistics |
| SequestScore | SEQUEST-like similarity score |
| SequestNormDeltaBestToRest | Normalized separation of best score from remaining candidates |
| SequestNormDeltaNext | Normalized separation of top two candidate scores |
| AndroScore | Andromeda-like match score |
| AndroNormDeltaBestToRest | Normalized best-to-rest Andromeda score separation |
| AndroNormDeltaNext | Normalized top-two Andromeda score separation |
| XTandemScore | X!Tandem-like match score |
| XtandemNormDeltaBestToRest | Normalized best-to-rest X!Tandem score separation |
| XtandemNormDeltaNext | Normalized top-two X!Tandem score separation |
| StringSequence | Matching peptide amino acid sequence in one-letter code |

### TIMs-specific columns (5)

| Column | Description |
| --- | --- |
| IonMobility | Inverse reduced mobility (1/K0) of the precursor; `NaN` if unavailable |
| Hyperscore | Exact b/y-fragment hyperscore |
| Expectscore | Expectation value estimated using the index-pass hyperscore distribution |
| MatchedIons | Count of matched b and y fragments used in the hyperscore |
| TotalIons | Count of theoretical peptide b and y fragments |

### Identifiers and logs

`ScanNr` is the MS2 position after ordering spectra by scan time. For `PSMId`, spectrum-id spaces become hyphens and the suffix `_ScanNr_charge_rank` is added. Rank is determined by hyperscore within the reported target or decoy group for that spectrum.

The output folder also contains `PeptideSpectrumMatchingTIMs_log.txt` and one `<run>_log.txt` for each run.

**Downstream compatibility:** PSMStatistics can read these PSM results directly and uses `Hyperscore` and `Expectscore`, when present. Its `.qpsm` output does not yet retain ion mobility. Furthermore, [PSMBasedQuantificationTIMs ↗](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantificationTIMs.html) currently expects the FragPipe-style columns generated by [MsFraggerToPSM ↗](https://csbiology.github.io/ProteomIQon/tools/MsFraggerToPSM.html). Quantification after identification with this tool therefore needs adaptations to those downstream tools.

## Installing and executing the tool

1. **Installing**

   Install the Tool with the following command inside the choosen Command-Line

   ```sh
   dotnet tool install --global ProteomIQon.PeptideSpectrumMatchingTIMs
   ```

2. **Execute**

     ```sh
     proteomiqon-peptidespectrummatchingtims -i path/to/run.mzlite -d path/to/AraTest.db -o path/to/output -p path/to/peptideSpectrumMatchingTIMsParams.json
     ```

3. **Running the tool for multiple runs**

   ```sh
   proteomiqon-peptidespectrummatchingtims -i path/to/run1.mzlite path/to/run2.mzlite path/to/run3.mzlite -d path/to/AraTest.db -o path/to/output -p path/to/peptideSpectrumMatchingTIMsParams.json -c 3
   ```

4. **With a large database, build the index in four slices and keep the spectra in memory between them:**

   ```sh
   proteomiqon-peptidespectrummatchingtims -i path/to/run.mzlite -d path/to/AraTest.db -o path/to/output -p path/to/peptideSpectrumMatchingTIMsParams.json -s 4 -k
   ```

5. **All Flags**

   ```sh
   proteomiqon-peptidespectrummatchingtims --help
   ```
