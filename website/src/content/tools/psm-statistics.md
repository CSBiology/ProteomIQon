---
name: "PSMStatistics"
phase: peptide-identification
subgroup: PSM statistics
order: 3
summary: "PSMStatistics learns one combined score from the per-spectrum scores reported by PeptideSpectrumMatching."
inputsFrom: [peptide-spectrum-matching, peptide-db]
---

## Tool Description

   [PeptideSpectrumMatching](./peptide-spectrum-matching.md) and [PeptideSpectrumMatchingTIMs](./peptide-spectrum-matching-tims.md) report several target and decoy candidates per spectrum, each with a `SEQUEST-like`, an `Andromeda-like` and an `X!Tandem-like` score. PSMStatistics combines these columns into one score with a `FastTree classifier (ML.NET)` trained to separate targets from decoys. It keeps the best candidate per scan by this score, computes its Q-value and PEP value, and writes the targets that pass both thresholds.

## How does PSMStats work?

1. **Evaluate peptide matches**

   Each MS/MS spectrum may have multiple peptide candidates, which are evaluated using three scoring methods: SEQUEST, Andromeda and X!Tandem.
   The tool also considers additional peptide characteristics, including precursor charge, mass accuracy, peptide length and missed cleavages.
   If available, Hyperscore and Expectscore provide additional scoring information. Retention time can also be included to improve the evaluation.

2. **Train a scoring model**

   PSMStatistics uses a FastTree classifier (ML.NET) to distinguish between: 
      - Targets: Matches to real peptide sequences.
      - Decoys: Matches to artificial sequences used to estimate incorrect identifications.
   
   The initial training uses highly confident target matches and decoys to teach the model which characteristics indicate a reliable identification.
   The model combines these characteristics into a single ModelScore for each peptide candidate.
   Training is repeated with updated high-confidence matches until the results no longer improve sufficiently or the maximum number of iterations is reached.

3. **Validate the model**
   By default, PSMStatistics uses three-fold cross-validation to reduce overfitting and improve the reliability of the scoring model.
   The spectra are divided into three groups. For each iteration of cross-validation, the model is trained using two groups and evaluates the remaining group.
   This ensures that each spectrum is evaluated by a model that was not trained on that spectrum.
   The resulting scores are normalized across the three groups to make them comparable.

4. **Estimate identification reliability**
   The tool selects the best match for each spectrum and calculates two values:
   - QValue: Estimates the proportion of incorrect identifications among matches accepted at a given score threshold.
   -  PEPValue: Estimates the probability that an individual peptide match is incorrect.
   Both values are calculated using the score distributions of target and decoy matches.
   Lower Q-Values and PEP-Values indicate more reliable peptide identifications.

5. **Filter and save results**
   The tool retains only the highest-scoring target match per spectrum if both its Q-Value and PEP-Value are below the specified thresholds.
   The filtered results are saved as `.qpsm `files, containing peptide information, search scores, the ModelScore and the estimated confidence values.
   These files can be used for subsequent protein inference and quantification.
   Alternative scoring method: Instead of training a model, PSMStatistics can apply fixed SEQUEST and Andromeda score thresholds. In this case, no Q-Values or PEP-Values are calculated.

![Protein Infrence](../../assets/images/SemiSupervisedScoring.png)

## Inputs

These inputs are needed to execute the tool:

| Flag | Meaning | Comes from |
| --- | --- | --- |
| `-i` | one or more `.psm` files, or a directory that is searched for `*.psm` (not recursive) | [PSMStatistics](./peptide-spectrum-matching.md) |
| `-d` | The Peptide Database | [PeptideDB](./peptide-db.md) |
| `-o` | the output directory, created when missing | - |
| `-p` | the parameter file in JSON | [Github ProteomIQon Core](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/pSMStatisticsParams.json) |
| `-dc` | switch: save diagnostic charts to the output directory |  |
| `-c` | number of files scored at the same time, default 1 | - |

- For every input `run.psm` the tool writes `run.qpsm` to the output directory.
- The file is tab separated with a header and holds one row per scan.
- The tool also creates a `run_plots` directory per input. It stays empty unless you pass `-dc`, which adds `Metrics.html` and `InitialSeparation.html`, plus one `separationAtIteration_n.html` per training iteration.

## Parameter file

PSMStatistics can either estimate a combined score using iterative learning or apply fixed score cutoffs. The default uses `Threshold.Estimate`.

### Threshold selection

| Parameter | Default | Meaning |
| --- | --- | --- |
| `Threshold` | `Threshold.Estimate estimationParams` | Learn a combined model score and retain PSMs according to q-value and PEP thresholds. The nested estimation settings are listed below. |
| `Threshold` (alternative) | `Threshold.Fixed { SequestLike = 0.; Andromeda = 0. }` | Skip model training and keep the best match per scan only if its `SequestScore` exceeds `SequestLike` and its `AndroScore` exceeds `Andromeda`. Decoy matches are removed. The resulting `ModelScore`, `QValue` and `PEPValue` are `nan`. Both zero cutoffs shown here are placeholders; choose values appropriate for your data. |
| ParseProteinIDRegexPattern | `id` | Regex applied to the protein names of the database to produce the ProteinNames column. "id" keeps them unchanged. The database holds the names as PeptideDB cut them, so keep "id". Another pattern is applied to those stored names and keeps the matched part (empty when it does not match). |
| StoreIsotopePatterns | `true` |  Store the theoretical isotope patterns of the ions in the .`.qpsm` files in the table IsotopePattern of the peptide database (-d), where the quantification tools read them. |
| KeppTemporaryFiles | `true` | Not used by PSMStatistics.

### Estimation parameters

These settings belong to the nested `Threshold.Estimate` record.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `QValueThreshold` | `0.01` | Retain PSMs whose q-value is below this threshold. |
| `PepValueThreshold` | `0.05` | Retain PSMs whose posterior error probability (PEP) is below this threshold. |
| `MaxIterations` | `15` | Maximum number of iterations used to retrain the scoring model. |
| `MinimumIncreaseBetweenIterations` | `0.005` | Stop retraining when the relative increase in positive matches at the q-value threshold is smaller than this fraction. |
| `PepValueFittingMethod` | `PepValueFittingMethod.IRLS` | Fit the PEP curve using iteratively reweighted least squares (`IRLS`); this is the only available fitting method. |
| CrossValidation | `true` | Score each scan with the model trained on the other two of three scan folds, and normalize the scores per fold. Each fold retrains on the targets with q < 0.01, independent of `QValueThreshold` and `PepValueThreshold`, which then only filter the output. The training runs three times per file, each time on two thirds of the scans, so a run takes about two to three times as long as without `CrossValidation`. Parameter files without this field use `false`. |
| RetentionTimeFeature | `true` |  	Add the features `RtDeviation` and `AbsRtDeviation`, the deviation of the scan time from the time predicted from the sequence. Parameter files without this field use `false`. |

### General settings

| Parameter | Default | Meaning |
| --- | --- | --- |
| `ParseProteinIDRegexPattern` | `"id"` | Regular expression used on database FASTA headers to populate `ProteinNames`. Use the same identifier pattern as in PeptideDB. |
| `KeepTemporaryFiles` | `true` | Field in the parameter record. Leave this setting at `true`. |

The default parameter file is [pSMStatisticsParams.json ↗](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/pSMStatisticsParams.json).

## Outputs

The output `.qpsm` file contains 27 columns. The columns are grouped by category below.

### Identification and protein assignment (6)

PSM identifiers, peptide sequences and corresponding proteins

| Column | Description |
| --- | --- |
| `PSMId` | Identifier of the MS/MS spectrum |
| `GlobalMod` | Indicator, if a peptide ion species is labeled or unlabeled |
| `PepSequenceID` | Unique identifier of the unmodified peptide sequence, which points to PeptideDB |
| `ModSequenceID` | Unique identifier of the modified peptide sequence (including PTMs, e.g., methylation), which points to PeptideDB |
| `StringSequence` | Amino Acid sequence (one-letter code) from PeptideDB which matches the psm candidate |
| `ProteinNames` | Protein Names out of the PeptideDB |

### Spectrum, precursor and peptide properties (8)

Scan details, precursor masses, peptide length and missed cleavages

| Column | Description |
| --- | --- |
| `ScanNr` | Scan identifier, combining the spectrum ID in the raw file with an ascending MS2 ID |
| `ScanTime` | Retention time (RT) in minutes of the MS/MS scan |
| `Charge` | Precursor ion charge state |
| `PrecursorMZ` | Precursor ion mass-to-charge ratio (m/z) |
| `TheoMass` | Theoretical peptide mass in the spectrum (based on amino acid composition) in Dalton |
| `AbsDeltaMass` | Absolute mass deviation between theoretical and measured mass (mass error) |
| `PeptideLength` | Peptide length in Amino Acid count |
| `MissCleavages` | Number of missed cleavages |

### SEQUEST scoring (3)

Similarity score and normalized differences between candidates

| Column | Description |
| --- | --- |
| `SequestScore` | SEQUEST similarity score (e.g., XCorr) quantifying agreement between theoretical and experimental spectra |
| `SequestNormDeltaBestToRest` | Normalized separation of the best SEQUEST score from the remaining candidate scores |
| `SequestNormDeltaNext` | Normalized separation between the best and second-best SEQUEST scores |

### Andromeda scoring (3)

Match score and normalized differences between candidates

| Column | Description |
| --- | --- |
| `AndroScore` | Andromeda score quantifying the match between theoretical and experimental spectra |
| `AndroNormDeltaBestToRest` | Normalized separation of the best Andromeda score from the remaining candidate scores |
| `AndroNormDeltaNext` | Normalized separation between the best and second-best Andromeda scores |

### X!Tandem scoring (3)

Match score and normalized differences between candidates

| Column | Description |
| --- | --- |
| `XTandemScore` | XTandem score quantifying the match between theoretical and experimental spectra |
| `XtandemNormDeltaBestToRest` | Normalized separation of the best XTandem score from the remaining candidate scores |
| `XtandemNormDeltaNext` | Normalized separation between the best and second-best XTandem scores |

### Target/decoy and statistical confidence (4)

Model score, FDR and posterior error probability

| Column | Description |
| --- | --- |
| `Label` | Target/decoy label: 1 = target, −1 = decoy |
| `ModelScore` | Best score value achieved by iterative model to distinguish between target & decoy |
| `QValue` | Q-value (False-Discovery-Rate) based on combined model scores |
| `PEPValue` | Posterior Error Probability |

## Installing and executing the tool

1. **Installing**

   Install the Tool with the following command inside the choosen Command-Line

   ```sh
   dotnet tool install --global ProteomIQon.ProteinInference
   ```

2. **Execute**
     ```sh
     proteomiqon-proteininference -i path/to/run.qpsm -d path/to/database.db -o path/to/output -p path/to/ProteinInferenceParams.json
     ```

3. **Several runs grouped together, with the annotation for the evidence class**

    ```sh
     proteomiqon-proteininference -i path/to/run1.qpsm path/to/run2.qpsm path/to/run3.qpsm -d path/to/database.db -g path/to/proteome.gff3 -o path/to/output -p path/to/ProteinInferenceParams.json
     ```

4. **A whole directory, with the q-value chart**

   ```sh
     proteomiqon-proteininference -i path/to/qpsmFolder -d path/to/database.db -o path/to/output -p path/to/ProteinInferenceParams.json -dc
     ```

5. **All Flags**

    ```sh
     proteomiqon-proteininference --help
     ```
