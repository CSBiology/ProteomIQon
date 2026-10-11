(**
---
title: PSMStatistics
category: Tools
categoryindex: 1
index: 7
---
*)

(*** hide ***)

(*** condition: prepare ***)
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

# PSMStatistics

[PeptideSpectrumMatching]({{root}}tools/PeptideSpectrumMatching.html) and [PeptideSpectrumMatchingTIMs]({{root}}tools/PeptideSpectrumMatchingTIMs.html) report several target and decoy candidates per spectrum, each with a SEQUEST-like, an Andromeda-like and an X!Tandem-like score.
PSMStatistics combines these columns into one score with a FastTree classifier (ML.NET) trained to separate targets from decoys. It keeps the best candidate per scan by this score, computes its [q-value](https://www.biofsharp.com/BioFSharp.Mz/04_01_fdr_control.html#Computing-q-values-by-counting-decoys) and [PEP value](https://www.biofsharp.com/BioFSharp.Mz/04_01_fdr_control.html#From-list-level-q-to-single-PSM-PEP), and writes the targets that pass both thresholds.

The model features are the charge (one-hot encoded), the precursor m/z, the theoretical mass, the absolute mass error, the peptide length, the number of missed cleavages (looked up in the peptide database), and the three scores with their `NormDeltaBestToRest` and `NormDeltaNext` columns. When every row of the input has finite `Hyperscore` and `Expectscore` values, as PeptideSpectrumMatchingTIMs writes them, the `Hyperscore` and the negative log10 of the `Expectscore` are added. The tool reads the 28-column `.psm` layout; `.psm` files with the older 23-column layout have to be regenerated.

The training is semi-supervised. The first model starts from a single score: the negative log10 of the `Expectscore` when the input carries it, otherwise the SEQUEST-like score. Among the best candidates per scan by this score, the targets with q < 0.001 are the positives, and all decoy candidates of the other scans are the negatives. Each model is fitted on a random 90 % of these candidates; the other 10 % give the accuracy, AUC and F1 in the log. FastTree runs on one thread with the leaf output capped at 10. The tool then scores all candidates with the model, takes the best candidate per scan again and retrains on the new positives: the targets with q < 0.01 with `CrossValidation`, otherwise those with a q-value below a tenth of `QValueThreshold` and a PEP below `PepValueThreshold`. After each retraining it counts the passing scans as described for `MinimumIncreaseBetweenIterations`. If the count drops, or the retraining fails, the training ends with the previous model. If the count grows by less than `MinimumIncreaseBetweenIterations`, it ends with the new model. Otherwise it continues, for at most `MaxIterations` retrainings. ML.NET and the assignment of scans to folds use the seed 1024.

<img src="{{root}}img/SemiSupervisedScoring.png" width="1000" height="750" />

With `CrossValidation`, the scans are shuffled and split into three folds whose scan counts differ by at most one. The training runs three times, each time on the candidates of two folds, and the resulting model scores the remaining fold, so no scan is scored by a model that was trained on it. As in Percolator, the scores of each fold are then shifted and scaled so that the lowest scoring target with q <= 0.01 of that fold gets 0 and the median decoy gets -1. The q-values and PEPs are computed on the pooled scores of the three folds. Without `CrossValidation`, one model is trained on all scans and scores them.

The q-values and PEPs are computed on the best candidate per scan, targets and decoys. The q-value of a score is the number of decoys divided by the number of targets among the best candidates that score at least as high, made monotone (`FDRControl.calculateQValueStorey` of BioFSharp.Mz, without a π0 estimate). For the PEP, the scores are binned with a width of one unit of the FastTree score (with `CrossValidation`, that unit converted to the normalized scale and averaged over the three folds), and the logit of the decoy share per bin is fitted with a smoothing spline by iteratively reweighted least squares, as in Percolator. The PEP of a target is the fitted decoy to target ratio at its score, capped at 1 and made non-increasing with the score.

The tool writes the best candidate of a scan when it is a target with a q-value below `QValueThreshold` and a PEP below `PepValueThreshold`. If no candidate passes, the `.qpsm` file holds only the header. If the estimation throws, the file also holds only the header and `run_log.txt` holds the exception.

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | one or more `.psm` files, or a directory that is searched for `*.psm` (not recursive) | [PeptideSpectrumMatching]({{root}}tools/PeptideSpectrumMatching.html), [PeptideSpectrumMatchingTIMs]({{root}}tools/PeptideSpectrumMatchingTIMs.html) |
| `-d` | the SQLite peptide database | [PeptideDB]({{root}}tools/PeptideDB.html) |
| `-o` | the output directory, created when missing | |
| `-p` | the parameter file in JSON | this page |
| `-dc` | switch: save diagnostic charts to the output directory | |
| `-c` | number of files scored at the same time, default 1 | |

For every input `run.psm` the tool writes `run.qpsm` to the output directory. The file is tab separated with a header and holds one row per scan: the identifiers of the peptide and its modification state, the search scores, the combined `ModelScore`, `QValue`, `PEPValue`, the sequence and the protein names.
[ProteinInference]({{root}}tools/ProteinInference.html), [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html) and [PSMBasedQuantificationTIMs]({{root}}tools/PSMBasedQuantificationTIMs.html) read `.qpsm` files.

The tool also creates a `run_plots` directory per input. It stays empty unless you pass `-dc`. Each training then writes `InitialSeparation.html` (target and decoy histograms of the starting score with its q-values), `Metrics.html` (the metrics on the 10 % test split of the last model fitted) and one `separationAtIteration_<n>.html` per retraining, counted from 0 (the chart of the first model is overwritten by retraining 0; target and decoy histograms of the model score with the q-value and PEP curves). With `CrossValidation` the names start with `Fold1_` to `Fold3_`, and the charts show the two training folds of that model. `Threshold.Fixed` writes no charts. Log output goes to `PSMStatistics_log.txt` and `run_log.txt` in the output directory.

## Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `Threshold` | `Threshold.Estimate estimationParams` | Learn the combined score and filter by q-value and PEP value. The nested record is listed below. |
| `Threshold` (alternative) | `Threshold.Fixed { SequestLike = 0.; Andromeda = 0. }` | Skip the learning. Among the matches whose `SequestScore` exceeds `SequestLike` and whose `AndroScore` exceeds `Andromeda`, keep the one with the highest `SequestScore` per scan, and drop it when it is a decoy. `ModelScore`, `QValue` and `PEPValue` are written as `nan`. Pick your own cutoffs; these values are placeholders. |
| `ParseProteinIDRegexPattern` | `"id"` | Regex applied to the protein names of the database to produce the `ProteinNames` column. `"id"` keeps them unchanged. The database holds the names as PeptideDB cut them, so keep `"id"`. Another pattern is applied to those stored names and keeps the matched part (empty when it does not match). |
| `KeepTemporaryFiles` | `true` | Not used by PSMStatistics. |

The fields of `Threshold.Estimate`:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `QValueThreshold` | `0.01` | Write PSMs with a q-value below this value. Without `CrossValidation`, the positives of each retraining are the targets with a q-value below a tenth of this value and a PEP below `PepValueThreshold`. |
| `PepValueThreshold` | `1.0` | Write PSMs with a PEP below this value. The default removes only PSMs with a PEP of 1. |
| `MaxIterations` | `15` | Upper limit for the retrainings after the first model, per fold with `CrossValidation`. |
| `MinimumIncreaseBetweenIterations` | `0.005` | Stop when a retraining raises the number of passing scans by less than this fraction of the new number. With `CrossValidation` this number counts the targets with q <= 0.01. Without it, it is the smaller of the numbers of best candidates per scan with q <= `QValueThreshold` and with PEP <= `PepValueThreshold`. |
| `PepValueFittingMethod` | `PepValueFittingMethod.IRLS` | How the PEP curve is fitted. `IRLS` (iteratively reweighted least squares) is the only case. |
| `CrossValidation` | `true` | Score each scan with the model trained on the other two of three scan folds, and normalize the scores per fold. Each fold retrains on the targets with q < 0.01, independent of `QValueThreshold` and `PepValueThreshold`, which then only filter the output. The training runs three times per file, each time on two thirds of the scans, so a run takes about two to three times as long as without `CrossValidation`. Parameter files without this field use `false`. |

The default file is [pSMStatisticsParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/pSMStatisticsParams.json).

## Writing a parameter file
*)

open ProteomIQon
open ProteomIQon.Domain

let psmStatisticsParams : Dto.PSMStatisticsParams =
    {
        Threshold =
            Threshold.Estimate
                {
                    QValueThreshold                  = 0.01
                    PepValueThreshold                = 1.0
                    MaxIterations                    = 15
                    MinimumIncreaseBetweenIterations = 0.005
                    PepValueFittingMethod            = PepValueFittingMethod.IRLS
                    CrossValidation                  = true
                }
        ParseProteinIDRegexPattern = "id"
        KeepTemporaryFiles         = true
    }

// Replace the temp folder with your project folder.
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "pSMStatisticsParams.json")

Json.serializeAndWrite outputPath psmStatisticsParams

(**
## Outputs
Columns of the `.qpsm` file:

| Column                     | Description                                                                                                      |
|----------------------------|------------------------------------------------------------------------------------------------------------------|
| PSMId                      | Identifier of the MS/MS spectrum                                                                                 |
| GlobalMod                  | 0 for an unlabeled, 1 for a labeled peptide ion species                                                        |
| PepSequenceID              | Identifier of the unmodified peptide sequence in the peptide database                                            |
| ModSequenceID              | Identifier of the modified peptide sequence in the peptide database                                              |
| Label                      | 1 (target); decoys are not written                                                                               |
| ScanNr                     | Number of the MS2 spectrum within the run                                                                        |
| ScanTime                   | Retention time of the MS/MS scan in minutes                                                                      |
| Charge                     | Precursor charge                                                                                                 |
| PrecursorMZ                | Precursor m/z                                                                                                    |
| TheoMass                   | Theoretical mass of the peptide in Da                                                                            |
| AbsDeltaMass               | Absolute difference between theoretical and measured mass                                                        |
| PeptideLength              | Number of amino acids                                                                                            |
| MissCleavages              | Number of missed cleavages, looked up in the peptide database                                                    |
| SequestScore               | SEQUEST-like score                                                                                               |
| SequestNormDeltaBestToRest | (best SEQUEST score of the spectrum - this score) / best score; 0 for the top candidate                          |
| SequestNormDeltaNext       | (this SEQUEST score - next lower SEQUEST score of the spectrum) / best score                                     |
| AndroScore                 | Andromeda-like score                                                                                             |
| AndroNormDeltaBestToRest   | (best Andromeda score of the spectrum - this score) / best score; 0 for the top candidate                        |
| AndroNormDeltaNext         | (this Andromeda score - next lower Andromeda score of the spectrum) / best score                                 |
| XtandemScore               | X!Tandem-like score                                                                                              |
| XtandemNormDeltaBestToRest | (best X!Tandem score of the spectrum - this score) / best score; 0 for the top candidate                         |
| XtandemNormDeltaNext       | (this X!Tandem score - next lower X!Tandem score of the spectrum) / best score                                   |
| ModelScore                 | Score of the learned model; with `CrossValidation` the fold-normalized score (0 at the lowest-scoring target with q <= 0.01 of its fold, -1 at the median decoy of the fold); without it the FastTree score |
| QValue                     | q-value of the `ModelScore`                                                                                      |
| PEPValue                   | Posterior error probability of the `ModelScore`                                                                  |
| StringSequence             | Peptide sequence of the match, as in the `.psm` input                                                            |
| ProteinNames               | Names of the proteins that contain the peptide, separated by `;`                                                 |

*)


(**
## Running the tool

Install with `dotnet tool install --global ProteomIQon.PSMStatistics`, then score one run:

```text
proteomiqon-psmstatistics -i path/to/run.psm -d path/to/database.db -o path/to/output -p path/to/pSMStatisticsParams.json
```

Several runs at once, three of them in parallel:

```text
proteomiqon-psmstatistics -i path/to/run1.psm path/to/run2.psm path/to/run3.psm -d path/to/database.db -o path/to/output -p path/to/pSMStatisticsParams.json -c 3
```

Add `-dc` to write the separation charts for each iteration:

```text
proteomiqon-psmstatistics -i path/to/psmFolder -d path/to/database.db -o path/to/output -p path/to/pSMStatisticsParams.json -c 3 -dc
```

All flags:

```text
proteomiqon-psmstatistics --help
```
*)
