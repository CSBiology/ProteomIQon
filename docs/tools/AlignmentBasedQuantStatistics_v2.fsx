(**
---
title: AlignmentBasedQuantStatistics
category: Tools
categoryindex: 1
index: 13
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

# AlignmentBasedQuantStatistics

A peptide ion transferred by [QuantBasedAlignment]({{root}}tools/QuantBasedAlignment_v2.html) and quantified by
[AlignmentBasedQuantification]({{root}}tools/AlignmentBasedQuantification_v2.html) has no identification in the run, so nothing yet says
whether the quantified peak is the ion. AlignmentBasedQuantStatistics gives every such ion two values:

- `AlignmentScore`: the probability that the quantified peak is the ion's peak, given that the analyte is present in the run.
- `AlignmentQValue`: the estimated share of false quantifications among the new ions of the run that score at least as well. A quantification is
  false when the chosen peak is not the ion's peak or, unless `AbsentAnalyteFdr` is `Off`, when the analyte is absent from the run.

New ions are the transferred ions without a PSM based quantification in the run.

[JoinQuantPepIonsWithProteins]({{root}}tools/JoinQuantPepIonsWithProteins.html) reads the output.
[LabelFreeProteinQuantification]({{root}}tools/LabelFreeProteinQuantification.html) and
[LabeledProteinQuantification]({{root}}tools/LabeledProteinQuantification.html) filter on the q-value only if their `Alignment_QValue` parameter
is set. They then drop transferred rows with a q-value at or above it. The default `None` keeps every transferred row.

<img src="{{root}}img/AlignmentV2Pipeline.png" width="1000" />

## Training data and label

The tool learns from overlap ions: ions that a run identified and quantified itself and that were also transferred to it. For them the run's
own PSM based peak is the reference. Let z be the distance between the transferred peak and the own peak, in units of the own peak's standard
deviation.

| Label | Condition | Use |
|-------|-----------|-----|
| right | z ≤ `PositiveRtSdFactor` (1) and &#124;log2 (transferred / own quantity)&#124; ≤ `MaxAbsLog2QuantRatio` (1) | training |
| ambiguous | `PositiveRtSdFactor` < z ≤ `NegativeRtSdFactor` (3), same quantity bound | left out of training |
| wrong | everything else | training |

Overlap ions whose own peak cannot serve as reference (no fitted peak, or a width or quantity that is not positive) get no label.

## Features

All features come from the transferred row of AlignmentBasedQuantification and the `.align` row. None uses the ion's own PSM row, so the
features of an overlap ion look like those of a new ion.

| Group | Features |
|-------|----------|
| scan time | residual to `PredictedScanTime` (signed, absolute, per peak sd, per local alignment error), the local alignment error, peak width and its ratio to the source peak |
| mass and isotopes | &#124;ppm&#124; of the monoisotopic peak, intensity weighted &#124;ppm&#124; over the isotope peaks, number of isotope peaks, similarity of the isotope pattern to the source, KL divergences to the theoretical pattern |
| chromatogram | correlation with the source trace (around the peak and over the whole XIC), DTW distance to the source trace, competing peaks in the XIC, share of the XIC maximum at the apex, XIC length, relative SEO |
| ion | source quantity and apex intensity, source scan time, m/z, charge |
| selection | from `<run>.alignfeatures`: `Sel_NCandidates`, `Sel_IsoCorr`, `Sel_IsoRatioLog2`, `Sel_ShapeCorr`, `Sel_RtZ`, `Sel_ImZ`, `Sel_LogRelHeight`, `Sel_AbsLog2FwhmRatio`, and `Sel_CandidateProbability`, `Sel_CandidateMargin` when all runs used the same selection model |
| ion mobility | &#124;`IonMobility` − `PredictedIonMobility`&#124; / `PredictedIonMobilitySd` (timsTOF only) |

The local alignment error of an ion is the median |own PSM scan time − `PredictedScanTime`| of the run's overlap ions predicted within ±2 min.
Every run needs `<run>.alignfeatures` next to its aligned `.quant`.

`UseIntensityRatioFeatures` adds four features: the ratios of target to source quantity and apex intensity, and the target quantity and apex
intensity. They are off by default: a model that prefers a ratio near 1 penalises peptides that really changed between the samples.

## Model, calibration and q-values

1. The overlap ions of all learning runs are pooled and split into 5 folds by peptide sequence, so the forms of a peptide stay together.
2. Per fold, a gradient boosted tree model (ML.NET FastTree) is trained on the right and wrong ions of the other four folds. Unless
   `AbsentAnalyteFdr` is `Off`, the decoys of the learning runs join this training as wrong examples, in the fold of their peptide.
3. Per fold, a Platt calibration (a logistic fit of the label on the raw score) is fitted on the right and wrong overlap ions of the held-out
   fold. It turns scores into probabilities.
4. The `AlignmentScore` p of a new ion is the mean of the calibrated probabilities of the 5 fold models.
5. The posterior error probability of a new ion is 1 − p. Per run, the new ions are sorted by p in descending order. The `AlignmentQValue` of
   an ion is the mean posterior error probability of all new ions with at least its p, made monotone.
6. Unless `AbsentAnalyteFdr` is `Off`, the absent-analyte rate at the ion's score (see Absent analytes) is added to the q-value, which is made
   monotone again.

The tree model rarely gives probabilities above 0.99, so a cutoff of q ≤ 0.01 can keep few or no transferred ions.

The log reports an out-of-fold check. Each overlap ion is scored by the fold model that did not see it and calibrated on the other folds.
The true share of wrong peaks is given at q ≤ 0.01, 0.05 and 0.1, once with the ambiguous ions counted as wrong and once without them.

## Absent analytes

The overlap ions are present in the run, so they cannot show how the model scores an ion whose analyte is absent. The decoy transfers
of AlignmentBasedQuantification (`<run>.aligndecoy`) are used for that. Per run, the absent-analyte rate at score t is

FDR_absent(t) = π · (decoys with p ≥ t / all decoys) · N_new / N_new(p ≥ t)

with N_new the number of new ions of the run and π the fraction of new ions whose analyte is absent. `Estimated` estimates π per run as
Storey does, with λ = 0.5. `Conservative` sets π to 1.

- A run with fewer than 200 decoys uses the decoys of all runs to score.
- Decoys without a fitted peak and new ions without an alignment row count as p = 0.
- When a run to score or to learn from has no `<run>.aligndecoy` next to its aligned `.quant` (AlignmentBasedQuantification writes none
  with `DecoyFraction` 0), or no decoy has a fitted peak, the tool warns and works as `Off`.

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | PSM based `.quant` of the runs to score | [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html) or [PSMBasedQuantificationTIMs]({{root}}tools/PSMBasedQuantificationTIMs.html) |
| `-ii` | `.align` of the runs to score | [QuantBasedAlignment]({{root}}tools/QuantBasedAlignment_v2.html) |
| `-iii` | aligned `.quant` of the runs to score | [AlignmentBasedQuantification]({{root}}tools/AlignmentBasedQuantification_v2.html) |
| `-l` | PSM based `.quant` of the runs to learn from | as `-i` |
| `-ll` | `.align` of the runs to learn from | as `-ii` |
| `-lll` | aligned `.quant` of the runs to learn from | as `-iii` |
| `-o` | output directory, created when missing | |
| `-p` | parameter file (JSON) | this page |
| `-mf` | switch: pair the files by name instead of by position | |
| `-c` | number of runs read and written in parallel, one core each, default 1 | |
| `-ts` | switch: write the training ions of each learning run | |
| `-dc` | switch: write `AlignmentScoreHistogram.html` | |

Each flag takes one file, several files, or a directory. The runs to learn from are usually the runs to score. Without `-mf` the files are
paired by position, the files of a directory sorted by name, and the numbers of files must agree.

For every scored run the tool writes `<run>.quant` with the layout of the input `.quant` files: the PSM based rows of `-i` and one row per new
ion with `AlignmentScore` and `AlignmentQValue`. When the PSM row of an overlap ion has no fitted peak, the transferred row replaces it and is
scored like a new ion. `-ts` writes `<run>.quanttestset` per learning run, with the overlap ions, their label, fold, raw and calibrated
out-of-fold score, out-of-fold q-value and features. Logs go to `AlignmentBasedQuantStatistics_log.txt`, with the features, the out-of-fold
check and the calibration per score bin and per source intensity quartile, and to one `<run>_log.txt` per run.

## Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| PositiveRtSdFactor | 1.0 | A peak within this many own peak sd is right. |
| NegativeRtSdFactor | 3.0 | Up to this many sd a peak is ambiguous, beyond it wrong. |
| MaxAbsLog2QuantRatio | 1.0 | Right and ambiguous peaks need &#124;log2 (transferred / own quantity)&#124; at most this. |
| UseIntensityRatioFeatures | false | Add target to source intensity ratios as features. |
| UseIonMobility | true | Add the ion mobility feature when the rows have a finite `IonMobility` and the `.align` file a `PredictedIonMobilitySd` above 0. |
| NumberOfFolds | 5 | Cross-validation folds. |
| Trainer | FastTree | `FastTree` or `LightGbm`. LightGbm needs the native library lib_lightgbm, which the package ships for win-x64, linux-x64 and osx-x64. |
| NumberOfLeaves, NumberOfTrees, MinimumExampleCountPerLeaf, LearningRate | null | null takes the trainer's defaults: 20, 100, 10, 0.2 for FastTree; 31, 100, 20, 0.1 for LightGbm. |
| AbsentAnalyteFdr | Estimated | `Off`, `Estimated` or `Conservative`: whether and how the absent-analyte rate enters the q-value (see Absent analytes). |
| Seed | 1 | Seed of the folds and of the trainer. |

See [the default parameter file](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/AlignmentBasedQuantStatisticsParams_v2.json).

## Writing a parameter file
*)

open ProteomIQon

let alignmentBasedQuantStatisticsParams : Dto.AlignmentBasedQuantStatisticsV2Params =
    {
        PositiveRtSdFactor         = 1.
        NegativeRtSdFactor         = 3.
        MaxAbsLog2QuantRatio       = 1.
        UseIntensityRatioFeatures  = false
        UseIonMobility             = true
        NumberOfFolds              = 5
        Trainer                    = "FastTree"
        // empty: the defaults of the trainer
        NumberOfLeaves             = System.Nullable()
        NumberOfTrees              = System.Nullable()
        MinimumExampleCountPerLeaf = System.Nullable()
        LearningRate               = System.Nullable()
        AbsentAnalyteFdr           = "Estimated"
        Seed                       = 1
    }

// Replace the temp folder with your project folder.
let outputPath =
    System.IO.Path.Combine(System.IO.Path.GetTempPath(), "AlignmentBasedQuantStatisticsParams.json")

Json.serializeAndWrite outputPath alignmentBasedQuantStatisticsParams

(**
## Running the tool

Install the tool with `dotnet tool install --global ProteomIQon.AlignmentBasedQuantStatistics`.
Score every run, learn from all of them, pair the files by name, and work on four runs at a time:

```text
proteomiqon-alignmentbasedquantstatistics -i path/to/quant -ii path/to/align -iii path/to/aligned -l path/to/quant -ll path/to/align -lll path/to/aligned -o path/to/output -p path/to/AlignmentBasedQuantStatisticsParams.json -mf -c 4
```

Add `-ts` to keep the training ions and `-dc` for the score histogram. Print the description of every argument:

```text
proteomiqon-alignmentbasedquantstatistics --help
```

## Limitations

- The false discovery rate of absent analytes is an estimate from decoys at a shifted scan time. Analytes near the detection limit are not
  modelled.
- The model needs overlap ions. With few identifications per run, or runs that share few ions, there is little to learn from.
*)
