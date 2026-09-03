(**
---
title: PSMStatistics
category: Tools
categoryindex: 1
index: 6
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

[PeptideSpectrumMatching]({{root}}tools/PeptideSpectrumMatching.html) reports several scores for every spectrum (SEQUEST-like, Andromeda-like, X!Tandem-like, plus mass error, peptide length and how far the best candidate is ahead of the next one).
PSMStatistics learns one combined score from these columns. It uses the target and decoy labels of the matches as training signal, retrains in iterations while the set of confident targets grows (semi supervised), and stops when an iteration adds too few new positives or the iteration limit is reached.
From the combined score it computes a [q-value](https://www.biofsharp.com/BioFSharp.Mz/04_01_fdr_control.html#Computing-q-values-by-counting-decoys) and a [PEP value](https://www.biofsharp.com/BioFSharp.Mz/04_01_fdr_control.html#From-list-level-q-to-single-PSM-PEP) per match and keeps the PSMs under both thresholds.

<img src="{{root}}img/SemiSupervisedScoring.png" width="1000" height="750" />

## Inputs and outputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | one or more `.psm` files, or a directory that is searched for `*.psm` (not recursive) | [PeptideSpectrumMatching]({{root}}tools/PeptideSpectrumMatching.html) |
| `-d` | the SQLite peptide database | [PeptideDB]({{root}}tools/PeptideDB.html) |
| `-o` | the output directory, created when missing | |
| `-p` | the parameter file in JSON | this page |
| `-dc` | switch: save diagnostic charts to the output directory | |
| `-c` | number of files scored at the same time, default 1 | |

For every input `run.psm` the tool writes `run.qpsm` to the output directory. The file is tab separated with a header and holds one row per scan: the identifiers of the peptide and its modification state, the search scores, the combined `ModelScore`, `QValue`, `PEPValue`, the sequence and the protein names.
[ProteinInference]({{root}}tools/ProteinInference.html) and [PSMBasedQuantification]({{root}}tools/PSMBasedQuantification.html) read `.qpsm` files.

The tool also creates a `run_plots` directory per input. It stays empty unless you pass `-dc`, which adds `Metrics.html` and `InitialSeparation.html`, plus one `separationAtIteration_<n>.html` per training iteration. Charts exist only for the estimated threshold. Log output goes to `PSMStatistics_log.txt` and `run_log.txt` in the output directory.

## Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `Threshold` | `Threshold.Estimate estimationParams` | Learn the combined score and filter by q-value and PEP value. The nested record is listed below. |
| `Threshold` (alternative) | `Threshold.Fixed { SequestLike = 0.; Andromeda = 0. }` | Skip the learning. Keep the best match per scan among those whose `SequestScore` exceeds `SequestLike` and whose `AndroScore` exceeds `Andromeda`, and drop it when it is a decoy. `ModelScore`, `QValue` and `PEPValue` are written as `nan`. Pick your own cutoffs, the values here are placeholders. |
| `ParseProteinIDRegexPattern` | `"id"` | Regex applied to the FASTA headers of the database to produce the `ProteinNames` column. Use the same pattern as in PeptideDB. |
| `KeepTemporaryFiles` | `true` | Part of the parameter record. Keep it at `true`. |

The fields of `Threshold.Estimate`:

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `QValueThreshold` | `0.01` | Keep PSMs with a q-value below this value. |
| `PepValueThreshold` | `0.05` | Keep PSMs with a PEP value below this value. |
| `MaxIterations` | `15` | Upper limit for the retraining iterations. |
| `MinimumIncreaseBetweenIterations` | `0.005` | Stop when the number of positives at the q-value threshold grows by less than this fraction from one iteration to the next. |
| `PepValueFittingMethod` | `PepValueFittingMethod.IRLS` | How the PEP curve is fitted. `IRLS` (iteratively reweighted least squares) is the only case. |

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
                    PepValueThreshold                = 0.05
                    MaxIterations                    = 15
                    MinimumIncreaseBetweenIterations = 0.005
                    PepValueFittingMethod            = PepValueFittingMethod.IRLS
                }
        ParseProteinIDRegexPattern = "id"
        KeepTemporaryFiles         = true
    }

// Replace the temp folder with your project folder.
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "pSMStatisticsParams.json")

Json.serializeAndWrite outputPath psmStatisticsParams

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
