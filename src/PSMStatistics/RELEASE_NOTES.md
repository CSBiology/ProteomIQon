#### Unreleased
* Add the parameter CrossValidation (default false): each scan is scored by the model trained on the other two of three scan folds, the scores are normalised per fold, and q-values and PEPs are computed on the pooled scores
* Give a decoy the pseudo-reversed sequence that PeptideSpectrumMatching and PeptideSpectrumMatchingTIMs score, without terminator tokens and with the flanking residues of its target, in the Peptide column of the training data instead of the reversed sequence. The column is not a model feature, so the results do not change
* Update BioFSharp.Mz to 0.2.3
* Cap the leaf output of the FastTree trees at 10 instead of 100, so that a few leaves no longer push near-empty spectra and decoys to the top of the ranking
* With CrossValidation, retrain each fold on the targets with q < 0.01 and accept a retrained model by its targets with q <= 0.01, independent of QValueThreshold and PepValueThreshold, which then only filter the output
* Train FastTree on one thread; the results do not change

#### 0.0.11 - Monday, September 28, 2026
* Use Hyperscore and Expectscore as features and seed the training with the expectation value when the input carries them
* Breaking: read the 28 column .psm layout of this release, .psm files of earlier releases have to be regenerated with PeptideSpectrumMatching
* Fix the PEPs: tied scores got a NaN PEP, and a few high scoring decoys could set almost every PEP to 1, so the PEP cutoff rejected many PSMs
* Update BioFSharp.Mz to 0.2.2

#### 0.0.10 - Wednesday, September 2, 2026
* Update to .NET 10
* Update BioFSharp to 2.0.0, BioFSharp.Mz to 0.2.1, FSharpAux to 2.1.0, FSharp.Stats to 0.6.0 and Plotly.NET to 6.0.0-preview.2

#### 0.0.9 - Monday, March 10, 2025
* Update to .NET 8

#### 0.0.8 - Thursday, August 2, 2022
* susbtitute PEP value calculation with iterative reweighted least squares based method

#### 0.0.7 - Friday, August 20, 2021
* rework fitting based PEP value calculation

#### 0.0.6 - Tuesday, July 13, 2021
* refined fitting based PEP value calculation

#### 0.0.5 - Friday, June 18, 2021
* add flag to trigger name based file matching
* add flag to trigger the creation of diagnostic charts

#### 0.0.4 - Friday, May 14, 2021
* update MzLite version to 0.1.1

#### 0.0.3 - Tuesday, May 5, 2021
* change relative file path handling

#### 0.0.2 - Tuesday, May 5, 2021
* Update PeptideSpectrumMatching Params
* Add tools

#### 0.0.1 - Tuesday, April 30, 2021
* Initial release
