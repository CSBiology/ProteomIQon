#### Unreleased
* Label overlap ions by their peak: positive within PositiveRtSdFactor sd of the own PSM peak and within MaxAbsLog2QuantRatio of its quant, ambiguous up to NegativeRtSdFactor sd (left out of training), else negative; overlap ions without truth are left out of training; their alignment row replaces the PSM row if that has no fitted scan time and is scored like a new ion, else it is not written
* Score with 26 features from the alignment row and the .align row: scan time residuals, isotope errors and pattern similarity, trace correlations with the source, competing peaks, source intensity; optional intensity ratio features and ion mobility residual
* Add the selection features of a <run>.alignfeatures file next to the aligned .quant as Sel_ features; CandidateProbability and CandidateMargin only when all runs used the same selection model
* Train one ML.NET FastTree model per cross-validation fold on the overlap ions of all learning runs, with folds by a seeded permutation of the peptides (PepSequenceID), and calibrate each fold on its held-out ions (Platt), on the held-out ions of all folds when its own are not of both classes
* Trainer selects the ML.NET model, FastTree (default) or LightGbm
* AlignmentScore is the mean calibrated probability over the folds; AlignmentQValue is the running mean of the posterior error probabilities per run
* Log an out-of-fold check of the scores and q-values
* Log the rows that are not used with their reason
* Stop with an error message on input paths that do not exist and on runs whose .quant files from different directories have the same name; match the files of a directory by position in name order
* -ts writes the training ions with truth, fold, out-of-fold score and features
* AbsentAnalyteFdr (default Estimated): with the <run>.aligndecoy files of AlignmentBasedQuantification, train with the decoys as extra negatives and add a per-run decoy based false discovery rate of absent analytes to AlignmentQValue, with the absent fraction estimated from the share of scores below 0.5 (Estimated) or set to 1 (Conservative); warn and use Off when the files are missing
* Breaking: .align files of earlier QuantBasedAlignment releases cannot be read, the runs have to be aligned again
* Add the dependency Microsoft.ML.LightGbm 1.7.1
* With -dc, save the chart AlignmentScoreHistogram with the out-of-fold AlignmentScore of the overlap ions by label, of the new ions and of the decoys; the charts ProbabilityHistogram and QValueDistribution of each run are gone
* Breaking: the parameter file has the parameters PositiveRtSdFactor (default 1), NegativeRtSdFactor (3), MaxAbsLog2QuantRatio (1), UseIntensityRatioFeatures (false), UseIonMobility (true), NumberOfFolds (5), Trainer, NumberOfLeaves, NumberOfTrees, MinimumExampleCountPerLeaf, LearningRate, AbsentAnalyteFdr and Seed (1) in place of PositiveQuantMzCutoff, NegativeQuantMzCutoff, PositiveQuantCutoff and NegativeQuantCutoff
