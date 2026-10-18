#### Unreleased
* Read the .align files of QuantBasedAlignment with the local error scale of the predicted scan time, the source run, the source model error and the predicted ion mobility
* Choose the peak among all detected peaks near the predicted scan time with a logistic model trained per run on the ions the run quantified itself (scan time distance in units of the local error scale, shape correlation with the source trace, relative height, width ratio, number of candidates, M+1 co-elution and ratio, m/z error); a fixed score when fewer than 200 ions can be used for training; cross-fitted over 5 seeded folds of peptides (PepSequenceID) so that the ions the run quantified itself are scored out-of-fold
* Drop the dynamic time warping of the XIC onto the source trace
* Drop the filters on the scan time difference (4 sd) and the quantity to apex ratio; ions that cannot be quantified are logged with their reason
* Measure the m/z correction in the MS1 spectra and set the m/z window to four standard deviations of the measured m/z (was 0.05 Da); honour a fixed MzWindow_Da
* Extract the XIC over the scan time window of the PSM based quantification and take candidates within max(MinScanTimeWindow, ScanTimeWindowSdFactor * PredictedScanTimeSd), capped at ScanTimeWindow
* Extract the XICs faster and with less memory
* Ion mobility: window of at least +-0.05 around the predicted ion mobility, sum per m/z over the mobility scans as PSMBasedQuantificationTIMs, ion mobility distance as candidate feature, ion mobility of the chosen peak in the IonMobility column
* Compare the isotope pattern on isotope traces extracted with the same windows
* Compute the ALS baseline and predict the isotope patterns as the other quantification tools
* Default to wavelet BorderPadMethod Zero, which gives identical results between runs; warn on Random
* Write the PSM based rows with a fitted peak followed by the alignment rows, one row per ion and source, and the candidate features per ion to <run>.alignfeatures; no .alignquantMetrics
* Drop an ion instead of the run when its processing fails, with the reason in the log
* Overwrite existing output files instead of appending to them
* Decoy transfers (DecoyFraction, default 0.2): a seeded fraction of the new ions is also quantified at a predicted scan time shifted by 4 to 6 min, through the same extraction, selection and fit, and written to <run>.aligndecoy for the absent-analyte false discovery rate of AlignmentBasedQuantStatistics; DecoyFraction 0 writes no decoys
* Breaking: remove the argument -iii (Metrics), the tool does not read .alignmetric files
* Breaking: remove the parameter PerformLocalWarp and add MinScanTimeWindow (default 1), ScanTimeWindowSdFactor (default 8) and Seed (default 1), so parameter files of earlier releases have to be rewritten
* Remove the dependency on Microsoft.ML and Microsoft.ML.FastTree
* Use the Wavelet settings of the parameter file in the wavelet peak detection; they were ignored before
* Stop with an error message on input paths that do not exist and on directories with different numbers of .mzlite, .align and .quant files; match the files of a directory by position in name order, or by file name with -mf, where a run without .align or .quant file is skipped with a warning; -i, -ii and -iv have to be all files or all directories
* Breaking: .align files of earlier QuantBasedAlignment releases cannot be read, the runs have to be aligned again
* Leave the labeled partner empty in the alignment row of an ion that has a fitted PSM based row of its own
* Update MzIO to 0.1.8, MzIO.Processing to 0.1.8, MzIO.SQL to 0.1.10 and MzIO.MzML to 0.1.11
* Read the isotope patterns from the table IsotopePattern of the peptide data base where PSMStatistics stored them and compute them otherwise
