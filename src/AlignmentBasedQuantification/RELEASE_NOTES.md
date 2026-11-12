#### Unreleased
* Fixed bug when using single file option caused by switched parameterss
* Read the .align files of QuantBasedAlignment with the local error scale of the predicted scan time, the source run, the source model error and the predicted ion mobility
* Choose the peak among all detected peaks near the predicted scan time with a logistic model trained per run on the ions the run quantified itself (scan time distance in units of the local error scale, shape correlation with the source trace, relative height, width ratio, number of candidates, M+1 co-elution and ratio, m/z error); a fixed score when fewer than 200 ions can be used for training; cross-fitted over 5 seeded folds of peptides (PepSequenceID) so that the ions the run quantified itself are scored out-of-fold
* Drop the dynamic time warping of the XIC onto the source trace
* Drop the filters on the scan time difference (4 sd) and the quantity to apex ratio; ions that cannot be quantified are logged with their reason
* Measure the m/z correction in the MS1 spectra in ppm over scan time and m/z, as PSMBasedQuantification does, and set the m/z window to four standard deviations of the residuals around it (was 0.05 Da); honour a fixed MzWindow_Da
* Extract the XIC over the scan time window of the PSM based quantification and take candidates within max(MinScanTimeWindow, ScanTimeWindowSdFactor * PredictedScanTimeSd), capped at ScanTimeWindow
* Extract the XICs faster and with less memory
* Ion mobility: an XIC point is the sum of all peaks inside the m/z window and the ion mobility window, as in PSMBasedQuantificationTIMs; IonMobilityWindow Estimate gives a window of at least +-0.05 around the predicted ion mobility, Fixed v sums the traces of an ion within +-v around its ion mobility centre, the intensity weighted ion mobility of the signal within v of the predicted ion mobility at the predicted scan time; ion mobility distance as candidate feature, ion mobility of the chosen peak in the IonMobility column, the ion mobility the traces are summed around in the IonMobilityCentre column
* Compare the isotope pattern on isotope traces extracted with the same windows
* Compute the ALS baseline and predict the isotope patterns as the other quantification tools
* Default to wavelet BorderPadMethod Zero, which gives identical results between runs; warn on Random
* Write the PSM based rows with a fitted peak followed by the alignment rows, one row per ion and source, and the candidate features per ion to <run>.alignfeatures; no .alignquantMetrics
* Drop an ion instead of the run when its processing fails, with the reason in the log
* Overwrite existing output files instead of appending to them
* Decoy transfers (DecoyFraction, default 0.2): a seeded fraction of the new ions is also quantified at a predicted scan time shifted by 4 to 6 min, through the same extraction, selection and fit, and written to <run>.aligndecoy for the absent-analyte false discovery rate of AlignmentBasedQuantStatistics; DecoyFraction 0 writes no decoys
* Breaking: remove the argument -iii (Metrics), the tool does not read .alignmetric files
* Breaking: remove the parameter PerformLocalWarp, add MinScanTimeWindow (default 1), ScanTimeWindowSdFactor (default 8) and Seed (default 1), and take a labeling in PerformLabeledQuantification as PSMBasedQuantification does instead of true or false, so parameter files of earlier releases have to be rewritten
* Remove the dependency on Microsoft.ML and Microsoft.ML.FastTree
* Use the Wavelet settings of the parameter file in the wavelet peak detection; they were ignored before
* Stop with an error message on input paths that do not exist and on directories with different numbers of .mzlite, .align and .quant files; match the files of a directory by position in name order, or by file name with -mf, where a run without .align or .quant file is skipped with a warning; -i, -ii and -iv have to be all files or all directories
* Breaking: .align files of earlier QuantBasedAlignment releases cannot be read, the runs have to be aligned again
* Leave the labeled partner empty in the alignment row of an ion that has a fitted PSM based row of its own
* Update MzIO to 0.1.8, MzIO.Processing to 0.1.8, MzIO.SQL to 0.1.10 and MzIO.MzML to 0.1.11
* Do not transfer an ion without a fitted PSM based row of its own whose labeled partner the run quantified through its PSM based row
* Read the isotope patterns from the table IsotopePattern of the peptide data base where PSMStatistics stored them and compute them otherwise
* Quantify ions with short XIC traces, on which the wavelet peak detection failed
* Update BioFSharp.Mz to 0.2.3
* Take the scan time of a fitted peak at the apex of its fitted shape and its width as the standard deviation of that shape, also for peaks with an exponential tail
* Add the parameter RequirePeakAtScanTime (default true): the labeled partner of a transferred ion is quantified only when a peak contains the fitted scan time of the ion, otherwise its quantity columns stay empty; false takes the peak that contains it or else the peak with the closest apex
* Breaking: the .quant files need the column IonMobilityCentre, so files of earlier PSMBasedQuantification and PSMBasedQuantificationTIMs releases have to be quantified again
* Add default parameters for timsTOF runs (AlignmentBasedQuantificationTIMsParams.json) with MzWindow_Da Fixed 0.008 and IonMobilityWindow Fixed 0.025

#### 0.0.4 - Wednesday, September 2, 2026
* Fall back to a median mz correction when the smoothing spline cannot be fitted on sparse data
* Update to .NET 10
* Update BioFSharp.Mz to 0.2.1, FSharpAux to 2.1.0, FSharp.Stats to 0.6.0 and Plotly.NET to 6.0.0-preview.2

#### 0.0.3 - Sunday, March 23, 2025
* Add Ion Mobility to output type
* Update to .NET 8

#### 0.0.2 - Monday, February 13, 2023
* Add Alignment Score fields to quantification results

#### 0.0.1 - Sunday, November 28, 2021
* Initial nuget release
