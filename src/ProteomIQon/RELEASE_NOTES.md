#### Unreleased
* Add StoreIsotopePatterns to the PSMStatistics parameters and the default parameters
* Add IsotopePatterns: the isotope pattern of a peptide at a charge (MIDA, resolution 0.01, minimum probability 0.001) for the quantification tools, its BLOB encoding, the table IsotopePattern of the peptide data base with a look up that reads the table and computes missing patterns, and the function that stores the patterns
* Add SpectrumIndex (the MS1 index and the MS2 headers from one pass over the spectrum descriptions, skipping with a warning the descriptions that cannot be read) and SortedSearch (binary searches and the closest scan time), and the m/z sort of the peaks of a spectrum in SortedPeaks, shared by the quantification and search tools
* Use CrossValidation and no PEP cutoff (PepValueThreshold 1.0) in the default PSMStatistics parameters
* Add Charts.useReferenceEqualityForJson, which makes saving Plotly.NET charts faster with the same output
* Add the local scan time error scale, the source run, the source model error and the predicted ion mobility with its error scale to AlignmentResult, the .align row of QuantBasedAlignment
* Breaking: change the parameters of AlignmentBasedQuantification (AlignmentBasedQuantificationParams) and its default parameters: PerformLocalWarp is removed, PerformLabeledQuantification is a Labeling as in QuantificationParams, and MinScanTimeWindow, ScanTimeWindowSdFactor, DecoyFraction, Seed and IonMobilityWindow in XicExtraction are added
* Breaking: replace the parameters of AlignmentBasedQuantStatistics (AlignmentBasedQuantStatisticsParams), PositiveQuantMzCutoff, NegativeQuantMzCutoff, PositiveQuantCutoff and NegativeQuantCutoff, by the label, feature, trainer, absent-analyte FDR and seed parameters, and change its default parameters
* Remove PSMStatisticsResultFragpipe. MsFraggerToPSM writes and PSMBasedQuantificationTIMs reads PSMStatisticsResult
* Write warnings to the log file of the tools
* Add CrossValidation to the PSMStatistics EstimationParams
* Add BaselineAls, the ALS baseline of FSharp.Stats baselineAls' computed faster with the same values, shared by PSMBasedQuantification, PSMBasedQuantificationTIMs and AlignmentBasedQuantification
* Update BioFSharp.Mz to 0.2.3
* Speed up the wavelet peak detection of FSharpStats'.Wavelet with the same results
* Update MzIO to 0.1.8, MzIO.Processing to 0.1.8, MzIO.SQL to 0.1.10 and MzIO.MzML to 0.1.11
* FSharpStats'.Wavelet.identify uses the parameters it is given, which were ignored before, and handles short traces, on which it failed
* Set MinPeakDistance 0.1 in the default parameters QuantificationParams
* Use the Wavelet peak detection with the ALS baseline correction in the default parameters QuantificationTIMsParams
* Add RetentionTimeFeature to the PSMStatistics EstimationParams, true in the default parameters
* Correct the default parameters QuantificationParams_NoBaseline, which could not be read, and make them equal to QuantificationParams except for the baseline correction
* Set Alignment_QValue 0.05 in the default parameters of LabelFreeProteinQuantification and LabeledProteinQuantification
* Set MinPepLength 6 in the default parameters of PeptideDB
* Add Xic, MzCalibration, PeakQuantification and PSMQuantification: the XIC of an ion (on spectra with ion mobility the sum of all peaks inside the m/z window and the ion mobility window, for a fixed window around the ion mobility centre of the ion), the m/z correction of a run, the peak detection and the choice of the peak, and the quantification of the PSMs of an ion with its quality filter, shared by PSMBasedQuantification, PSMBasedQuantificationTIMs and AlignmentBasedQuantification
* Breaking: change QuantificationParams, the parameters of PSMBasedQuantification and PSMBasedQuantificationTIMs, and the default parameters: remove FragPipe and add QualityFilterLog2Limit, CentrePSMs, RequirePeakAtScanTime and IonMobilityWindow in XicExtraction
* Remove DTW', the dynamic time warping of the alignment tools
* Add PeakShape, the apex and the standard deviation of a fitted peak from its Gaussian or exponentially modified Gaussian parameters, used by the alignment tools
* Add the column IonMobilityCentre as the last column of the .quant row (QuantificationResult): the ion mobility the XIC of the row is summed around, NaN without an ion mobility window. Breaking: QuantBasedAlignment, AlignmentBasedQuantification and AlignmentBasedQuantStatistics do not read .quant files without it
* Set MzWindow_Da Fixed 0.008 and IonMobilityWindow Fixed 0.025 in the default parameters QuantificationTIMsParams and add AlignmentBasedQuantificationTIMsParams, the default parameters of AlignmentBasedQuantification with the same two values

#### 0.0.12 - Monday, September 28, 2026
* Extend PeptideSpectrumMatchingResult with IonMobility, Hyperscore, Expectscore, MatchedIons and TotalIons. Breaking: .psm files written by earlier releases have 23 columns and have to be regenerated before PSMStatistics of this release reads them
* Add the PeptideSpectrumMatchingTIMs parameter record and default parameters, with MinFragmentsModelling 1
* Keep one spline knot per score in the PEP estimation, so tied scores no longer get a NaN PEP
* Make the PEPs monotone by pooling adjacent violators, so a few high scoring decoys no longer set almost every PEP to 1
* Update BioFSharp.Mz to 0.2.2

#### 0.0.11 - Wednesday, September 2, 2026
* Add Gabor3D parameters and default params for PSMBasedQuantificationTIMs
* Write a header only csv instead of failing when a result sequence is empty
* Remove modules that moved upstream into BioFSharp.Mz (SparsePeakArray, Fitting, SearchDB, ProteinInference, FDRControl, Fragmentation)
* Update to .NET 10
* Update System.Data.SQLite.Core to 1.0.119 for the memory access fix
* Update BioFSharp to 2.0.0, BioFSharp.Mz to 0.2.1, FSharpAux to 2.1.0, FSharp.Stats to 0.6.0, Plotly.NET to 6.0.0-preview.2 and DynamicObj to 7.0.1

#### 0.0.10 - Monday, June 02, 2025
* Fix transform for edge cases in pep value calculation (by @paulineHans)

#### 0.0.9 - Sunday, March 23, 2025
* Upgrade to .NET 8
* Add Ion Mobility to relevant output types

#### 0.0.8 - Monday, February 13, 2023
* Add Alignment Score fields to quantification results

#### 0.0.7 - Monday, September 13, 2021
* initial nuget release of
	* LabelFreeProteinQuantification
	* LabeledProteinQuantification
	* JoinQuantPepIonsWithProteins

#### 0.0.6 - Friday, August 20, 2021
* update PSMStats Params
* update PSMBasedQuantification Params
* minor changes

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
