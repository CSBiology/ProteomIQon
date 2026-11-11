#### Unreleased
* Read the isotope patterns of the precursor ions from the table IsotopePattern of the peptide data base where PSMStatistics stored them and compute them otherwise; the results do not change
* Skip with a warning the spectrum descriptions that cannot be read and the MS1 spectra without a finite scan time or an ID
* Update BioFSharp.Mz to 0.2.2
* Read the .qpsm of PSMStatistics, which MsFraggerToPSM now writes too, and fail with the missing columns when a .qpsm has another layout. Breaking: .qpsm files of earlier MsFraggerToPSM releases have to be regenerated
* Take the ion mobility of a PSM from its MS2 spectrum in the mzlite, as the .qpsm does not carry it; PSMs without a readable spectrum with ion mobility are skipped with a warning; the tool stops when none has one
* Fill MeanPercolatorScore, QValue and PEPValue from ModelScore, QValue and PEPValue of the PSMs, and pick the top PSMs of TopKPSMs by ModelScore
* Speed up the reading of the spectra and the XIC extraction with the same results
* Write the diagnostic charts faster with the same content
* Take the ion mobility and scan time of the wavelet path from the top PSMs with an ion mobility
* Take QValue, PEPValue and ProteinNames from the PSM with the lowest QValue that is not NaN
* Update BioFSharp.Mz to 0.2.3
* Compute the ALS baseline faster with the same values
* Speed up the wavelet peak detection with the same results
* Update MzIO to 0.1.8, MzIO.Processing to 0.1.8, MzIO.SQL to 0.1.10 and MzIO.MzML to 0.1.11
* Lower the memory with the same results
* Write AlignmentQValue 0 instead of NaN for the .quant rows
* Use the Wavelet settings of the parameter file in the wavelet peak detection; they were ignored before
* Default to the Wavelet peak detection with the ALS baseline correction
* Quantify ions with short XIC traces, on which the wavelet peak detection failed
* Remove the parameter FragPipe, which had no effect
* Add the parameter QualityFilterLog2Limit (default 2): the quality filter keeps the rows with -limit < log2 ((area / median area) / (apex / median apex)) < limit
* Add the parameter CentrePSMs (default 3): the number of PSMs with the highest summed matched fragment intensity whose weighted scan time and ion mobility centre the XIC
* Add IonMobilityWindow to XicExtraction: Fixed v sets the half width of the ion mobility window in 1/K0, which was 0.05; Estimate is not available in this tool
* Breaking: parameter files of earlier releases stop the tool; they need CentrePSMs, QualityFilterLog2Limit and IonMobilityWindow in XicExtraction
* Add the parameter RequirePeakAtScanTime (default true): an ion is quantified only when a peak contains the weighted scan time of its PSMs, and its labeled partner only when a peak of the partner trace contains the expected scan time; false takes the peak with the closest apex when no peak contains the scan time
* Fit the m/z correction of a run in ppm over scan time and m/z; MzWindow_Da Estimate stays four standard deviations of the precursor m/z offsets
* Sum all peaks inside the m/z window and the ion mobility window to an XIC point, and read the measured isotope pattern with the same windows
* Sum the traces of an ion within +-IonMobilityWindow around its ion mobility centre, the intensity weighted ion mobility of the signal within the window around the ion mobility of the MS/MS spectra at the weighted scan time of the identifications
* Add the column IonMobilityCentre as the last column of the .quant: the ion mobility the XIC of the row is summed around

#### 0.0.13 - Wednesday, September 2, 2026
* Add Gabor3D algorithm for peak detection
* Fall back to a median mz correction when the smoothing spline cannot be fitted on sparse data
* Make the light and heavy quality filters robust against empty results and non finite intensities
* Update to .NET 10
* Update System.Data.SQLite.Core to 1.0.119 for the memory access fix
* Update BioFSharp to 2.0.0, BioFSharp.Mz to 0.2.1, FSharpAux to 2.1.0, FSharp.Stats to 0.6.0 and Plotly.NET to 6.0.0-preview.2

#### 0.0.11 - Sunday, March 23, 2025
* Add Ion Mobility to output type

#### 0.0.10 - Monday, March 10, 2025
* Update to .NET 8

#### 0.0.9 - Monday, February 13, 2023
* Add Alignment Score fields to quantification results

#### 0.0.8 - Tuesday, April 26, 2022
* Constrain delta to be < 0.1 for mz error correction

#### 0.0.7 - Friday, August 20, 2021
* refine error logging
* improve mz refinement and added estimation for mz Window
* improve precursor accuracy calculation

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
