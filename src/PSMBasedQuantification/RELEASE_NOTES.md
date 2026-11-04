#### Unreleased
* Read the isotope patterns of the precursor ions from the table IsotopePattern of the peptide data base where PSMStatistics stored them and compute them otherwise; the results do not change
* Skip with a warning the spectrum descriptions that cannot be read and the MS1 spectra without a finite scan time or an ID
* Update BioFSharp.Mz to 0.2.2
* Read the spectra faster with the same results
* Compute the ALS baseline faster with the same values
* Update BioFSharp.Mz to 0.2.3
* Speed up the wavelet peak detection with the same results
* Update MzIO to 0.1.8, MzIO.Processing to 0.1.8, MzIO.SQL to 0.1.10 and MzIO.MzML to 0.1.11
* Write AlignmentQValue 0 instead of NaN for the .quant rows
* Write the diagnostic charts faster with the same content
* Write an empty .quant with a warning for a run without PSMs instead of stopping
* Use the Wavelet settings of the parameter file in the wavelet peak detection; they were ignored before
* Default to MinPeakDistance 0.1 min
* Pick the top PSMs of TopKPSMs by ModelScore instead of SequestScore

#### 0.0.12 - Wednesday, September 2, 2026
* Fall back to a median mz correction when the smoothing spline cannot be fitted on sparse data
* Write a header only output file instead of failing when no peptide ion could be quantified
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
