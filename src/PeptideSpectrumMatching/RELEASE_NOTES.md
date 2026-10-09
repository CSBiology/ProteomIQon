#### Unreleased
* About twice as fast and with a third of the memory
* Read the spectrum headers of mzlite files from the description text, one spectrum at a time, instead of deserializing every spectrum
* Find the MS1 spectrum of every MS2 spectrum by binary search
* Read the peptide data base once for all files into a mass sorted table instead of copying it into memory and querying it per spectrum
* Parse and fragment every candidate peptide once while neighbouring spectra of a charge state reach it
* Write the result rows directly instead of through the reflection based CSV writer, with floats at full precision
* Report the X!Tandem-like scores of the peptide of a row; they were paired with the Andromeda-like results by their rank, so the X!Tandem-like columns could belong to another peptide
* Number spectra with equal scan times in file order and score candidates in the order of RoundedMass and ModSequenceID, which decides the order of equal scores
* Skip spectra whose description cannot be read or lacks the ID, the scan time or the precursor m/z, with a warning
* Log failed spectra as warnings and their number per file
* Log the time spent per phase of every file
* Overwrite an existing result file instead of appending to it

#### 0.0.10 - Monday, September 28, 2026
* Write the ion mobility and hyperscore columns of the shared result record as NaN and 0
* Update BioFSharp.Mz to 0.2.2

#### 0.0.9 - Wednesday, September 2, 2026
* Update to .NET 10
* Update BioFSharp to 2.0.0, BioFSharp.Mz to 0.2.1 and FSharpAux to 2.1.0

#### 0.0.8 - Monday, March 10, 2025
* Update to .NET 8

#### 0.0.7 - Friday, July 9, 2021
* improve output formatting

#### 0.0.6 - Friday, July 9, 2021
* improve performance of charge state determination

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
