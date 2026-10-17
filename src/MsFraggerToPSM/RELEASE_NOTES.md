#### Unreleased
* Update BioFSharp.Mz to 0.2.2
* Write the .qpsm in the column layout of PSMStatistics: ModelScore is the FragPipe probability, PEPValue is 1 - probability, QValue is the mean PEP of all PSMs with an equal or lower PEP, XtandemScore is the hyperscore, the other scores are NaN. Breaking: the ion mobility, expectation value and hyperscore columns are gone
* Assign a PSM only to spectra of its charge or without a charge state and prefer the spectra with its ion mobility (within 0.001) before scan time and precursor m/z. Rows without a spectrum of their charge are dropped
* Exclude a PSM with a NaN PEP from the q-values of the others and give it the q-value NaN. An empty Probability cell is not supported
* Number ScanNr with equal scan times in file order and without unreadable spectra, as PeptideSpectrumMatchingTIMs does
* Update BioFSharp.Mz to 0.2.3
* Update NLog to 4.7.11

#### 0.0.4 - Wednesday, September 2, 2026
* Update to .NET 10
* Update System.Data.SQLite.Core to 1.0.119 for the memory access fix
* Update BioFSharp.Mz to 0.2.1 and FSharp.Stats to 0.6.0

#### 0.0.3 - Sunday, March 23, 2025
* Fix input file pattern matching

#### 0.0.2 - Sunday, March 23, 2025
* Fix input file assignment

#### 0.0.1 - Wednesday, March 19, 2025
* Initial release
