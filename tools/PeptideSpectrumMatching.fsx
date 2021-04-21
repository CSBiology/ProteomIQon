(**
// can't yet format YamlFrontmatter (["title: PeptideSpectrumMatching"; "category: Tools"; "categoryindex: 1"; "index: 3"], Some { StartLine = 2 StartColumn = 0 EndLine = 6 EndColumn = 8 }) to pynb markdown

# Peptide Spectrum Matching

This tool selects potential peptides for each measured spectrum and computes theoretical spectra for them. Those theoretical spectra are then compared to the measured specctrum and scored according to their similarity. 
This way, the most likely peptide from which the measured spectrum originated from can be identified using the similarity scores.

## Parameters

| **Parameter**                  | **Default Value**                                                                                                                         | **Description**                                                    |
|--------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------|
| ChargeStateDeterminationParams | {ExpectedMinimalCharge = 2; ExpectedMaximumCharge = 5; Width = 1.1; MinIntensity = 0.15; DeltaMinIntensity = 0.3; NrOfRndSpectra = 10000} | Parameters used for the charge state determination of the peptides |
| LookUpPPM                      | 30                                                                                                                                        | Mass range in Da in which potential peptides are selected          |
| MS2ScanRange                   | 100.,2000.                                                                                                                                | m/z range for MS2 spectra                                          |
| nTerminalSeries                | NTerminalSeries.B                                                                                                                         | Considered ions starting from the N-Terminus                       |
| cTerminalSeries                | CTerminalSeries.Y                                                                                                                         | Considered ions starting from the C-Terminus                       |
| Andromeda                      | {PMinPMax = 4,10; MatchingIonTolerancePPM = 100.}                                                                                         | Andromeda scoring parameters                                       |

*)