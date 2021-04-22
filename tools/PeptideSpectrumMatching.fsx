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

## Parameter Generation

Parameters are handed to the cli tool as a .json file. you can download the default file [here](https://github.com/CSBiology/ProteomIQon/blob/master/src/ProteomIQon/defaultParams/peptideSpectrumMatchingParams.json), 
or use an F# script, which can be downloaded or run in Binder at the top of the page, to write your own parameter file:

*)
#r "nuget: BioFSharp.Mz, 0.1.5-beta"
#r "nuget: Newtonsoft.Json, 12.0.3"
#r "nuget: ProteomIQon, 0.0.1"

open BioFSharp.Mz.SearchDB
open Newtonsoft.Json
open ProteomIQon
open ProteomIQon.Domain
open BioFSharp.Mz

let chargeDetermParams :ChargeState.ChargeDetermParams =
    {
        ExpectedMinimalCharge = 2
        ExpectedMaximumCharge = 5
        Width                 = 1.1
        MinIntensity          = 0.15
        DeltaMinIntensity     = 0.3
        NrOfRndSpectra        = 10000
    }

let andromedaParams: AndromedaParams =
    {
        PMinPMax                = 4,10
        MatchingIonTolerancePPM = 100.
    }

let peptideSpectrumMatchingParams :Dto.PeptideSpectrumMatchingParams =
    {
        ChargeStateDeterminationParams = chargeDetermParams 
        LookUpPPM                      = 30.
        MS2ScanRange                   = 100.,2000.
        nTerminalSeries                = NTerminalSeries.B
        cTerminalSeries                = CTerminalSeries.Y
        Andromeda                      = andromedaParams
    }

let serialized =
    peptideSpectrumMatchingParams
    |> JsonConvert.SerializeObject

System.IO.File.WriteAllText("YourPathHere",serialized)
(**
If you are running this tool in Binder, you can copy the output of the following codeblock and save it in a JSON file.

*)
serialized