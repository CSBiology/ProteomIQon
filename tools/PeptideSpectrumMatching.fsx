(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/PeptideSpectrumMatching.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.ipynb)

# PeptideSpectrumMatching

PeptideSpectrumMatching identifies the peptides behind the MS/MS spectra of a run. For each MS2 scan it first determines the charge of the precursor ion from the isotope pattern in the preceding MS1 scan ([charge state determination](https://www.biofsharp.com/BioFSharp.Mz/01_03_charge_state_determination.html)).
With the charge, the precursor m/z becomes a peptide mass, and the tool fetches every peptide from the [PeptideDB](https://csbiology.github.io/ProteomIQon/tools/PeptideDB.html) database whose mass lies within LookUpPPM of it.
For each candidate it computes the theoretical fragment ions ([fragmentation](https://www.biofsharp.com/BioFSharp.Mz/02_01_fragmentation.html)) and compares them with the measured spectrum. Each candidate gets a [SEQUEST-like](https://www.biofsharp.com/BioFSharp.Mz/02_03_sequest_like_scoring.html) score and an [Andromeda-like](https://www.biofsharp.com/BioFSharp.Mz/02_04_andromeda_like_scoring.html) score. An X!Tandem-like score is computed alongside the Andromeda one.
The candidates include reversed decoy peptides, which [PSMStatistics](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html) later uses to estimate the false discovery rate ([what the decoys are for](https://www.biofsharp.com/BioFSharp.Mz/02_03_sequest_like_scoring.html#What-the-decoys-are-for)).

<img src="https://csbiology.github.io/ProteomIQon/img/PSM.png" width="1000" height="750" />
## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | one or more `.mzlite` or `.mzML` files, or a directory that is searched for `*.mzlite`, `*.mzML` and `*.mzml` (not recursive) | [MzMLToMzLite](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLite.html)
`-d` | the SQLite peptide database | [PeptideDB](https://csbiology.github.io/ProteomIQon/tools/PeptideDB.html)
`-o` | the output directory, created when missing | &#32;
`-p` | the parameter file in JSON | this page
`-c` | number of runs scored at the same time, default 1 | &#32;


All four flags are mandatory. The reader opens files by extension and knows `.mzlite` and `.mzML` only, so an mzML file has to carry the extension in exactly that spelling. The tool stops when the database file does not exist.

The tool writes one `<run>.psm` per input into the output directory. It is a tab separated table with a header. Each row is one candidate peptide for one spectrum with its charge, precursor m/z, theoretical mass, the three scores with their delta values to the next candidates, the peptide sequence and a Label of 1 for a target peptide and -1 for a decoy.
The file is opened in append mode, so running the tool twice into the same directory adds a second header and a second set of rows. Delete the old file first.
The `.psm` is read by [PSMStatistics](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html). The output directory also receives `PeptideSpectrumMatching_log.txt` and one `<run>_log.txt` per input.

## Parameters

Parameter | Default | Meaning
--- | --- | ---
ChargeStateDeterminationParams | see below | Settings of the charge state determination.
LookUpPPM | 30.0 | Mass window, in ppm around the precursor mass, from which candidate peptides are fetched.
nTerminalSeries | NTerminalSeries.B | Fragment ion series counted from the N terminus. Keep B, the scoring builds the theoretical spectra from the b and y series.
cTerminalSeries | CTerminalSeries.Y | Fragment ion series counted from the C terminus. Keep Y.
Andromeda | { PMinPMax = 4, 10; MatchingIonTolerancePPM = 100.0 } | Settings of the Andromeda-like and X!Tandem-like scoring.


ChargeStateDeterminationParams is a BioFSharp.Mz ChargeState.ChargeDetermParams record.

Parameter | Default | Meaning
--- | --- | ---
ExpectedMinimalCharge | 2 | Lowest precursor charge that is tested.
ExpectedMaximumCharge | 5 | Highest precursor charge that is tested.
Width | 1.1 | Width of the m/z window around the precursor in which the MS1 isotope pattern is examined.
MinIntensity | 0.15 | Minimum relative intensity of a peak to count as part of the isotope pattern.
DeltaMinIntensity | 0.3 | Minimum relative intensity difference used when comparing neighbouring peaks of the pattern.
NrOfRndSpectra | 10000 | Number of random spectra used to estimate how likely a matching pattern arises by chance.


Andromeda is a ProteomIQon.Domain.AndromedaParams record.

Parameter | Default | Meaning
--- | --- | ---
PMinPMax | 4, 10 | Lowest and highest number of most intense peaks kept per 100 Da window. Every count in between is tried and the best score is kept.
MatchingIonTolerancePPM | 100.0 | Tolerance in ppm for matching a theoretical fragment to a measured peak.


The default file is [peptideSpectrumMatchingParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/peptideSpectrumMatchingParams.json).

## Writing a parameter file

*)
open BioFSharp.Mz
open ProteomIQon

let chargeDetermParams : ChargeState.ChargeDetermParams =
    {
        ExpectedMinimalCharge = 2
        ExpectedMaximumCharge = 5
        Width                 = 1.1
        MinIntensity          = 0.15
        DeltaMinIntensity     = 0.3
        NrOfRndSpectra        = 10000
    }

let andromedaParams : Domain.AndromedaParams =
    {
        PMinPMax                = 4, 10
        MatchingIonTolerancePPM = 100.0
    }

let peptideSpectrumMatchingParams : Dto.PeptideSpectrumMatchingParams =
    {
        ChargeStateDeterminationParams = chargeDetermParams
        LookUpPPM                      = 30.0
        nTerminalSeries                = NTerminalSeries.B
        cTerminalSeries                = CTerminalSeries.Y
        Andromeda                      = andromedaParams
    }

// Replace the temp folder with your project folder.
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "peptideSpectrumMatchingParams.json")

Json.serializeAndWrite outputPath peptideSpectrumMatchingParams
(**
## Running the tool

Install the tool with `dotnet tool install --global ProteomIQon.PeptideSpectrumMatching`, then score one run:

```text
proteomiqon-peptidespectrummatching -i path/to/run.mzlite -d path/to/AraTest.db -o path/to/output -p path/to/peptideSpectrumMatchingParams.json
```

Several runs, three of them scored at the same time:

```text
proteomiqon-peptidespectrummatching -i path/to/run1.mzlite path/to/run2.mzlite path/to/run3.mzlite -d path/to/AraTest.db -o path/to/output -p path/to/peptideSpectrumMatchingParams.json -c 3
```

All flags:

```text
proteomiqon-peptidespectrummatching --help
```

*)