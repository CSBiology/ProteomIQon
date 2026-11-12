(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/MsFraggerToPSM.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/MsFraggerToPSM.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/MsFraggerToPSM.ipynb)

# MsFraggerToPSM

MsFraggerToPSM turns the `psm.tsv` result of an MSFragger search, as written by FragPipe, into the `.qpsm` files that [PSMBasedQuantificationTIMs](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantificationTIMs.html) reads.
With it, an ion mobility run can be identified in FragPipe and quantified in ProteomIQon, in place of the [PeptideSpectrumMatching](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.html) and [PSMStatistics](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html) pair.
For each row of the tsv the tool rebuilds the modified peptide sequence in ProteomIQon notation, looks it up in the [PeptideDB](https://csbiology.github.io/ProteomIQon/tools/PeptideDB.html) database to get the peptide and modified sequence identifiers, and assigns the row to an MS2 spectrum of the mzlite. Among the spectra of the same charge (or without charge) it takes the one whose scan time and precursor m/z lie closest to the values in the tsv, and it looks first among those whose ion mobility is within 0.001 of the tsv value. A row without any spectrum of its charge is dropped.
The tool has no parameter file.

## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | one or more `psm.tsv` files, or a directory that is searched for `*.tsv` | MSFragger (FragPipe)
`-ii` | the matching `.mzlite` files, or a directory that is searched for `*.mzlite` | [MzMLToMzLiteIonMobility](https://csbiology.github.io/ProteomIQon/tools/MzMLToMzLiteIonMobility.html)
`-o` | the output directory, created when missing | &#32;
`-d` | the SQLite peptide database | [PeptideDB](https://csbiology.github.io/ProteomIQon/tools/PeptideDB.html)
`-c` | number of pairs converted at the same time, default 1 | &#32;


All four flags are mandatory. The tsv needs the columns Assigned Modifications, Peptide, Retention, Charge, Calibrated Observed M/Z, Calculated Peptide Mass, Delta Mass, Peptide Length, Number of Missed Cleavages, Hyperscore, Probability and Protein, and should have Ion Mobility. Retention is read in seconds and stored in minutes.
Files are paired by position, the first tsv with the first mzlite and so on. Give both lists in the same order and with the same length, otherwise the tool converts nothing.
The tool stops when the database file does not exist. The lookup uses the unlabeled sequences of the database. Rows whose modified sequence has no entry in the database are dropped, so the database has to be built with the same modifications as the FragPipe search.

From the Assigned Modifications column the tool understands oxidized methionine and carbamidomethylated cysteine, plus N-terminal acetylation. A row with any other modification is printed to the console and dropped.

The tool writes one `<run>.qpsm` per mzlite into the output directory, named after the mzlite. It has the column layout of the `.qpsm` of [PSMStatistics](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html), with one row per identified spectrum and the spectrum identifier of the mzlite as `PSMId`. `ScanNr` is the position of that spectrum among the MS2 spectra in scan time order, `Label` is 1 and `GlobalMod` is 0. `ModelScore` is the FragPipe probability, `PEPValue` is 1 minus the probability and `QValue` is the mean `PEPValue` of all rows of the tsv with an equal or lower `PEPValue`. `XtandemScore` is the MSFragger hyperscore, and the other search scores are `NaN`. PSMBasedQuantificationTIMs takes the ion mobility from the spectrum.
An existing `.qpsm` of the same name is appended to. Delete the old file before a rerun.
The output directory also receives `MSFraggerToPSM_log.txt`.

## Running the tool

Install the tool with `dotnet tool install --global ProteomIQon.MsFraggerToPSM`, then convert one run:

```text
proteomiqon-msfraggertopsm -i path/to/run/psm.tsv -ii path/to/run.mzlite -d path/to/AraTest.db -o path/to/output
```

Several runs, paired by position and converted three at a time:

```text
proteomiqon-msfraggertopsm -i path/to/run1/psm.tsv path/to/run2/psm.tsv path/to/run3/psm.tsv -ii path/to/run1.mzlite path/to/run2.mzlite path/to/run3.mzlite -d path/to/AraTest.db -o path/to/output -c 3
```

All flags:

```text
proteomiqon-msfraggertopsm --help
```

*)