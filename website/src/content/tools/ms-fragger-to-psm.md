---
name: "MSFragger to PSM"
phase: peptide-identification
subgroup: MSFragger import
order: 4
summary: "MsFraggerToPSM turns the psm.tsv result of an MSFragger search, as written by FragPipe, into the .qpsm files that [PSMBasedQuantificationTIMs](./psm-based-quantification-tims.md) reads."
inputsFrom: [mzml-to-mzlite-ion-mobility, peptide-db]
---

## Tool Description

MsFraggerToPSM turns the `psm.tsv` result of an MSFragger search, as written by FragPipe, into the `.qpsm` files that [PSMBasedQuantificationTIMs](./psm-based-quantification-tims.md) reads.

## How does MsFraggerToPSM work?

1. **The advantage**

   With this tool, an ion mobility run can be identified in FragPipe and quantified in ProteomIQon, in place of the [PeptideSpectrumMatching](./peptide-spectrum-matching.md) and [PSMStatistics](./psm-statistics.md) pair.

2. **How does this work?**

   For each row of the tsv the tool rebuilds the modified peptide sequence in ProteomIQon notation, looks it up in the [PeptideDB](./peptide-db.md) database to get the peptide and modified sequence identifiers, and looks for the corresponding MS2 spectrum in the `.mzlite`. Candidates must have the same precursor charge as the TSV row or have no recorded charge. Among candidates within 0.001 of the reported ion mobility, it minimizes the sum of the absolute scan time and precursor m/z differences. If no candidate meets the mobility criterion or the TSV omits ion mobility, it uses the closest charge-compatible candidate. A row without a charge-compatible candidate is dropped.

## Inputs

These inputs are needed to execute the tool:

### Input flags

| Flag | Meaning | Comes from |
| --- | --- | --- |
| `-i` | one or more `psm.tsv` files, or a directory that is searched for `*.tsv` | MSFragger (FragPipe) |
| `-ii` | the matching `.mzlite` files, or a directory searched for `*.mzlite` | [MzMLToMzLiteIonMobility](./mzml-to-mzlite-ion-mobility.md) |
| `-o` | the output directory, created when missing | - |
| `-d` | the SQLite peptide database | [PeptideDB](./peptide-db.md) |
| `-c` | number of pairs converted at the same time, default 1 | - |

`-i`, `-ii`, `-o` and `-d` are required. This tool takes no parameter file.

### Required TSV columns

The tsv needs the columns:

- `Assigned Modifications`
- `Peptide`
- `Retention`
- `Charge`
- `Calibrated Observed M/Z`
- `Calculated Peptide Mass`
- `Delta Mass`
- `Peptide Length`
- `Number of Missed Cleavages`
- `Hyperscore`
- `Probability`
- `Protein`

`Ion Mobility` is optional. Retention is read in seconds and stored in minutes.

### File pairing and database requirements

The files are paired by position, the first tsv with the first mzlite and so on. Give both lists in the same order and with the same length, otherwise the tool converts nothing. The tool stops when the database file does not exist. The lookup uses the unlabeled sequences of the database. Rows whose modified sequence has no entry in the database are dropped, so the database has to be built with the same modifications as the FragPipe search.

### Supported modifications

From the Assigned Modifications column the tool understands oxidized methionine and carbamidomethylated cysteine, plus N-terminal acetylation. A row with any other modification is printed to the console and dropped.

### Generated files and logs

The tool writes one `run.qpsm` per mzlite into the output directory, named after the mzlite. It uses the PSMStatistics column layout. `PSMId` identifies the matched `.mzlite` spectrum and `ScanNr` gives its position in scan-time order among the MS2 spectra. The import sets `Label` to `1` and `GlobalMod` to `0`. `ModelScore` is the FragPipe probability and `PEPValue` is one minus that probability. `QValue` is the mean PEP across TSV rows with a PEP at or below that of the current row. `XtandemScore` holds the MSFragger hyperscore. The other search scores are `NaN`. [PSMBasedQuantificationTIMs](./psm-based-quantification-tims.md) retrieves ion mobility from the matched spectrum.

An existing `.qpsm` of the same name is appended to. Delete the old file before a rerun. The output directory also receives `MSFraggerToPSM_log.txt`.

## Installing and executing the tool

### Installing

Install the tool globally with the .NET CLI.

```sh
dotnet tool install --global ProteomIQon.MsFraggerToPSM
```

### Convert one run

```sh
proteomiqon-msfraggertopsm -i path/to/run/psm.tsv -ii path/to/run.mzlite -d path/to/AraTest.db -o path/to/output
```

### Convert multiple runs

```sh
proteomiqon-msfraggertopsm -i path/to/run1/psm.tsv path/to/run2/psm.tsv path/to/run3/psm.tsv -ii path/to/run1.mzlite path/to/run2.mzlite path/to/run3.mzlite -d path/to/AraTest.db -o path/to/output -c 3
```

### All flags

```sh
proteomiqon-msfraggertopsm --help
```
