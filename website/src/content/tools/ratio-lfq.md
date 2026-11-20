---
name: "Ratio LFQ"
phase: protein-quantification
subgroup: Ratio LFQ
order: 4
summary: "RatioLFQ turns the light to heavy ratios of different runs into one intensity per run and protein, the way label free quantification would report it."
inputsFrom: [labeled-protein-quantification]
---

## Tool Description

A labeled experiment gives one light-to-heavy ratio per protein and run. RatioLFQ uses the shared heavy reference to convert these ratios into one comparable intensity per protein and run.

## How does RatioLFQ work?

1. **Fit the differences between runs**

   For each protein, the tool takes the log₂ of the available ratios and calculates the difference for every pair of runs. It fits one log₂ intensity per run by least squares.

2. **Set the intensity scale**

   The differences leave one free parameter. The tool constrains the mean fitted log₂ intensity to equal the mean log₂ reference intensity for that protein. It exponentiates each fitted value with base 2 to obtain the LFQ intensity.

## Inputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | One tab-separated protein table, typically `LabeledQuant.txt` | [LabeledProteinQuantification](./labeled-protein-quantification.md) |
| `-o` | Path to the output file. Its parent directory must exist. | |
| `-pc` | Protein identifier column, typically `ProteinGroup` | |
| `-rac` | suffix that selects the ratio columns, `.Ratio_LightByHeavy` for that file | |
| `-rfc` | suffix that selects the reference columns, `.Quant_Heavy` for that file | |

`-rac` and `-rfc` match column names by suffix. For `LabeledQuant.txt`, use the leading dot in both suffixes to match `<run>.Ratio_LightByHeavy` and `<run>.Quant_Heavy`. A reference suffix of `Quant_Heavy` also matches columns ending in `CV_Quant_Heavy`, `StDev_Quant_Heavy` and `SEM_Quant_Heavy` and includes them in the reference mean.

## Outputs

The output is a tab-separated table with a `ProteinGroup` column and one intensity column per input ratio column. Each intensity column keeps the name of its input column and adds the suffix `_LFQ`. A protein needs ratios from at least two runs to appear in the table. A value of `0` marks a missing ratio for that run.

## Installing and executing the tool

1. **Installing**

   ```sh
   dotnet tool install --global ProteomIQon.RatioLFQ
   ```

2. **Calculate LFQ intensities**

   ```sh
   proteomiqon-ratiolfq -i path/to/LabeledQuant.txt -o path/to/LabeledQuant_LFQ.txt -pc ProteinGroup -rac .Ratio_LightByHeavy -rfc .Quant_Heavy
   ```

3. **All flags**

   ```sh
   proteomiqon-ratiolfq --help
   ```
