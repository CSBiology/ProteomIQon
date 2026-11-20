---
name: "RatioLFQ"
phase: protein-quantification
subgroup: Ratio LFQ
order: 4
summary: "RatioLFQ turns the peptide ion table of LabeledProteinQuantification into one light-to-heavy ratio and one intensity per protein group and run."
inputsFrom: [labeled-protein-quantification]
---

## Tool Description

A labeled experiment gives a light and a heavy intensity for every peptide ion and run. One channel is the spike-in, which is the same in every run, and the other is the sample. RatioLFQ turns the peptide ion table of [LabeledProteinQuantification](./labeled-protein-quantification.md) into one light-to-heavy ratio and one intensity per protein group and run.

## How does RatioLFQ work?

1. **Compare the runs pairwise**

   For every protein group, the tool compares each pair of runs over the ions that both runs share. The difference of two runs is the median difference of the log₂ light-to-heavy ratios of these ions. An ion specific offset of the ratio cancels in a pairwise difference, so the ratios of different runs stay comparable when different ions were quantified in different runs.

2. **Fit one ratio per run**

   The tool finds by least squares one log₂ ratio per run that fits all pairwise differences. A comparison counts with the number of ions the two runs share. Runs that share fewer than `MinSharedIons` ions with every other run, directly or through other runs, form their own group and keep their own level.

3. **Set the level of the ratios**

   The tool places the result so that the median difference between the ion ratios and the result is zero. Every ion measurement counts once in this step. A protein group seen in one run keeps the median ratio of its ions.

4. **Calculate the intensity**

   The intensity is the ratio, oriented sample over reference, times one reference intensity per protein group. That reference intensity is the mean of the up to three largest per-ion medians of the reference channel.

## Inputs

| Flag | Meaning | Comes from |
|------|---------|------------|
| `-i` | The ion table `GlobModAggregation.txt`, one row per peptide ion and run | [LabeledProteinQuantification](./labeled-protein-quantification.md) |
| `-o` | The output directory, created when missing | |
| `-p` | The parameter file in JSON | [Default parameter file](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/RatioLFQParams.json) |

The table needs the columns `FileName`, `ProteinGroup`, `StringSequence`, `Charge`, `Quant_Light` and `Quant_Heavy`. Sequence and charge together identify an ion.

The intensities have to be linear, so the LabeledProteinQuantification parameter file must not set a `LabeledTransform`. The table also has to be written with a light-to-heavy correlation filter in LabeledProteinQuantification. Without the filter, ions whose light and heavy peaks do not belong together enter the comparison of runs. [LabeledQuantificationParams_withCorrFilter.json](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/LabeledQuantificationParams_withCorrFilter.json) is an example that has one and no transform.

Rows without two positive intensities are ignored. The tool stops with a message when the table has more than one row for the same ion and run or when a required column is missing.

## Parameter file

| Parameter | Default | Description |
| --- | --- | --- |
| `Reference` | `Heavy` | `Heavy` or `Light`. The channel that is spiked in at the same amount in every run. The other channel is the sample. |
| `MinSharedIons` | `1` | The number of ions that two runs must share so that their comparison is used. |
| `CentreRuns` | `false` | Shifts the log₂ ratios of every run so that the median over all its ions is 0. This removes differences in the mixing of sample and standard and is only meaningful when most proteins do not change between runs. |

The default file is [RatioLFQParams.json](https://github.com/CSBiology/ProteomIQon/blob/main/src/ProteomIQon/defaultParams/RatioLFQParams.json):

```json
{"Reference":{"Case":"Heavy"},"MinSharedIons":1,"CentreRuns":false}
```

## Outputs

The output is `RatioLFQ.txt` in the output directory, a tab-separated table with the column `ProteinGroup` and, for every run, the columns `<run>.Ratio_LightByHeavy` and `<run>.LFQ_Intensity`. The ratio is always light over heavy, and the intensity is that of the sample. An empty cell means that the protein group has no ion in that run. The log is `RatioLFQ_log.txt` in the output directory.

## Installing and executing the tool

1. **Installing**

   ```sh
   dotnet tool install --global ProteomIQon.RatioLFQ
   ```

2. **Calculate ratios and intensities**

   ```sh
   proteomiqon-ratiolfq -i path/to/GlobModAggregation.txt -o path/to/out -p path/to/RatioLFQParams.json
   ```

3. **All flags**

   ```sh
   proteomiqon-ratiolfq --help
   ```
