(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/RatioLFQ.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/RatioLFQ.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/RatioLFQ.ipynb)

# RatioLFQ

A labeled experiment gives one light to heavy ratio per protein and run. Ratios of different runs are only comparable through their shared reference, so RatioLFQ turns them into one intensity per run and protein, the way label free quantification would report it. For every protein it takes the log2 ratio of each run, forms the difference of the log2 ratios for every pair of runs, and solves a least squares problem for one log2 value per run that reproduces these differences. The system is underdetermined by one degree of freedom, and the tool fixes it with the constraint that the sum of the run values equals the mean log2 reference intensity of the protein times the number of runs with a ratio. The solution is exponentiated back and written as the LFQ value of each run.

## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | one tab separated protein table, typically LabeledQuant.txt | [LabeledProteinQuantification](https://csbiology.github.io/ProteomIQon/tools/LabeledProteinQuantification.html)
`-o` | the path of the output file (not a directory), its folder has to exist | &#32;
`-pc` | name of the column that identifies the protein, ProteinGroup for that file | &#32;
`-rac` | suffix that selects the ratio columns, `.Ratio_LightByHeavy` for that file | &#32;
`-rfc` | suffix that selects the reference columns, `.Quant_Heavy` for that file | &#32;


`-rac` and `-rfc` are plain suffix tests on the column name, so in LabeledQuant.txt, where the columns are named `<run>.Ratio_LightByHeavy` and `<run>.Quant_Heavy`, pass the suffixes with the leading dot. Without the dot, `Quant_Heavy` also matches CV_Quant_Heavy, StDev_Quant_Heavy and SEM_Quant_Heavy, and those would go into the reference mean.

The output file is a tab separated table with a ProteinGroup column followed by one column per ratio column, named like the input column with the suffix `_LFQ`. Proteins with a ratio in only one run are left out because there is no pair of runs to compare. Runs without a ratio for a protein get a 0 in that cell, which stands for no ratio, not for a zero intensity.

## Running the tool

The tool installs with `dotnet tool install --global ProteomIQon.RatioLFQ`. A call on the labeled protein table:

```text
proteomiqon-ratiolfq -i path/to/LabeledQuant.txt -o path/to/LabeledQuant_LFQ.txt -pc ProteinGroup -rac .Ratio_LightByHeavy -rfc .Quant_Heavy
```

All flags:

```text
proteomiqon-ratiolfq --help
```

*)