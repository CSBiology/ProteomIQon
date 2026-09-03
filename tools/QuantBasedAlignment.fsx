(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/QuantBasedAlignment.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/QuantBasedAlignment.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/QuantBasedAlignment.ipynb)

# QuantBasedAlignment

In data dependent acquisition the instrument picks the ions it fragments more or less at random. A peptide identified and quantified in one run is
therefore often present in another run without ever having been fragmented there, so that run has no identification and no quantification for it.
Run alignment maps the scan time of one run onto another, so that a later tool can look for the peptide where it should elute.
QuantBasedAlignment reads the quantified peptide ions of a target run (`-i`) and of one or more source runs (`-ii`). On the peptide ions a source run
shares with the target it fits a smoothing spline from source scan times to target scan times, and with that spline it predicts the scan time of
every source peptide ion in the target run. [AlignmentBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantification.html) then extracts and quantifies the
peptide ions at the predicted scan times.

The tool links Intel MKL, so its package is platform specific: install ProteomIQon.QuantBasedAlignment_win_x64 on Windows or
ProteomIQon.QuantBasedAlignment_linux_x64 on Linux. Both install the command proteomiqon-quantbasedalignment. The package
ProteomIQon.QuantBasedAlignment on nuget.org is the library both wrap and installs no command.

## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | one or more `.quant` files of the target runs, or a directory that is searched for `*.quant` | [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html) or [PSMBasedQuantificationTIMs](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantificationTIMs.html)
`-ii` | one or more `.quant` files of the source runs, or a directory that is searched for `*.quant` | [PSMBasedQuantification](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.html) or [PSMBasedQuantificationTIMs](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantificationTIMs.html)
`-o` | the output directory, created when missing | &#32;
`-dc` | switch: write alignment fit charts | &#32;
`-c` | number of targets processed in parallel, default 1, only matters with several targets or a target directory | &#32;


There is no parameter file. Before fitting, the tool keeps only peptide ions with a successful peak fit whose fitted scan time lies within two
standard deviations of the scan time of the identification.

For every target the tool writes two tab separated files into the output directory. `<target>.align` has one row per peptide ion taken from a
source run: the peptide ion, its m/z, the predicted scan time in the target run (PredictedScanTime) and the scan time, apex intensity,
quantification and traces the ion had in the source run. The file covers every peptide ion any source run quantified, including ions the target
itself has already quantified. [AlignmentBasedQuantStatistics](https://csbiology.github.io/ProteomIQon/tools/AlignmentBasedQuantStatistics.html) uses that overlap as training
material. When several source runs contain the same peptide ion, the source with the best fit (highest R squared) supplies the row.
`<target>.alignmetric` collects the quality metrics of the fit against each source run, one block per source. AlignmentBasedQuantification reads
both files. With `-dc` the tool also writes `<target>_Metrics.html`, a chart of the fits. Delete existing `.align` and `.alignmetric` files of
the target in the output directory before a rerun, the tool appends to them. Logs go to `QuantBasedAlignment_log.txt` and one
`<target>_log.txt` per target.

With a single target the tool aligns it against all source files. With several targets, or a target directory, every target is aligned against
every source file except itself, so the same directory can be passed to `-i` and `-ii`.

## Running the tool

Align one target run against two source runs.

```text
proteomiqon-quantbasedalignment -i path/to/target.quant -ii path/to/source1.quant path/to/source2.quant -o path/to/output
```

Align every run of an experiment against all others, four at a time.

```text
proteomiqon-quantbasedalignment -i path/to/quant -ii path/to/quant -o path/to/output -c 4
```

Print the description of every argument.

```text
proteomiqon-quantbasedalignment --help
```

*)