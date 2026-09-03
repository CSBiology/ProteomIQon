(**
[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/ProteinInference.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/ProteinInference.ipynb)

# ProteinInference

A peptide identified by [PSMStatistics](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html) often fits more than one protein, for example two splice variants of the same gene. ProteinInference maps the identified peptides back to proteins and reports [protein groups](https://www.biofsharp.com/BioFSharp.Mz/04_02_protein_inference.html#Grouping-the-evidence-into-protein-groups): sets of proteins that share their peptide evidence.
Each group gets a score from its peptides and a decoy score from the reversed proteins, and from these a [protein level q-value](https://www.biofsharp.com/BioFSharp.Mz/04_02_protein_inference.html#Estimating-the-protein-level-FDR).

Two parameters decide how groups are formed and which peptides count for them. Take two proteins pA and pB, each found by one unique peptide and both by a third shared peptide. `Protein` (the integration strictness) decides whether pA, pB and the group `pA;pB` are all reported or whether overlapping groups are merged into the smallest set that explains every peptide. `Peptide` (the peptide usage) decides whether the abundance of `pA;pB` later comes from the shared peptide only or also from the peptides unique to pA and pB. The figure shows the outcomes for every combination.

<img src="https://csbiology.github.io/ProteomIQon/img/ProteinInference.png" width="1200" height="1000" />
## Inputs and outputs

Flag | Meaning | Comes from
--- | --- | ---
`-i` | one or more `.qpsm` files, or a directory that is searched for `*.qpsm` | [PSMStatistics](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html)
`-d` | the SQLite peptide database | [PeptideDB](https://csbiology.github.io/ProteomIQon/tools/PeptideDB.html)
`-g` | optional GFF3 annotation of the proteome | your genome annotation
`-o` | the output directory, created when missing | &#32;
`-p` | the parameter file in JSON | this page
`-dc` | switch: save diagnostic charts to the output directory | &#32;


The GFF3 file tells the tool which proteins are transcripts of the same gene. With it the `Class` column (the peptide evidence class) separates groups made of isoforms of one gene from groups that mix unrelated proteins. Without `-g` every database protein counts as its own gene and the class carries no such information.

For every input `run.qpsm` the tool writes `run.prot`, tab separated with a header and the columns `ProteinGroup` (protein identifiers joined by `;`), `PeptideSequence` (the peptides of the group found in this run, joined by `;`), `Class`, `TargetScore`, `DecoyScore` and `QValue`.
[JoinQuantPepIonsWithProteins](https://csbiology.github.io/ProteomIQon/tools/JoinQuantPepIonsWithProteins.html) and [AddDeducedPeptides](https://csbiology.github.io/ProteomIQon/tools/AddDeducedPeptides.html) read `.prot` files.

With `-dc` the tool saves a chart of scores against q-values, as `QValueGraph.html` in the output directory when the files are grouped and as `run.prot_QValueGraph.html` per run otherwise. Logs go to `ProteinInference_log.txt` and a few `ProteinInference_<step>_log.txt` files in the output directory. There is no parallelism flag.

## Parameters

Parameter | Default | Meaning
--- | --- | ---
`ProteinIdentifierRegex` | `"Cre\S+"` | Regex that extracts the protein identifier from the protein names in the database and, with `-g`, from the GFF3 entries. The default matches Chlamydomonas identifiers. Put a pattern for your own proteome here.
`Protein` | `ProteinInference.IntegrationStrictness.Maximal` | `Maximal` keeps every protein group intact. `Minimal` merges overlapping groups into the smallest set of proteins that explains all peptides.
`Peptide` | `ProteinInference.PeptideUsageForQuantification.Minimal` | `Minimal` uses only the best matching peptides of a group. `Maximal` also uses peptides that point to a larger group containing this one. `MaximalInverse` also uses peptides that point to a part of this group.
`GroupFiles` | `true` | Infer all input files together: the peptides of every run go into one inference and one q-value calculation, and every output file reports the same protein groups, each restricted to the peptides seen in that run. AddDeducedPeptides needs this. `false` infers each run on its own.
`GetQValue` | `QValueMethod.LogisticRegression FDRMethod.MAYU` | Estimate the protein FDR with MAYU and turn it into q-values by logistic regression on the target and decoy scores. The other FDR estimates are `FDRMethod.DecoyTargetRatio` and `FDRMethod.Conservative` (an FDR of 1). `QValueMethod.Storey` uses the Storey estimate instead, `QValueMethod.NoQValue` writes `nan`. If the chosen calculation fails the tool falls back to `nan` and logs it.


The default file is [ProteinInferenceParams.json](https://github.com/CSBiology/ProteomIQon/blob/dev/src/ProteomIQon/defaultParams/ProteinInferenceParams.json).

## Writing a parameter file

*)
open ProteomIQon
open ProteomIQon.Domain
open BioFSharp.Mz

let proteinInferenceParams : Dto.ProteinInferenceParams =
    {
        ProteinIdentifierRegex = @"Cre\S+"
        Protein                = ProteinInference.IntegrationStrictness.Maximal
        Peptide                = ProteinInference.PeptideUsageForQuantification.Minimal
        GroupFiles             = true
        GetQValue              = QValueMethod.LogisticRegression FDRMethod.MAYU
    }

// Replace the temp folder with your project folder.
let outputPath = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "ProteinInferenceParams.json")

Json.serializeAndWrite outputPath proteinInferenceParams
(**
## Running the tool

Install with `dotnet tool install --global ProteomIQon.ProteinInference`, then infer the proteins of one run:

```text
proteomiqon-proteininference -i path/to/run.qpsm -d path/to/database.db -o path/to/output -p path/to/ProteinInferenceParams.json
```

Several runs grouped together, with the annotation for the evidence class:

```text
proteomiqon-proteininference -i path/to/run1.qpsm path/to/run2.qpsm path/to/run3.qpsm -d path/to/database.db -g path/to/proteome.gff3 -o path/to/output -p path/to/ProteinInferenceParams.json
```

A whole directory, with the q-value chart:

```text
proteomiqon-proteininference -i path/to/qpsmFolder -d path/to/database.db -o path/to/output -p path/to/ProteinInferenceParams.json -dc
```

All flags:

```text
proteomiqon-proteininference --help
```

*)