---
name: "Add Deduced Peptides"
phase: alignment
order: 4
summary: "AddDeducedPeptides assigns the protein groups of the combined inference to every peptide present in each aligned .quant file."
inputsFrom: [alignment-based-quantification, protein-inference]
---

## Tool Description

AddDeducedPeptides takes the protein groups of the combined inference and assigns them to every peptide present in each aligned `.quant` file, so that JoinQuantPepIonsWithProteins can join the transferred quantifications to proteins.

After alignment and AlignmentBasedQuantification, a run's `.quant` file holds peptides that were never identified in that run, so the run's `.prot` file from ProteinInference does not list them.

## How does AddDeducedPeptides work?

1. **Combine protein inference results**

   The tool reads all `.prot` files as one combined inference and does not pair them with `.quant` files by name.

   This requires ProteinInference to run over all runs together (`GroupFiles = true`). Each run's `.prot` file must carry identical scores and Q-values for a protein group.

2. **Validate protein group scores**

   If a protein group has different scores, the tool stops with the message:

   > If you are running this tool please infer proteins over all files

3. **Generate run-specific protein files**

   For every input `.quant` file, the tool creates `<quant>.prot` in the output directory. For example, `run.quant` produces `run.prot`.

   Each output is a tab-separated table with the columns of a ProteinInference result.

4. **Filter peptide associations**

   Each protein group retains its scores and Q-value. The `PeptideSequence` column lists only peptides present in the corresponding `.quant` file.

   Protein groups without matching peptides are removed.

5. **Continue downstream processing**

   JoinQuantPepIonsWithProteins reads the newly generated `.prot` file together with the same `.quant` file used here. When scoring transfers, run [AlignmentBasedQuantStatistics](./alignment-based-quant-statistics.md) first and use its output in both steps.

## Inputs

These inputs are needed to execute the tool:

| Flag | Meaning | Comes from |
| --- | --- | --- |
| `-i` | One or more aligned `.quant` files, or a directory searched for `*.quant` | [AlignmentBasedQuantification](./alignment-based-quantification.md) or [AlignmentBasedQuantStatistics](./alignment-based-quant-statistics.md) |
| `-ii` | One or more `.prot` files, or a directory searched for `*.prot` | [Protein Inference ↗](./protein-inference.md) |
| `-o` | The output directory, created when missing | — |

- For AddDeducedPeptides there is no parameter file needed.
- Delete an existing output `<quant>.prot` file before a rerun, because the tool appends to it.
- Logs go to AddDeducedPeptides_log.txt.

## Outputs

Each `<quant>.prot` file contains these columns.

| Column | Description |
| --- | --- |
| ProteinGroup | Protein Identifier Collection |
| PeptideSequence | Peptide sequences assigned to this group and present in the corresponding `.quant`, separated by semicolons |
| Class | Peptide Evidence Classes (C1A, C1B, C2A, C2B, C3A, C3B) which define how a peptide has been assigned to its corresponding protein |
| TargetScore | Target Score calculated based on the individual peptide scores of PSM |
| DecoyScore | Decoy Score calculated based on the individual peptide scores of PSM |
| QValue | Q-value (False-Discovery-Rate) for every identified protein |

## Running the tool 
1. **Installing**

   Install the tool globally with the .NET CLI.

   ```sh
   dotnet tool install --global ProteomIQon.AddDeducedPeptides
   ```

2. **A single run**

   ```sh
   proteomiqon-adddeducedpeptides -i path/to/run.quant -ii path/to/run.prot -o path/to/output
   ```

3. **All runs of an experiment**

   The tool reads every `.prot` file from the input directory and writes a new `.prot` file for each `.quant` file.

   ```sh
   proteomiqon-adddeducedpeptides -i path/to/quant -ii path/to/prot -o path/to/output
   ```

4. **Help**

   ```sh
   proteomiqon-adddeducedpeptides --help
   ```
