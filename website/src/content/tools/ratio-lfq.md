---
name: "Ratio LFQ"
phase: protein-quantification
subgroup: Ratio LFQ
order: 4
summary: "RatioLFQ turns the light to heavy ratios of different runs into one intensity per run and protein, the way label free quantification would report it."
inputsFrom: []
---

## Tool Description

A labeled experiment gives one light to heavy ratio per protein and run. Ratios of different runs are only comparable through their shared reference, so RatioLFQ turns them into one intensity per run and protein, the way label free quantification would report it. of each run.

## What's happening?

For every protein it takes the log2 ratio of each run, forms the difference of the log2 ratios for every pair of runs, and solves a least squares problem for one log2 value per run that reproduces these differences. The system is underdetermined by one degree of freedom, and the tool fixes it with the constraint that the sum of the run values equals the mean log2 reference intensity of the protein times the number of runs with a ratio. The solution is exponentiated back and written as the LFQ value
