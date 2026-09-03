(**
# How to document your work

The pages on this site are literate F# scripts under `docs/`, rendered with the
[fsdocs](https://fsprojects.github.io/FSharp.Formatting/) tool. Prose lives in `(** ... *)`
comment blocks as markdown, code between the blocks is shown with syntax highlighting and
tooltips.

## Building the docs locally

From the repository root, restore the tools once and then run the docs target. It builds the
core project in Release and renders the site into `output/`:

```text
dotnet tool restore
.\build.cmd builddocs
```

On Linux or macOS use `./build.sh builddocs`. For editing with hot reload use the watch target.
It serves the site on a local port and re-renders a page whenever you save it:

```text
.\build.cmd watchdocs
```

The rendered pages link their stylesheet with the absolute site root. Open them through the
watch server. Opened from disk they load without a stylesheet.

## Adding a page for a tool

Copy an existing file from `docs/tools/`, for example `PeptideDB.fsx`, rename it after the tool
folder in `src/`, and adjust the frontmatter at the top. `index` sets the position in the
sidebar. The pages are ordered along the processing chain, so pick the slot where the new tool
fits and shift the later ones if needed.

Each tool page has the same four parts: a short description of what the tool does and where it
sits in the chain, an inputs and outputs section, the parameter table with the defaults from
`src/ProteomIQon/defaultParams/`, and the script that writes a parameter file followed by the
command line calls. Tools without a parameter file skip the table and the script. Keep the parameter script runnable. It references the core assembly built
by `builddocs`, so this proves it before you publish:

	dotnet fsi docs/tools/YourTool.fsx

BioFSharp.Mz already explains the theory, and the pages link to it. The
[BioFSharp.Mz documentation](https://www.biofsharp.com/BioFSharp.Mz/) covers signal detection,
search databases, scoring, FDR control, quantification and protein inference.

## Publishing

Pushing changes under `docs/`, `build/` or `src/ProteomIQon/` to the `dev` branch triggers the `Deploy Docs` GitHub Action.
It runs the same `builddocs` target and pushes `output/` to the `gh-pages` branch, which
GitHub Pages serves at [csbiology.github.io/ProteomIQon](https://csbiology.github.io/ProteomIQon/).
The workflow can also be started by hand from the Actions tab.

*)