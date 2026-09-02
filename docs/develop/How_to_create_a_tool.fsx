(**
---
title: How to create your own tool
category: Developer Notes
categoryindex: 2
index: 1
---
*)
(**
# How to create a tool

Every ProteomIQon tool is a small console project under `src/` with the same three files:
`CLIArgumentParsing.fs` declares the command line arguments with
[Argu](https://fsprojects.github.io/Argu/), `Program.fs` parses them, reads the parameter JSON
and loops over the input files, and a third file holds the actual work. Parameter records and
result row types shared between tools live in the core project `src/ProteomIQon/`.

## Scaffolding from the template

The repository ships a `dotnet new` template that creates this layout. Install it once from
the repository root:
*)

(**
	dotnet new install ./ConsoleTemplate/template
*)

(**
`dotnet new list` should now show `pct` (prototypical proteomiqon console project). Create the
project inside `src/` and add it to the solution:
*)

(**
	cd src
	dotnet new pct -n MyTool --force
	cd ..
	dotnet sln ProteomIQon.sln add src/MyTool/MyTool.fsproj
*)

(**
Compare the new `.fsproj` with a current tool such as `src/PeptideDB/PeptideDB.fsproj` and take
the target framework, the package versions, the `PackAsTool` and `ToolCommandName` properties
from there.

## Wiring the tool in

Four more places need to know about the tool:

1. A parameter record in `src/ProteomIQon/DTO.fs` with a `toDomain` conversion, plus a default
   JSON file in `src/ProteomIQon/defaultParams/`. The tools read parameters with
   `Json.ReadAndDeserialize`.
2. A `RELEASE_NOTES.md` in the project folder. The build reads the package version from it.
3. An entry in the `projects` list of `build/ProjectInfo.fs`, so `Build` and `Pack` include it.
4. A page in `docs/tools/`, see [How to document your work]({{root}}develop/How_to_document_your_work.html).

Then open a pull request against `dev`.
*)
