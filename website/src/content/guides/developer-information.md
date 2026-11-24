---
title: "Developer Information"
eyebrow: "Development"
intro: "Build, test and package ProteomIQon tools, then create and document your own."
description: "Developer documentation for building, packaging, testing and documenting ProteomIQon tools."
order: 5
---

## How to build binaries or packages

The build is a FAKE project in the `./build` folder with the `build.cmd` and the `build.sh` run it with the target name as the first argument. Without a target it builds the whole solution.

## Building

The Build target compiles every released tool listed in build/ProjectInfo.fs.

\- For Windows

```sh
.\build.cmd build
```

\- For Linux:

```sh
.\build.sh build
```

The Preprocessing project references packages from the CSBiology feed on GitHub Packages. nuget.config reads the credentials for that feed from the environment variables GithubUser and PackageToken. Set both before building the full solution, which includes that project. The released tools restore from nuget.org alone. To build a single tool, call dotnet on its project file:

```sh
dotnet build src/PeptideDB/PeptideDB.fsproj -c Release
```

## Packing

The Pack target creates a NuGet package for every released project in pkg/. The version comes from the top entry of each project's RELEASE_NOTES.md that has a version, so bump that file first. A leading #### Unreleased section without a version collects changes until the next release and is skipped; replace it with the version and the date when you release: For Windows:

```sh
.\build.cmd pack
```

For Linux:

```sh
.\build.sh pack
```

PackPrerelease does the same with a suffix it asks for on the console. Every tool packs as a .NET tool with the package id ProteomIQon.`<ToolName>` and the command name proteomiqon-`<toolname>`. A package from pkg/ can be installed locally to try it out before publishing:

```sh
dotnet tool install --global ProteomIQon.PeptideDB --add-source ./pkg --version 0.0.10
```

Published versions are on [Nuget.org](https://www.nuget.org/)

## Tests

RunTests builds and runs tests/ProteomIQon.Tests.fsproj. The end to end pipeline tests start the tool assemblies, so this target needs the full build to succeed first: For Windows:

```sh
.\build.cmd runtests
```

For Linux:

```sh
.\build.sh runtests
```

## How to create your own tool

Every ProteomIQon tool is a small console project located in `src/`, containing three files with distinct responsibilities.

1. `CLIArgumentParsing.fs`

   Declares command-line arguments using Argu.

2. `Program.fs`

   Parses arguments, reads the parameter JSON and loops over input files.

3. Tool-specific implementation

   Contains the actual processing logic.

**Shared components**

Parameter records and result row types shared between tools are located in `src/ProteomIQon/` in the core project.

## Scaffolding from the template

The repository ships a `dotnet new` template that creates this layout. Follow these steps to set up your new tool.

1. **Install the template**

   Run the following command once from the repository root:

   ```text
   dotnet new install ./ConsoleTemplate/template
   ```

2. **Create your project**

   Run `dotnet new list` to verify that `pct` (the prototypical ProteomIQon console project) is available.

   Then create your project inside `src/` and add it to the solution:

   ```text
   cd src
   dotnet new pct -n MyTool --force
   cd ..
   dotnet sln ProteomIQon.sln add src/MyTool/MyTool.fsproj
   ```

**Information**

Compare the new `.fsproj` with an existing tool such as `src/PeptideDB/PeptideDB.fsproj`. Use its target framework, package versions, `PackAsTool` and `ToolCommandName` properties as references.

## Wiring the tool in

Four more places need to know about the tool:

1. A parameter record in `src/ProteomIQon/DTO.fs` with a `toDomain` conversion, plus a default JSON file in `src/ProteomIQon/defaultParams/`. The tools read parameters with `Json.ReadAndDeserialize`.
2. A `RELEASE_NOTES.md` in the project folder. The build reads the package version from it.
3. An entry in the projects list of `build/ProjectInfo.fs`, so Build and Pack include it.
4. A page in `docs/tools/`; see *How to document your work*.

**Next step**

Then open a pull request against `dev`.

## How to document your work

1. The pages on this site are literate F# scripts under `docs/`, rendered with the fsdocs tool.
2. Prose lives in `(** ... *)` comment blocks as Markdown; code between the blocks is shown with syntax highlighting and tooltips.

## Building the docs locally

1. From the repository root, restore the tools once and then run the docs target. It builds the core project in Release and renders the site into `output/`:

   ```text
   dotnet tool restore
   .\build.cmd builddocs
   ```

2. On Linux or macOS use `./build.sh builddocs`. For editing with hot reload use the watch target. It serves the site on a local port and re-renders a page whenever you save it:

   ```text
   .\build.cmd watchdocs
   ```

## Adding a page for a tool

1. Copy an existing file from `docs/tools/`, for example `PeptideDB.fsx`, rename it after the tool folder in `src/`, and adjust the frontmatter at the top. The `index` value sets the position in the sidebar. The pages are ordered along the processing chain, so pick the slot where the new tool fits and shift the later ones if needed.

2. Each tool page has the same four parts: a short description of what the tool does and where it sits in the chain, an inputs and outputs section, the parameter table with the defaults from `src/ProteomIQon/defaultParams/`, and the script that writes a parameter file followed by the command line calls. Tools without a parameter file skip the table and the script.

3. Keep the parameter script runnable. It references the core assembly built by builddocs, so this proves it before you publish:

   ```text
   dotnet fsi docs/tools/YourTool.fsx
   ```
