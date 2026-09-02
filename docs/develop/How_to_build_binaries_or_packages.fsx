(**
---
title: How to build binaries or packages
category: Developer Notes
categoryindex: 2
index: 3
---
*)
(**
# How to build binaries or packages

The build is a [FAKE](https://fake.build/) project under `build/`. `build.cmd` and `build.sh`
run it with the target name as the first argument. Without a target it builds the whole
solution.

## Building

The `Build` target compiles every released tool listed in `build/ProjectInfo.fs`:
*)

(**
```text
.\build.cmd build
```
*)

(**
The `Preprocessing` project references packages from the CSBiology feed on GitHub Packages.
`nuget.config` reads the credentials for that feed from the environment variables `GithubUser`
and `PackageToken`. Set both before building the full solution, which includes that project. The
released tools restore from nuget.org alone.

To build a single tool, call `dotnet` on its project file:
*)

(**
	dotnet build src/PeptideDB/PeptideDB.fsproj -c Release
*)

(**
## Packing

The `Pack` target creates a NuGet package for every released project in `pkg/`. The version
comes from the top entry of each project's `RELEASE_NOTES.md`, so bump that file first:
*)

(**
```text
.\build.cmd pack
```
*)

(**
`PackPrerelease` does the same with a suffix it asks for on the console. Every tool packs as a
.NET tool with the package id `ProteomIQon.<ToolName>` and the command name
`proteomiqon-<toolname>`. A package from `pkg/` can be installed locally to try it out before
publishing:
*)

(**
	dotnet tool install --global ProteomIQon.PeptideDB --add-source ./pkg --version 0.0.10
*)

(**
Published versions are on [nuget.org](https://www.nuget.org/profiles/CSBiology).

## Tests

`RunTests` builds and runs `tests/ProteomIQon.Tests.fsproj`. The end to end pipeline tests
start the tool assemblies, so this target needs the full build to succeed first:
*)

(**
```text
.\build.cmd runtests
```
*)
