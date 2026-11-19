---
title: "Installation Guide"
eyebrow: "Set up the toolchain"
intro: "Install the ProteomIQon library or choose individual tools."
description: "Install the ProteomIQon library or individual command-line tools."
order: 2
---

## Installation of ProteomIQon as a Library

The [ProteomIQon package on NuGet](https://www.nuget.org/packages/ProteomIQon/) installs only the library. To load it in an F# script or F# Interactive, use the following directive.

```fsharp
#r "nuget: ProteomIQon, 0.1.0"
```

This approach is for Linux and Windows suitable.

## Installation of Individual Tools

To install individual tools of ProteomIQon you can use for example, Nuget.org [PeptideSpectrumMatching](https://www.nuget.org/packages/ProteomIQon.PeptideSpectrumMatching). Alternatively, you can install individual tools from the command line using the following command:

```sh
dotnet tool install --global ProteomIQon.[ToolName] --version [x.x.x]
```

This approach is for Linux and Windows suitable.

## Important to know

Most of the ProteomIQon tools need a param file (in .JSON format) to run. The param file contains all the parameters for the tool. Detailed information about the single parameters can be found in the respective tool documentation. [ProteomIQon Tools](../tools/)
