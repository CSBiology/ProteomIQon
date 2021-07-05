(**
---
title: How to build binaries or packages
category: Developer Notes
categoryindex: 2
index: 3
---
*)
(**
# How to build the binaries

If you want to build the binaries for all tools at one, you can simply navigate to the project root in the console and call:

	.\build.cmd

If you want to build only a singe project, then you can use the "BuildProj" target like this:

	.\build.cmd -t BuildProj

After the "Clean" target you will then be prompted to enter the name of the project you wish to build.

# How to pack a nuget package

If you want to pack a nuget package or build a dotnet tool, you can do this for all tools at once by invoking the "Pack" or "PackPrerelease" target.

	.\build.cmd -t Pack

Packing only a single project works similar to the same case for the binaries. The targets for this are calles "PackProj" and "PackPrereleaseProj":

	.\build.cmd -t PackProj

After the "Clean" target you will then be prompted to enter the name of the project you wish to build.
*)
