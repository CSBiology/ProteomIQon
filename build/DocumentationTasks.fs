module DocumentationTasks

open Helpers
open ProjectInfo
open BasicTasks

open BlackFox.Fake
open Fake.DotNet

/// The docs pages only reference the core assembly, so the docs targets build just that
/// project instead of every released tool. That keeps the docs workflow on GitHub Actions
/// short and independent of the tools.
let buildCore =
    BuildTask.create "BuildCore" [] {
        CoreProject.ProjFile
        |> DotNet.build (fun p ->
            { p with
                Configuration = DotNet.BuildConfiguration.Release
                MSBuildParams = { p.MSBuildParams with DisableInternalBinLog = true }
            }
        )
    }

let buildDocs =
    BuildTask.create "BuildDocs" [ buildCore ] {
        printfn "building docs with stable version %s" stableDocsVersionTag

        runDotNet
            (sprintf
                "fsdocs build --clean --noapidocs --properties Configuration=Release --parameters fsdocs-package-version %s"
                stableDocsVersionTag)
            "./"
    }

let buildDocsPrerelease =
    BuildTask.create "BuildDocsPrerelease" [ setPrereleaseTag; buildCore ] {
        printfn "building docs with prerelease version %s" prereleaseTag

        runDotNet
            (sprintf
                "fsdocs build --clean --noapidocs --properties Configuration=Release --parameters fsdocs-package-version %s"
                prereleaseTag)
            "./"
    }

let watchDocs =
    BuildTask.create "WatchDocs" [ buildCore ] {
        printfn "watching docs with stable version %s" stableDocsVersionTag

        runDotNet
            (sprintf
                "fsdocs watch --clean --noapidocs --properties Configuration=Release --parameters fsdocs-package-version %s"
                stableDocsVersionTag)
            "./"
    }

let watchDocsPrerelease =
    BuildTask.create "WatchDocsPrerelease" [ setPrereleaseTag; buildCore ] {
        printfn "watching docs with prerelease version %s" prereleaseTag

        runDotNet
            (sprintf
                "fsdocs watch --clean --noapidocs --properties Configuration=Release --parameters fsdocs-package-version %s"
                prereleaseTag)
            "./"
    }
