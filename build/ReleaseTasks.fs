module ReleaseTasks

open MessagePrompts
open ProjectInfo
open BasicTasks
open TestTasks
open PackageTasks
open DocumentationTasks

open BlackFox.Fake
open Fake.Core
open Fake.DotNet
open Fake.Api
open Fake.Tools
open Fake.IO
open Fake.IO.Globbing.Operators


let createTag =
    BuildTask.create "CreateTag" [ clean; build; runTests; pack; buildDocs.IfNeeded ] {
        if promptYesNo (sprintf "tagging branch with %s OK?" branchTag) then
            Git.Branches.tag "" branchTag
            Git.Branches.pushTag "" projectRepo branchTag
        else
            failwith "aborted"
    }

let createPrereleaseTag =
    BuildTask.create
        "CreatePrereleaseTag"
        [
            setPrereleaseTag
            clean
            build
            runTests
            packPrerelease
            buildDocsPrerelease.IfNeeded
        ] {
        if promptYesNo (sprintf "tagging branch with %s OK?" prereleaseTag) then
            Git.Branches.tag "" prereleaseTag
            Git.Branches.pushTag "" projectRepo prereleaseTag
        else
            failwith "aborted"
    }


let publishNuget =
    BuildTask.create "PublishNuget" [ clean; build; runTests; pack; createTag.IfNeeded ] {
        let targets =
            (!!(sprintf "%s/*.*pkg" pkgDir))

        printfn "package files:" 

        for target in targets do
            printfn "%A" target

        printfn "package versions to release:" 

        projects
        |> List.iter (fun p ->  
            printfn $"{p.Name} @ {p.PackageVersionTag}"
        )

        if promptYesNo "OK?" then
            let source =
                "https://api.nuget.org/v3/index.json"

            let apikey =
                Environment.environVar "NUGET_KEY"

            if System.String.IsNullOrWhiteSpace apikey then
                failwith "NUGET_KEY environment variable is not set - cannot publish"

            for artifact in targets do
                let result =
                    DotNet.exec id "nuget" (sprintf "push -s %s -k %s %s --skip-duplicate" source apikey artifact)

                if not result.OK then
                    failwith "failed to push packages"
        else
            failwith "aborted"
    }

let publishNugetPrerelease =
    BuildTask.create
        "PublishNugetPrerelease"
        [
            clean
            build
            runTests
            packPrerelease
            createPrereleaseTag.IfNeeded
        ] {
        let targets =
            (!!(sprintf "%s/*.*pkg" pkgDir))

        printfn "package files:" 

        for target in targets do
            printfn "%A" target

        printfn "package versions to release:" 

        projects
        |> List.iter (fun p ->  
            printfn $"{p.Name} @ {p.PackagePrereleaseTag}"
        )

        if promptYesNo "OK?" then
            let source =
                "https://api.nuget.org/v3/index.json"

            let apikey =
                Environment.environVar "NUGET_KEY"

            if System.String.IsNullOrWhiteSpace apikey then
                failwith "NUGET_KEY environment variable is not set - cannot publish"

            for artifact in targets do
                let result =
                    DotNet.exec id "nuget" (sprintf "push -s %s -k %s %s --skip-duplicate" source apikey artifact)

                if not result.OK then
                    failwith "failed to push packages"
        else
            failwith "aborted"
    }

let releaseDocs =  BuildTask.create "ReleaseDocs" [buildDocs; publishNuget.IfNeeded] {
    let msg = 
        sprintf "release docs for version %s?" stableDocsVersionTag
    if promptYesNo msg then
        Shell.deleteDir "temp"
        Git.CommandHelper.runSimpleGitCommand "." (sprintf "clone %s temp/gh-pages --depth 1 -b gh-pages" projectRepo) |> ignore
        Shell.copyRecursive "output" "temp/gh-pages" true |> printfn "%A"
        Git.CommandHelper.runSimpleGitCommand "temp/gh-pages" "add ." |> printfn "%s"
        let cmd = sprintf """commit -a -m "Update generated documentation for version %s""" stableDocsVersionTag
        Git.CommandHelper.runSimpleGitCommand "temp/gh-pages" cmd |> printfn "%s"
        Git.Branches.push "temp/gh-pages"
    else failwith "aborted"
}

let prereleaseDocs =  BuildTask.create "PrereleaseDocs" [buildDocsPrerelease; publishNugetPrerelease.IfNeeded] {
    let msg = sprintf "release docs for version %s?" prereleaseTag
    if promptYesNo msg then
        Shell.deleteDir "temp"
        Git.CommandHelper.runSimpleGitCommand "." (sprintf "clone %s temp/gh-pages --depth 1 -b gh-pages" projectRepo) |> ignore
        Shell.copyRecursive "output" "temp/gh-pages" true |> printfn "%A"
        Git.CommandHelper.runSimpleGitCommand "temp/gh-pages" "add ." |> printfn "%s"
        let cmd = sprintf """commit -a -m "Update generated documentation for version %s""" prereleaseTag
        Git.CommandHelper.runSimpleGitCommand "temp/gh-pages" cmd |> printfn "%s"
        Git.Branches.push "temp/gh-pages"
    else failwith "aborted"
}
