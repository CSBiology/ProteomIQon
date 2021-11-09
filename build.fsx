#r "paket:
nuget BlackFox.Fake.BuildTask
nuget Fake.Core.Target
nuget Fake.Core.Process
nuget Fake.Core.ReleaseNotes
nuget Fake.IO.FileSystem
nuget Fake.DotNet.Cli
nuget Fake.DotNet.MSBuild
nuget Fake.DotNet.AssemblyInfoFile
nuget Fake.DotNet.Paket
nuget Fake.DotNet.FSFormatting
nuget Fake.DotNet.Fsi
nuget Fake.DotNet.NuGet
nuget Fake.Api.Github
nuget Fake.DotNet.Testing.Expecto 
nuget Fake.Tools.Git //"

#if !FAKE
#load "./.fake/build.fsx/intellisense.fsx"
#r "netstandard" // Temp fix for https://github.com/dotnet/fsharp/issues/5216
#endif

open BlackFox.Fake
open System.IO
open Fake.Core
open Fake.DotNet
open Fake.IO
open Fake.IO.FileSystemOperators
open Fake.IO.Globbing.Operators
open Fake.Tools

[<AutoOpen>]
/// user interaction prompts for critical build tasks where you may want to interrupt when you see wrong inputs.
module MessagePrompts =

    let prompt (msg:string) =
        System.Console.Write(msg)
        System.Console.ReadLine().Trim()
        |> function | "" -> None | s -> Some s
        |> Option.map (fun s -> s.Replace ("\"","\\\""))

    let rec promptYesNo msg =
        match prompt (sprintf "%s [Yn]: " msg) with
        | Some "Y" | Some "y" -> true
        | Some "N" | Some "n" -> false
        | _ -> System.Console.WriteLine("Sorry, invalid answer"); promptYesNo msg

    let promptProj (msg: string) =
        match prompt (sprintf "%s: " msg) with
        | Some x -> x
        | _ -> failwith ("Sorry, invalid answer")

    let releaseMsg = """This will stage all uncommitted changes, push them to the origin and bump the release version to the latest number in the RELEASE_NOTES.md file. 
        Do you want to continue?"""

    let releaseDocsMsg = """This will push the docs to gh-pages. Remember building the docs prior to this. Do you want to continue?"""

    let projMsg = """Name of the project you want to build"""

/// Executes a dotnet command in the given working directory
let runDotNet cmd workingDir =
    let result =
        DotNet.exec (DotNet.Options.withWorkingDirectory workingDir) cmd ""
    if result.ExitCode <> 0 then failwithf "'dotnet %s' failed in %s" cmd workingDir

/// Metadata about the project
module ProjectInfo = 

    let project = "ProteomIQon"

    let summary = "dotnet proteomics tools"

    let configuration = "Release"

    // Git configuration (used for publishing documentation in gh-pages branch)
    // The profile where the project is posted
    let gitOwner = "CSBiology"
    let gitName = "ProteomIQon"

    let gitHome = sprintf "%s/%s" "https://github.com" gitOwner

    let projectRepo = sprintf "%s/%s/%s" "https://github.com" gitOwner gitName

    let website = "/ProteomIQon"

    let pkgDir = "pkg"

    //let release = ReleaseNotes.load "RELEASE_NOTES.md"

    //let stableVersion = SemVer.parse release.NugetVersion

    //let stableVersionTag = (sprintf "%i.%i.%i" stableVersion.Major stableVersion.Minor stableVersion.Patch )

    let mutable prereleaseSuffix = ""

    let mutable prereleaseTag = ""

    let mutable isPrerelease = false

    let mutable projPattern = !! "src/**/*.*proj"

    let testProject = "tests/ProteomIQon.Tests/ProteomIQon.Tests.fsproj"

/// Barebones, minimal build tasks
module BasicTasks = 

    open ProjectInfo

    let setPrereleaseTag = BuildTask.create "SetPrereleaseTag" [] {
        printfn "Please enter pre-release package suffix"
        let suffix = System.Console.ReadLine()
        prereleaseSuffix <- suffix
        isPrerelease <- true
    }

    let clean = BuildTask.create "Clean" [] {
        !! "src/**/bin"
        ++ "src/**/obj"
        ++ "pkg"
        ++ "bin"
        |> Shell.cleanDirs 
    }

    let build = BuildTask.create "Build" [clean] {
        !! "src/**/*.*proj"
        |> Seq.iter (DotNet.build id)
    }

    let buildProj =
        BuildTask.create "BuildProj" [clean] {
            let proj =
                let rec loop (acc: IGlobbingPattern) =
                    if acc |> Seq.isEmpty then
                        printfn "Project doesn't exist. Try again."
                        let projName = promptProj projMsg
                        loop (!! (sprintf "src/**/%s.*proj" projName))
                    else
                        acc
                loop (!! (sprintf "src/**/%s.*proj" (promptProj projMsg)))
            projPattern <- proj
            projPattern
            |> Seq.iter (DotNet.build id)
        }


    let copyBinariesProj = BuildTask.create "CopyBinariesProj" [clean; buildProj] {
        let targets = 
            projPattern
            -- "src/**/*.shproj"
            |>  Seq.map (fun f -> ((Path.getDirectory f) </> "bin" </> configuration, "bin" </> (Path.GetFileNameWithoutExtension f)))
        targets
        |>  Seq.iter (fun (fromDir, toDir) -> Shell.copyDir toDir fromDir (fun _ -> true)
        )
    }

    let copyBinaries = BuildTask.create "CopyBinaries" [clean; build] {
        let targets = 
            !! "src/**/*.??proj"
            -- "src/**/*.shproj"
            |>  Seq.map (fun f -> ((Path.getDirectory f) </> "bin" </> configuration, "bin" </> (Path.GetFileNameWithoutExtension f)))
        for i in targets do printfn "%A" i
        targets
        |>  Seq.iter (fun (fromDir, toDir) -> Shell.copyDir toDir fromDir (fun _ -> true))
    }

/// Test executing build tasks
module TestTasks = 

    open ProjectInfo
    open BasicTasks

    let runTests = BuildTask.create "RunTests" [clean; build; copyBinaries] {
        let standardParams = Fake.DotNet.MSBuild.CliArguments.Create ()
        Fake.DotNet.DotNet.test(fun testParams ->
            {
                testParams with
                    Logger = Some "console;verbosity=detailed"
            }
        ) testProject
    }

    let runTestsProj = BuildTask.create "RunTestsProj" [clean; buildProj; copyBinariesProj] {
        let standardParams = Fake.DotNet.MSBuild.CliArguments.Create ()
        Fake.DotNet.DotNet.test(fun testParams ->
            {
                testParams with
                    Logger = Some "console;verbosity=detailed"
            }
        ) testProject
    }

/// Package creation
module PackageTasks = 

    open ProjectInfo

    open BasicTasks
    open TestTasks

    let pack = BuildTask.create "Pack" [clean; build; copyBinaries] {
        if promptYesNo (sprintf "creating stable package OK?") 
            then
                !! "src/**/*.*proj"
                |> Seq.iter (fun projPath -> 
                    let releaseNotePath = Path.Combine (Path.getDirectory projPath, "RELEASE_NOTES.md")
                    let release = ReleaseNotes.load releaseNotePath
                    let stableVersion = SemVer.parse release.NugetVersion
                    let stableVersionTag = (sprintf "%i.%i.%i" stableVersion.Major stableVersion.Minor stableVersion.Patch )
                    (Fake.DotNet.DotNet.pack (fun p ->
                        let msBuildParams =
                            {p.MSBuildParams with 
                                Properties = ([
                                    "Version",stableVersionTag
                                    "PackageReleaseNotes",  (release.Notes |> String.concat "\r\n")
                                ] @ p.MSBuildParams.Properties)
                            }
                        {
                            p with 
                                MSBuildParams = msBuildParams
                                OutputPath = Some pkgDir
                        }
                    )projPath
                    )
                )
        else failwith "aborted"
    }

    let packProj = BuildTask.create "PackProj" [clean; buildProj; copyBinariesProj] {
        if promptYesNo (sprintf "creating stable package OK?")
            then
                projPattern
                |> Seq.iter (fun projPath -> 
                    let releaseNotePath = Path.Combine (Path.getDirectory projPath, "RELEASE_NOTES.md")
                    let release = ReleaseNotes.load releaseNotePath
                    let stableVersion = SemVer.parse release.NugetVersion
                    let stableVersionTag = (sprintf "%i.%i.%i" stableVersion.Major stableVersion.Minor stableVersion.Patch )
                    (Fake.DotNet.DotNet.pack (fun p ->
                        let msBuildParams =
                            {p.MSBuildParams with 
                                Properties = ([
                                    "Version",stableVersionTag
                                    "PackageReleaseNotes",  (release.Notes |> String.concat "\r\n")
                                ] @ p.MSBuildParams.Properties)
                            }
                        {
                            p with 
                                MSBuildParams = msBuildParams
                                OutputPath = Some pkgDir
                        }
                    )projPath
                    )
                )
        else failwith "aborted"
        }

    let packPrerelease = BuildTask.create "PackPrerelease" [setPrereleaseTag; clean; build; copyBinaries] {
        if promptYesNo (sprintf "package suffix will be %s OK?" prereleaseSuffix )
            then 
                !! "src/**/*.*proj"
                //-- "src/**/Plotly.NET.Interactive.fsproj"
                |> Seq.iter (fun projPath -> 
                    let releaseNotePath = Path.Combine (Path.getDirectory projPath, "RELEASE_NOTES.md")
                    let release = ReleaseNotes.load releaseNotePath
                    prereleaseTag <- (sprintf "%s-%s" release.NugetVersion prereleaseSuffix)
                    (Fake.DotNet.DotNet.pack (fun p ->
                        let msBuildParams =
                            {p.MSBuildParams with 
                                Properties = ([
                                    "Version", prereleaseTag
                                    "PackageReleaseNotes",  (release.Notes |> String.toLines )
                                ] @ p.MSBuildParams.Properties)
                            }
                        {
                            p with 
                                VersionSuffix = Some prereleaseSuffix
                                OutputPath = Some pkgDir
                                MSBuildParams = msBuildParams
                        }
                ) projPath
                )
                )
        else
            failwith "aborted"
    }

    let packPrereleaseProj = BuildTask.create "PackPrereleaseProj" [setPrereleaseTag; clean; buildProj; copyBinariesProj] {
        if promptYesNo (sprintf "package tag will be %s OK?" prereleaseTag )
            then 
                projPattern
                |> Seq.iter (fun projPath -> 
                    let releaseNotePath = Path.Combine (Path.getDirectory projPath, "RELEASE_NOTES.md")
                    let release = ReleaseNotes.load releaseNotePath
                    prereleaseTag <- (sprintf "%s-%s" release.NugetVersion prereleaseSuffix)
                    (Fake.DotNet.DotNet.pack (fun p ->
                        let msBuildParams =
                            {p.MSBuildParams with 
                                Properties = ([
                                    "Version", prereleaseTag
                                    "PackageReleaseNotes",  (release.Notes |> String.toLines )
                                ] @ p.MSBuildParams.Properties)
                            }
                        {
                            p with 
                                VersionSuffix = Some prereleaseSuffix
                                OutputPath = Some pkgDir
                                MSBuildParams = msBuildParams
                        }
                ) projPath
                )
                )
        else
            failwith "aborted"
    }

    let pusProj = 
        BuildTask.create "PushNupkg" [] {
            let proj =
                let rec loop (acc: IGlobbingPattern) =
                    if acc |> Seq.isEmpty then
                        printfn "Package doesn't exist. Try again."
                        let projName = promptProj projMsg
                        loop (!! (sprintf "pkg/%s*.nupkg" projName))
                    else
                        acc
                loop (!! (sprintf "pkg/%s*.nupkg" (promptProj projMsg)))
            projPattern <- proj
            projPattern
            |> Seq.iter (fun pkgPath -> 
                let source = "https://api.nuget.org/v3/index.json"
                let apikey =  Environment.environVar "NUGET_KEY"
                let result = DotNet.exec id "nuget" (sprintf "push -s %s -k %s %s --skip-duplicate" source apikey pkgPath)
                if not result.OK then failwith "failed to push packages"
                )
        }
    //let packagedToolPath = Path.getFullName "./pkg/aglet/" 
    
    //Target.create "InstallLocalTool" (fun _ ->
    //    runDotNet "tool install --add-source nupkg aglet" packagedToolPath
    //)
open BasicTasks
BuildTask.runOrDefault copyBinaries