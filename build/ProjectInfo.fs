module ProjectInfo

open Fake.Core


/// Contains relevant information about a project (e.g. version info, project location)
type ProjectInfo = {
    Name: string
    ProjFile: string
    ReleaseNotes: ReleaseNotes.ReleaseNotes Option
    PackageVersionTag: string
    mutable PackagePrereleaseTag: string
    AssemblyVersion: string
    AssemblyInformationalVersion: string
} with 
    /// creates a ProjectInfo given a name, project file path, and release notes file path.
    /// version info is created from the version header of the uppermost release notes entry.
    /// Assembly version is set to X.0.0, where X is the major version from the releas enotes.
    static member create(
        name: string,
        projFile: string,
        releaseNotesPath: string
    ): ProjectInfo = 
        let release = releaseNotesPath |> ReleaseNotes.load
        let stableVersion = release.NugetVersion |> SemVer.parse
        let stableVersionTag = $"{stableVersion.Major}.{stableVersion.Minor}.{stableVersion.Patch}"
        let assemblyVersion = $"{stableVersion.Major}.0.0"
        let assemblyInformationalVersion = stableVersionTag
        {
            Name = name
            ProjFile = projFile
            ReleaseNotes = Some release
            PackagePrereleaseTag = ""
            PackageVersionTag = stableVersionTag
            AssemblyVersion = assemblyVersion
            AssemblyInformationalVersion = assemblyInformationalVersion
        }    
    static member create(
        name: string,
        projFile: string
    ): ProjectInfo = 
        {
            Name = name
            ProjFile = projFile
            ReleaseNotes = None
            PackagePrereleaseTag = ""
            PackageVersionTag = ""
            AssemblyVersion = ""
            AssemblyInformationalVersion = ""
        }

// adapt this to reflect the core project in your repository. The only effect this will have is the version displayed in the docs, as it is currently only possible to have one version displayed there.
let CoreProject = ProjectInfo.create("ProteomIQon", "src/ProteomIQon/ProteomIQon.fsproj", "src/ProteomIQon/RELEASE_NOTES.md")

let projects = 
    [
        // add relative paths (from project root) to your projects here, including individual reslease notes files
        // e.g. ProjectInfo.create("MyProject", "src/MyProject/MyProject.fsproj", "src/MyProject/RELEASE_NOTES.md")
        CoreProject
        ProjectInfo.create("TableSort", "src/TableSort/TableSort.fsproj", "src/TableSort/RELEASE_NOTES.md")
        ProjectInfo.create("AddDeducedPeptides", "src/AddDeducedPeptides/AddDeducedPeptides.fsproj", "src/AddDeducedPeptides/RELEASE_NOTES.md")
        ProjectInfo.create("AlignmentBasedQuantification", "src/AlignmentBasedQuantification/AlignmentBasedQuantification.fsproj", "src/AlignmentBasedQuantification/RELEASE_NOTES.md")
        ProjectInfo.create("AlignmentBasedQuantStatistics", "src/AlignmentBasedQuantStatistics/AlignmentBasedQuantStatistics.fsproj", "src/AlignmentBasedQuantStatistics/RELEASE_NOTES.md")
        ProjectInfo.create("ConsensusSpectralLibrary", "src/ConsensusSpectralLibrary/ConsensusSpectralLibrary.fsproj", "src/ConsensusSpectralLibrary/RELEASE_NOTES.md")
        ProjectInfo.create("JoinQuantPepIonsWithProteins", "src/JoinQuantPepIonsWithProteins/JoinQuantPepIonsWithProteins.fsproj", "src/JoinQuantPepIonsWithProteins/RELEASE_NOTES.md")
        ProjectInfo.create("LabeledProteinQuantification", "src/LabeledProteinQuantification/LabeledProteinQuantification.fsproj", "src/LabeledProteinQuantification/RELEASE_NOTES.md")
        ProjectInfo.create("LabelEfficiencyCalculator", "src/LabelEfficiencyCalculator/LabelEfficiencyCalculator.fsproj", "src/LabelEfficiencyCalculator/RELEASE_NOTES.md")
        ProjectInfo.create("LabelFreeProteinQuantification", "src/LabelFreeProteinQuantification/LabelFreeProteinQuantification.fsproj", "src/LabelFreeProteinQuantification/RELEASE_NOTES.md")
        ProjectInfo.create("MzliteToMzML", "src/MzliteToMzML/MzliteToMzML.fsproj", "src/MzliteToMzML/RELEASE_NOTES.md")
        ProjectInfo.create("MzMLToMzLite", "src/MzMLToMzLite/MzMLToMzLite.fsproj", "src/MzMLToMzLite/RELEASE_NOTES.md")
        ProjectInfo.create("MzMLToMzLiteIonMobility", "src/MzMLToMzLiteIonMobility/MzMLToMzLiteIonMobility.fsproj", "src/MzMLToMzLiteIonMobility/RELEASE_NOTES.md")
        ProjectInfo.create("MzTAB", "src/MzTAB/MzTAB.fsproj", "src/MzTAB/RELEASE_NOTES.md")
        ProjectInfo.create("PeptideDB", "src/PeptideDB/PeptideDB.fsproj", "src/PeptideDB/RELEASE_NOTES.md")
        ProjectInfo.create("PeptideSpectrumMatching", "src/PeptideSpectrumMatching/PeptideSpectrumMatching.fsproj", "src/PeptideSpectrumMatching/RELEASE_NOTES.md")
        ProjectInfo.create("Preprocessing", "src/Preprocessing/Preprocessing.fsproj", "src/Preprocessing/RELEASE_NOTES.md")
        ProjectInfo.create("ProteinInference", "src/ProteinInference/ProteinInference.fsproj", "src/ProteinInference/RELEASE_NOTES.md")
        ProjectInfo.create("PSMBasedQuantification", "src/PSMBasedQuantification/PSMBasedQuantification.fsproj", "src/PSMBasedQuantification/RELEASE_NOTES.md")
        ProjectInfo.create("PSMBasedQuantificationTIMs", "src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fsproj", "src/PSMBasedQuantificationTIMs/RELEASE_NOTES.md")
        ProjectInfo.create("PSMStatistics", "src/PSMStatistics/PSMStatistics.fsproj", "src/PSMStatistics/RELEASE_NOTES.md")
        ProjectInfo.create("QuantBasedAlignment", "src/QuantBasedAlignment/QuantBasedAlignment.fsproj", "src/QuantBasedAlignment/RELEASE_NOTES.md")
        ProjectInfo.create("QuantBasedAlignment_linux-x64", "src/QuantBasedAlignment_linux-x64/QuantBasedAlignment_linux-x64.fsproj", "src/QuantBasedAlignment_linux-x64/RELEASE_NOTES.md")
        ProjectInfo.create("QuantBasedAlignment_win-x64", "src/QuantBasedAlignment_win-x64/QuantBasedAlignment_win-x64.fsproj", "src/QuantBasedAlignment_win-x64/RELEASE_NOTES.md")
        ProjectInfo.create("SpectralLibrary", "src/SpectralLibrary/SpectralLibrary.fsproj", "src/SpectralLibrary/RELEASE_NOTES.md")
        ProjectInfo.create("SWATHAnalysis", "src/SWATHAnalysis/SWATHAnalysis.fsproj", "src/SWATHAnalysis/RELEASE_NOTES.md")
    ]


let project = "ProteomIQon"

let testProjects = 
    [
        // add relative paths (from project root) to your testprojects here
        // e.g. ProjectInfo.create("MyTestProject", "tests/MyTestProject/MyTestProject.fsproj")
        ProjectInfo.create("ProteomIQon.Tests", "tests/ProteomIQon.Tests.fsproj")
    ]

let solutionFile  = $"{project}.sln"

let configuration = "Release"

let gitOwner = "CSBiology"

let gitHome = $"https://github.com/{gitOwner}"

let projectRepo = $"https://github.com/{gitOwner}/{project}"

let pkgDir = "pkg"


/// docs are always targeting the version of the core project
let stableDocsVersionTag = CoreProject.PackageVersionTag

/// branch tag is always the version of the core project
let branchTag = CoreProject.PackageVersionTag

/// prerelease suffix used by prerelease buildtasks
let mutable prereleaseSuffix = ""

/// prerelease tag used by prerelease buildtasks
let mutable prereleaseTag = ""

/// mutable switch used to signal that we are building a prerelease version, used in prerelease buildtasks
let mutable isPrerelease = false
