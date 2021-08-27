module Pipeline

open Expecto
open System
open System.IO
open ProteomIQon.Core.InputPaths
open Fake.DotNet
open BioFSharp.Mz
open System.Data.SQLite

let runDotNet cmd workingDir =
    let result =
        DotNet.exec (DotNet.Options.withWorkingDirectory workingDir) cmd ""
    if result.ExitCode <> 0 then failwithf "'dotnet %s' failed in %s" cmd workingDir

let selectAllModSequence cn =
    let querystring = "SELECT * FROM ModSequence"
    use cmd = new SQLiteCommand(querystring, cn)
    use reader = cmd.ExecuteReader()
    [
        while reader.Read() do 
            yield (reader.GetValue(0), reader.GetValue(1),reader.GetValue(2), reader.GetValue(3), reader.GetValue(4), reader.GetValue(5))
    ]

let selectAllCleavageIndex cn =
    let querystring = "SELECT * FROM CleavageIndex"
    use cmd = new SQLiteCommand(querystring, cn)
    use reader = cmd.ExecuteReader()
    [
        while reader.Read() do 
            yield (reader.GetValue(0), reader.GetValue(1),reader.GetValue(2), reader.GetValue(3), reader.GetValue(4), reader.GetValue(5))
    ]

let selectAllPepSequence cn =
    let querystring = "SELECT * FROM PepSequence"
    use cmd = new SQLiteCommand(querystring, cn)
    use reader = cmd.ExecuteReader()
    [
        while reader.Read() do 
            yield (reader.GetValue(0), reader.GetValue(1))
    ]

let selectAllProtein cn =
    let querystring = "SELECT * FROM Protein"
    use cmd = new SQLiteCommand(querystring, cn)
    use reader = cmd.ExecuteReader()
    [
        while reader.Read() do 
            yield (reader.GetValue(0), reader.GetValue(1))
    ]

[<Tests>]
let pipelineTests =
    testList "Pipeline" [
        testCase "PeptideDB" <| fun _ ->
            let relToDirectory = getRelativePath Environment.CurrentDirectory
            let fastaPath = relToDirectory "../../../data/PeptideDB/in/example.fasta"
            let peptideDBParams = relToDirectory "../../../data/PeptideDB/in/peptideDBParams.json"
            let outDirectory = relToDirectory "../../../data/PeptideDB/out/"
            let pepDBExe = relToDirectory "../../../../../bin/PeptideDB/net5.0/ProteomIQon.PeptideDB.dll"
            // run tool
            runDotNet (sprintf "%s -i %s -o %s -p %s" pepDBExe fastaPath outDirectory peptideDBParams) Environment.CurrentDirectory
            let referenceDB =
                let dbPath = relToDirectory "../../../data/PeptideDB/out/MinimalReference.db"
                use cn = SearchDB.getDBConnection dbPath
                let res = selectAllCleavageIndex cn, selectAllModSequence cn, selectAllPepSequence cn, selectAllProtein cn
                res
            let testDB =
                let dbPath = relToDirectory "../../../data/PeptideDB/out/Minimal.db"
                use cn = SearchDB.getDBConnection dbPath
                let res = selectAllCleavageIndex cn, selectAllModSequence cn, selectAllPepSequence cn, selectAllProtein cn
                res
            let compare = referenceDB = testDB
            // cleanup
            File.Delete (relToDirectory "../../../data/PeptideDB/out/Minimal.db")
            File.Delete (relToDirectory "../../../data/PeptideDB/out/PeptideDB_log.txt")
            File.Delete (relToDirectory "../../../data/PeptideDB/out/PeptideDB_Minimal_log.txt")
            Expect.isTrue compare "Peptide databases are different"
    
        testCase "PeptideSpectrumMatching" <| fun _ ->
            let relToDirectory = getRelativePath Environment.CurrentDirectory
            let db = relToDirectory "../../../data/PeptideSpectrumMatching/in/Minimal.db"
            let mzlite = relToDirectory "../../../data/PeptideSpectrumMatching/in/minimal.mzlite"
            let psmParams = relToDirectory "../../../data/PeptideSpectrumMatching/in/defaultParams.json"
            let outDirectory = relToDirectory "../../../data/PeptideSpectrumMatching/out/"
            let psmExe = relToDirectory "../../../../../bin/PeptideSpectrumMatching/net5.0/ProteomIQon.PeptideSpectrumMatching.dll"
            // run tool
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmExe mzlite outDirectory psmParams db) Environment.CurrentDirectory
            let referencePSM = 
                let psmPath = relToDirectory "../../../data/PeptideSpectrumMatching/out/minimalReference.psm"
                File.ReadAllLines psmPath
            let testPSM = 
                let psmPath = relToDirectory "../../../data/PeptideSpectrumMatching/out/minimal.psm"
                File.ReadAllLines psmPath
            let compare = referencePSM = testPSM
            // cleanup
            File.Delete (relToDirectory "../../../data/PeptideSpectrumMatching/out/minimal.psm")
            File.Delete (relToDirectory "../../../data/PeptideSpectrumMatching/out/minimal_log.txt")
            File.Delete (relToDirectory "../../../data/PeptideSpectrumMatching/out/PeptideSpectrumMatching_log.txt")
            Expect.isTrue compare "PSMs are different"

        testCase "PSMStatistics" <| fun _ ->
            let subject = true
            Expect.isTrue subject "I compute, therefore I am."

        testCase "PSMBasedQuantification" <| fun _ ->
            let subject = true
            Expect.isTrue subject "I compute, therefore I am."

        testCase "ProteinInference" <| fun _ ->
            let subject = true
            Expect.isTrue subject "I compute, therefore I am."
    ]
