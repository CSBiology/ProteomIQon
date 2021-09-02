module Pipeline

open Expecto
open System
open System.IO
open ProteomIQon.Core.InputPaths
open Fake.DotNet
open BioFSharp.Mz
open System.Data.SQLite
open FSharpAux.IO

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
            let relToDirectory = getRelativePath Environment.CurrentDirectory
            let dbEstimate = relToDirectory "../../../data/PSMStatistics/in/MinimalEstimate.db"
            let dbFixed = relToDirectory "../../../data/PSMStatistics/in/MinimalFixed.db"
            let psmEstimate = relToDirectory "../../../data/PSMStatistics/in/minimalEstimate.psm"
            let psmFixed = relToDirectory "../../../data/PSMStatistics/in/minimalFixed.psm"
            let psmStatsParamsEstimate = relToDirectory "../../../data/PSMStatistics/in/pSMStatisticsParamsEstimate.json"
            let psmStatsParamsFixed = relToDirectory "../../../data/PSMStatistics/in/pSMStatisticsParamsFixed.json"
            let outDirectoryEstimate = relToDirectory "../../../data/PSMStatistics/out/estimateOut"
            let outDirectoryFixed = relToDirectory "../../../data/PSMStatistics/out/fixedOut"
            let psmStatsExe = relToDirectory "../../../../../bin/PSMStatistics/net5.0/ProteomIQon.PSMStatistics.dll"
            // run tool
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmStatsExe psmEstimate outDirectoryEstimate psmStatsParamsEstimate dbEstimate) Environment.CurrentDirectory
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" psmStatsExe psmFixed outDirectoryFixed psmStatsParamsFixed dbFixed) Environment.CurrentDirectory
            let referenceQPSMEstimate =
                let qpsmPath = relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalReference.qpsm"
                File.ReadAllLines qpsmPath
            let testQPSMEstimate =
                let qpsmPath = relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalEstimate.qpsm"
                File.ReadAllLines qpsmPath
            let referenceQPSMFixed =
                let qpsmPath = relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalReference.qpsm"
                File.ReadAllLines qpsmPath
            let testQPSMFixed =
                let qpsmPath = relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalFixed.qpsm"
                File.ReadAllLines qpsmPath
            let compare =
                referenceQPSMEstimate = testQPSMEstimate && referenceQPSMFixed = testQPSMFixed
            // cleanup
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalFixed.qpsm")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalFixed_log.txt")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/fixedOut/PSMStatistics_log.txt")
            Directory.Delete (relToDirectory "../../../data/PSMStatistics/out/fixedOut/minimalFixed_plots")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalEstimate.qpsm")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalEstimate_log.txt")
            File.Delete (relToDirectory "../../../data/PSMStatistics/out/estimateOut/PSMStatistics_log.txt")
            Directory.Delete (relToDirectory "../../../data/PSMStatistics/out/estimateOut/minimalEstimate_plots")
            Expect.isTrue compare "QPSMs are different"

        testCase "PSMBasedQuantification" <| fun _ ->
            let relToDirectory = getRelativePath Environment.CurrentDirectory
            let db = relToDirectory "../../../data/PSMBasedQuantification/in/Minimal.db"
            let qpsm = relToDirectory "../../../data/PSMBasedQuantification/in/minimal.qpsm"
            let mzlite = relToDirectory "../../../data/PSMBasedQuantification/in/minimal.mzlite"
            let quantParams = relToDirectory "../../../data/PSMBasedQuantification/in/QuantificationParams.json"
            let outDirectory = relToDirectory "../../../data/PSMBasedQuantification/out"
            let quantExe = relToDirectory "../../../../../bin/PSMBasedQuantification/net5.0/ProteomIQon.PSMBasedQuantification.dll"
            // cleanup
            try
                File.Delete (relToDirectory "../../../data/PSMBasedQuantification/out/minimal.quant")
                File.Delete (relToDirectory "../../../data/PSMBasedQuantification/out/minimal_log.txt")
                File.Delete (relToDirectory "../../../data/PSMBasedQuantification/out/PSMBasedQuantification_log.txt")
                Directory.Delete (relToDirectory "../../../data/PSMBasedQuantification/out/minimal_plots")
            with
            | _ -> ()
            // run tool
            runDotNet (sprintf "%s -i %s -ii %s -o %s -p %s -d %s -dc" quantExe mzlite qpsm outDirectory quantParams db) Environment.CurrentDirectory
            let referenceQuant =
                let quantPath = relToDirectory "../../../data/PSMBasedQuantification/out/minimalReference.quant"
                FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.QuantificationResult>().ReadFile(quantPath,'\t',false,1)
                |> Array.ofSeq
            let testQuant =
                let quantPath = relToDirectory "../../../data/PSMBasedQuantification/out/minimal.quant"
                FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.QuantificationResult>().ReadFile(quantPath,'\t',false,1)
                |> Array.ofSeq
            let compare,fields = 
                let unequalFields =
                    Array.map2 (fun (reference: ProteomIQon.Dto.QuantificationResult) (test: ProteomIQon.Dto.QuantificationResult) ->
                        [|
                            "StringSequence",
                                test.StringSequence = reference.StringSequence
                            "GlobalMod",
                                test.GlobalMod = reference.GlobalMod
                            "Charge",
                                test.Charge = reference.Charge
                            "PepSequenceID",
                                test.PepSequenceID = reference.PepSequenceID
                            "ModSequenceID",
                                test.ModSequenceID = reference.ModSequenceID
                            "PrecursorMZ",
                                test.PrecursorMZ >= reference.PrecursorMZ * 0.99 &&
                                test.PrecursorMZ <= reference.PrecursorMZ * 1.01
                            "MeasuredMass",
                                test.MeasuredMass >= reference.MeasuredMass * 0.99 &&
                                test.MeasuredMass <= reference.MeasuredMass * 1.01
                            "TheoMass"
                                ,test.TheoMass = reference.TheoMass
                            "AbsDeltaMass",
                                test.AbsDeltaMass >= reference.AbsDeltaMass * 0.7 &&
                                test.AbsDeltaMass <= reference.AbsDeltaMass * 1.3
                            "ProteinNames",
                                test.ProteinNames = reference.ProteinNames
                            "QuantMz_Light",
                                test.QuantMz_Light >= reference.QuantMz_Light * 0.99 &&
                                test.QuantMz_Light <= reference.QuantMz_Light * 1.01
                            "Quant_Light",
                                test.Quant_Light >= reference.Quant_Light * 0.99 &&
                                test.Quant_Light <= reference.Quant_Light * 1.01
                            "MeasuredApex_Light",
                                test.MeasuredApex_Light >= reference.MeasuredApex_Light * 0.99 &&
                                test.MeasuredApex_Light <= reference.MeasuredApex_Light * 1.01
                            //"Seo_Light",
                            //    test.Seo_Light >= reference.Seo_Light * 0.99 &&
                            //    test.Seo_Light <= reference.Seo_Light * 1.01
                            "Difference_SearchRT_FittedRT_Light",
                                abs test.Difference_SearchRT_FittedRT_Light >= abs reference.Difference_SearchRT_FittedRT_Light * 0.9 &&
                                abs test.Difference_SearchRT_FittedRT_Light <= abs reference.Difference_SearchRT_FittedRT_Light * 1.1
                            "KLDiv_CorrectedObserved_Theoretical_Light",
                                test.KLDiv_CorrectedObserved_Theoretical_Light >= reference.KLDiv_CorrectedObserved_Theoretical_Light * 0.9 &&
                                test.KLDiv_CorrectedObserved_Theoretical_Light <= reference.KLDiv_CorrectedObserved_Theoretical_Light * 1.1
                            "KLDiv_Observed_Theoretical_Light",
                                test.KLDiv_Observed_Theoretical_Light >= reference.KLDiv_Observed_Theoretical_Light * 0.99 &&
                                test.KLDiv_Observed_Theoretical_Light <= reference.KLDiv_Observed_Theoretical_Light * 1.01
                            "QuantMz_Heavy",
                                test.QuantMz_Heavy >= reference.QuantMz_Heavy * 0.99 &&
                                test.QuantMz_Heavy <= reference.QuantMz_Heavy * 1.01
                            "Quant_Heavy",
                                test.Quant_Heavy >= reference.Quant_Heavy * 0.99 &&
                                test.Quant_Heavy <= reference.Quant_Heavy * 1.01
                            "MeasuredApex_Heavy",
                                test.MeasuredApex_Heavy >= reference.MeasuredApex_Heavy * 0.99 &&
                                test.MeasuredApex_Heavy <= reference.MeasuredApex_Heavy * 1.01
                            //"Seo_Heavy",
                            //    test.Seo_Heavy >= reference.Seo_Heavy * 0.99 &&
                            //    test.Seo_Heavy <= reference.Seo_Heavy * 1.01
                            "Difference_SearchRT_FittedRT_Heavy",
                                abs test.Difference_SearchRT_FittedRT_Heavy >= abs reference.Difference_SearchRT_FittedRT_Heavy * 0.9 &&
                                abs test.Difference_SearchRT_FittedRT_Heavy <= abs reference.Difference_SearchRT_FittedRT_Heavy * 1.1
                            "KLDiv_CorrectedObserved_Theoretical_Heavy",
                                test.KLDiv_CorrectedObserved_Theoretical_Heavy >= reference.KLDiv_CorrectedObserved_Theoretical_Heavy * 0.9 &&
                                test.KLDiv_CorrectedObserved_Theoretical_Heavy <= reference.KLDiv_CorrectedObserved_Theoretical_Heavy * 1.1
                            "KLDiv_Observed_Theoretical_Heavy",
                                test.KLDiv_CorrectedObserved_Theoretical_Heavy >= reference.KLDiv_CorrectedObserved_Theoretical_Heavy * 0.9 &&
                                test.KLDiv_CorrectedObserved_Theoretical_Heavy <= reference.KLDiv_CorrectedObserved_Theoretical_Heavy * 1.1
                            "Correlation_Light_Heavy",
                                test.Correlation_Light_Heavy >= reference.Correlation_Light_Heavy * 0.99 &&
                                test.Correlation_Light_Heavy <= reference.Correlation_Light_Heavy * 1.01
                            "QuantificationSource",
                                test.QuantificationSource = reference.QuantificationSource
                        |]
                    ) referenceQuant testQuant
                    |> Array.collect (
                        Array.filter (fun (field,equal) -> equal = false)
                    )
                    |> Array.distinct
                    |> Array.map fst
                if unequalFields |> Array.isEmpty then true , ""
                else false, unequalFields |> String.concat ";"
            Expect.isTrue compare (sprintf "Quants are different in the following fields: %s" fields)
        testCase "ProteinInference" <| fun _ ->
            let subject = true
            Expect.isTrue subject "I compute, therefore I am."
    ]
