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

        testCase "QuantBasedAlignment" <| fun _ ->
            let relToDirectory = getRelativePath Environment.CurrentDirectory
            let target = relToDirectory "../../../data/QuantBasedAlignment/in"
            let sources = relToDirectory "../../../data/QuantBasedAlignment/in"
            let outDirectory = relToDirectory "../../../data/QuantBasedAlignment/out"
            let relQuantPath =
                if System.Runtime.InteropServices.RuntimeInformation.IsOSPlatform(System.Runtime.InteropServices.OSPlatform.Linux) then
                    "../../../../../bin/QuantBasedAlignment_linux-x64/net5.0/ProteomIQon.QuantBasedAlignment_linux-x64.dll"
                elif System.Runtime.InteropServices.RuntimeInformation.IsOSPlatform(System.Runtime.InteropServices.OSPlatform.Windows) then
                    "../../../../../bin/QuantBasedAlignment_win-x64/net5.0/ProteomIQon.QuantBasedAlignment_win_x64.dll"
                else    
                    failwith "not supported OS to test QuantBasedAlignment"                
            let quantExe = relToDirectory relQuantPath
            // cleanup
            try
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc1.align")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc1.alignmetric")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc1_log.txt")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/QuantBasedAlignment_log.txt")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc1_Metrics.html")   
                
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc2.align")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc2.alignmetric")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc2_log.txt")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc2_Metrics.html")   

                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc3.align")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc3.alignmetric")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc3_log.txt")
                File.Delete (relToDirectory "../../../data/QuantBasedAlignment/out/trunc3_Metrics.html")   
            with
            | _ -> ()
            // run tool
            runDotNet (sprintf "%s -i %s -ii %s -o %s" quantExe target sources outDirectory) Environment.CurrentDirectory
            let referenceAlign =
                let quantPath = relToDirectory "../../../data/QuantBasedAlignment/out/trunc1_Reference.align"
                FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.AlignmentResult>().ReadFile(quantPath,'\t',false,1)
                |> Array.ofSeq
            let testAlign =
                let quantPath = relToDirectory "../../../data/QuantBasedAlignment/out/trunc1.align"
                FSharpAux.IO.SchemaReader.Csv.CsvReader<ProteomIQon.Dto.AlignmentResult>().ReadFile(quantPath,'\t',false,1)
                |> Array.ofSeq
            let compare = 
                    referenceAlign
                    |> Array.filter (fun (reference: ProteomIQon.Dto.AlignmentResult) ->
                            testAlign
                            |> Array.exists (fun (test: ProteomIQon.Dto.AlignmentResult) ->
                                test.StringSequence = reference.StringSequence
                                test.GlobalMod = reference.GlobalMod
                                test.Charge = reference.Charge
                                test.PepSequenceID = reference.PepSequenceID
                                test.ModSequenceID = reference.ModSequenceID
                            )
                    )  
                    |> Array.length
                    |> (=) referenceAlign.Length           
            Expect.isTrue compare (sprintf "Quants are different in the following fields")


        testCase "ProteinInference" <| fun _ ->
            let relToDirectory = getRelativePath Environment.CurrentDirectory
            let db = relToDirectory "../../../data/ProteinInference/in/Minimal.db"
            let qpsm = relToDirectory "../../../data/ProteinInference/in/minimal.qpsm"
            let proteinInferenceParamsStorey = relToDirectory "../../../data/ProteinInference/in/ProteinInferenceParamsStorey.json"
            let proteinInferenceParamsMAYU = relToDirectory "../../../data/ProteinInference/in/ProteinInferenceParamsMAYU.json"
            let outDirectoryStorey = relToDirectory "../../../data/ProteinInference/out/storeyOut"
            let outDirectoryMAYU = relToDirectory "../../../data/ProteinInference/out/mayuOut"
            let protInfExe = relToDirectory "../../../../../bin/ProteinInference/net5.0/ProteomIQon.ProteinInference.dll"
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" protInfExe qpsm outDirectoryStorey proteinInferenceParamsStorey db) Environment.CurrentDirectory
            runDotNet (sprintf "%s -i %s -o %s -p %s -d %s" protInfExe qpsm outDirectoryMAYU proteinInferenceParamsMAYU db) Environment.CurrentDirectory
            let referenceProtStorey =
                let protPath = relToDirectory "../../../data/ProteinInference/out/storeyOut/minimalReference.prot"
                File.ReadAllLines protPath
            let protStorey =
                let protPath = relToDirectory "../../../data/ProteinInference/out/storeyOut/minimal.prot"
                File.ReadAllLines protPath
            let referenceProtMAYU =
                let protPath = relToDirectory "../../../data/ProteinInference/out/mayuOut/minimalReference.prot"
                File.ReadAllLines protPath
            let protMAYU =
                let protPath = relToDirectory "../../../data/ProteinInference/out/mayuOut/minimal.prot"
                File.ReadAllLines protPath
            let compare = referenceProtStorey = protStorey && referenceProtMAYU = protMAYU
            // cleanup
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/minimal.prot")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/ProteinInference_createClassItemCollection_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/ProteinInference_inferProteins_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/ProteinInference_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/ProteinInference_readAndInferFile_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/storeyOut/QValueGraph.html")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/minimal.prot")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/ProteinInference_createClassItemCollection_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/ProteinInference_inferProteins_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/ProteinInference_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/ProteinInference_readAndInferFile_log.txt")
            File.Delete (relToDirectory "../../../data/ProteinInference/out/mayuOut/QValueGraph.html")
            Expect.isTrue compare "Prots are different"
        testCase "LabelFreeProteinQuantification" <| fun _ ->
            let relToDirectory = getRelativePath Environment.CurrentDirectory
            let quantAndProt = relToDirectory "../../../data/LabelFreeProteinQuantification/in/minimal.quantAndProt"
            let labelFreeQuantificationParams = "../../../data/LabelFreeProteinQuantification/in/LabelFreeQuantificationParams.json"
            let labelFreeQuantificationParamsChargeAgg = "../../../data/LabelFreeProteinQuantification/in/LabelFreeQuantificationParams_ChargeAgg.json"
            let labelFreeQuantificationParamsChargeAggModAgg = "../../../data/LabelFreeProteinQuantification/in/LabelFreeQuantificationParams_ChargeAgg_ModAgg.json"
            let labelFreeQuantificationParamsTransformFilterSum = "../../../data/LabelFreeProteinQuantification/in/LabelFreeQuantificationParams_Transform_Filter_Sum.json"
            let outDirectory = "../../../data/LabelFreeProteinQuantification/out/normal"
            let outDirectoryChargeAgg = "../../../data/LabelFreeProteinQuantification/out/chargeAgg"
            let outDirectoryChargeAggModAgg = "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg"
            let outDirectoryTransformFilterSum = "../../../data/LabelFreeProteinQuantification/out/transformFilterSum"
            let labelFreeExe = relToDirectory "../../../../../bin/LabelFreeProteinQuantification/net5.0/ProteomIQon.LabelFreeProteinQuantification.dll"
            runDotNet (sprintf "%s -i %s -o %s -p %s" labelFreeExe quantAndProt outDirectory labelFreeQuantificationParams) Environment.CurrentDirectory
            runDotNet (sprintf "%s -i %s -o %s -p %s" labelFreeExe quantAndProt outDirectoryChargeAgg labelFreeQuantificationParamsChargeAgg) Environment.CurrentDirectory
            runDotNet (sprintf "%s -i %s -o %s -p %s" labelFreeExe quantAndProt outDirectoryChargeAggModAgg labelFreeQuantificationParamsChargeAggModAgg) Environment.CurrentDirectory
            runDotNet (sprintf "%s -i %s -o %s -p %s" labelFreeExe quantAndProt outDirectoryTransformFilterSum labelFreeQuantificationParamsTransformFilterSum) Environment.CurrentDirectory
            let referenceLabelFree =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/minimalReference.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/minimalProteinReference.txt"
                File.ReadAllLines labelFreeProtein
            let labelFree =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/LabelFreeQuant.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/ProteinAggregation.txt"
                File.ReadAllLines labelFreeProtein
            let referenceLabelFreeChargeAgg =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/minimalReference.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/minimalProteinReference.txt"
                File.ReadAllLines labelFreeProtein,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/minimalChargeReference.txt"
                File.ReadAllLines labelFreeCharge
            let labelFreeChargeAgg =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/LabelFreeQuant.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/ProteinAggregation.txt"
                File.ReadAllLines labelFreeProtein,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/ChargeAggregation.txt"
                File.ReadAllLines labelFreeCharge
            let referenceLabelFreeChargeAggModAgg =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/minimalReference.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/minimalProteinReference.txt"
                File.ReadAllLines labelFreeProtein,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/minimalChargeReference.txt"
                File.ReadAllLines labelFreeCharge,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/minimalModificationReference.txt"
                File.ReadAllLines labelFreeCharge
            let labelFreeChargeAggModAgg =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/LabelFreeQuant.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ProteinAggregation.txt"
                File.ReadAllLines labelFreeProtein,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ChargeAggregation.txt"
                File.ReadAllLines labelFreeCharge,
                let labelFreeCharge = relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ModificationAggregation.txt"
                File.ReadAllLines labelFreeCharge
            let referenceLabelFreeTransformFilterSum =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/minimalReference.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/minimalProteinReference.txt"
                File.ReadAllLines labelFreeProtein
            let labelFreeTransformFilterSum =
                let labelFree = relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/LabelFreeQuant.txt"
                File.ReadAllLines labelFree,
                let labelFreeProtein = relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/ProteinAggregation.txt"
                File.ReadAllLines labelFreeProtein
            let compare = referenceLabelFree = labelFree && referenceLabelFreeChargeAgg = labelFreeChargeAgg && referenceLabelFreeChargeAggModAgg = labelFreeChargeAggModAgg && referenceLabelFreeTransformFilterSum = labelFreeTransformFilterSum
            // cleanup
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/LabelFreeQuant.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/normal/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/LabelFreeQuant.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAgg/ChargeAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/LabelFreeQuant.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ChargeAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/chargeAggModAgg/ModificationAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/LabelFreeQuant.txt")
            File.Delete (relToDirectory "../../../data/LabelFreeProteinQuantification/out/transformFilterSum/ProteinAggregation.txt")
            Expect.isTrue compare "Output files are not identical"

        testCase "LabeledProteinQuantification" <| fun _ ->
            let relToDirectory = getRelativePath Environment.CurrentDirectory
            let quantAndProt = relToDirectory "../../../data/LabeledProteinQuantification/in/minimal.quantAndProt"
            let labeledQuantificationParams = "../../../data/LabeledProteinQuantification/in/LabeledQuantificationParams.json"
            let labeledQuantificationParamsCorrChargeAgg = "../../../data/LabeledProteinQuantification/in/LabeledQuantificationParams_CorrFilter_ChargeAgg.json"
            let labeledQuantificationParamsCorrChargeAggModAgg = "../../../data/LabeledProteinQuantification/in/LabeledQuantificationParams_CorrFilter_ChargeAgg_ModAgg.json"
            let labeledQuantificationParamsTransformFilterSum = "../../../data/LabeledProteinQuantification/in/LabeledQuantificationParams_Transform_Filter_Sum.json"
            let labeledQuantificationParamsCorr = "../../../data/LabeledProteinQuantification/in/LabeledQuantificationParams_withCorrFilter.json"
            let outDirectory = "../../../data/LabeledProteinQuantification/out/normal"
            let outDirectoryCorrChargeAgg = "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg"
            let outDirectoryCorrChargeAggModAgg = "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg"
            let outDirectoryTransformFilterSum = "../../../data/LabeledProteinQuantification/out/transformFilterSum"
            let outDirectoryCorr = "../../../data/LabeledProteinQuantification/out/corrFilter"
            let labeledExe = relToDirectory "../../../../../bin/LabeledProteinQuantification/net5.0/ProteomIQon.LabeledProteinQuantification.dll"
            runDotNet (sprintf "%s -i %s -o %s -p %s" labeledExe quantAndProt outDirectory labeledQuantificationParams) Environment.CurrentDirectory
            runDotNet (sprintf "%s -i %s -o %s -p %s" labeledExe quantAndProt outDirectoryCorrChargeAgg labeledQuantificationParamsCorrChargeAgg) Environment.CurrentDirectory
            runDotNet (sprintf "%s -i %s -o %s -p %s" labeledExe quantAndProt outDirectoryCorrChargeAggModAgg labeledQuantificationParamsCorrChargeAggModAgg) Environment.CurrentDirectory
            runDotNet (sprintf "%s -i %s -o %s -p %s" labeledExe quantAndProt outDirectoryTransformFilterSum labeledQuantificationParamsTransformFilterSum) Environment.CurrentDirectory
            runDotNet (sprintf "%s -i %s -o %s -p %s" labeledExe quantAndProt outDirectoryCorr labeledQuantificationParamsCorr) Environment.CurrentDirectory
            let referenceLabeled =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/minimalReference.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/minimalProteinReference.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/minimalGlobModReference.txt"
                File.ReadAllLines labeledGlobMod
            let labeled =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/LabeledQuant.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/ProteinAggregation.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/normal/GlobModAggregation.txt"
                File.ReadAllLines labeledGlobMod
            let referenceLabeledCorrChargeAgg =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/minimalReference.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/minimalProteinReference.txt"
                File.ReadAllLines labeledProtein,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/minimalChargeReference.txt"
                File.ReadAllLines labeledCharge,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/minimalGlobModReference.txt"
                File.ReadAllLines labeledGlobMod
            let labeledCorrChargeAgg =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/LabeledQuant.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/ProteinAggregation.txt"
                File.ReadAllLines labeledProtein,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/ChargeAggregation.txt"
                File.ReadAllLines labeledCharge,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/GlobModAggregation.txt"
                File.ReadAllLines labeledGlobMod
            let referenceLabeledCorrChargeAggModAgg =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/minimalReference.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/minimalProteinReference.txt"
                File.ReadAllLines labeledProtein,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/minimalChargeReference.txt"
                File.ReadAllLines labeledCharge,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/minimalModificationReference.txt"
                File.ReadAllLines labeledCharge,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/minimalGlobModReference.txt"
                File.ReadAllLines labeledGlobMod
            let labeledCorrChargeAggModAgg =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/LabeledQuant.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ProteinAggregation.txt"
                File.ReadAllLines labeledProtein,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ChargeAggregation.txt"
                File.ReadAllLines labeledCharge,
                let labeledCharge = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ModificationAggregation.txt"
                File.ReadAllLines labeledCharge,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/GlobModAggregation.txt"
                File.ReadAllLines labeledGlobMod
            let referenceLabeledTransformFilterSum =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/minimalReference.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/minimalProteinReference.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/minimalGlobModReference.txt"
                File.ReadAllLines labeledGlobMod
            let labeledTransformFilterSum =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/LabeledQuant.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/ProteinAggregation.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/GlobModAggregation.txt"
                File.ReadAllLines labeledGlobMod
            let referenceLabeledCorr =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/minimalReference.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/minimalProteinReference.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/minimalGlobModReference.txt"
                File.ReadAllLines labeledGlobMod
            let labeledCorr =
                let labeled = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/LabeledQuant.txt"
                File.ReadAllLines labeled,
                let labeledProtein = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/ProteinAggregation.txt"
                File.ReadAllLines labeledProtein,
                let labeledGlobMod = relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/GlobModAggregation.txt"
                File.ReadAllLines labeledGlobMod
            let compare = referenceLabeled = labeled && referenceLabeledCorrChargeAgg = labeledCorrChargeAgg && referenceLabeledCorrChargeAggModAgg = labeledCorrChargeAggModAgg && referenceLabeledTransformFilterSum = labeledTransformFilterSum && referenceLabeledCorr = labeledCorr
            // cleanup
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/normal/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/normal/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/normal/LabeledQuant.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/normal/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/normal/GlobModAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/LabeledQuant.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilter/GlobModAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/LabeledQuant.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/GlobModAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAgg/ChargeAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/LabeledQuant.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/GlobModAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ChargeAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/corrFilterChargeAggModAgg/ModificationAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/LabeledProteinQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/LabeledQuantification_log.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/LabeledQuant.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/ProteinAggregation.txt")
            File.Delete (relToDirectory "../../../data/LabeledProteinQuantification/out/transformFilterSum/GlobModAggregation.txt")
            Expect.isTrue compare "Output files are not identical"
        testCase "AddDeducedPeptides" <| fun _ ->
            let relToDirectory = getRelativePath Environment.CurrentDirectory
            let quantDirectory = relToDirectory "../../../data/AddDeducedPeptides/in/quant"
            let protDirectory = relToDirectory "../../../data/AddDeducedPeptides/in/prot"
            let outDirectory = relToDirectory "../../../data/AddDeducedPeptides/out"
            let addDeducedPeptidesExe = relToDirectory "../../../../../bin/AddDeducedPeptides/net5.0/ProteomIQon.AddDeducedPeptides.dll"
            // run tool
            runDotNet (sprintf "%s -i %s -ii %s -o %s" addDeducedPeptidesExe quantDirectory protDirectory outDirectory) Environment.CurrentDirectory
            let reference1 =
                let qpsmPath = relToDirectory "../../../data/AddDeducedPeptides/out/minimalReference1.prot"
                File.ReadAllLines qpsmPath
            let test1 =
                let qpsmPath = relToDirectory "../../../data/AddDeducedPeptides/out/minimalQuant1.prot"
                File.ReadAllLines qpsmPath
            let reference2 =
                let qpsmPath = relToDirectory "../../../data/AddDeducedPeptides/out/minimalReference2.prot"
                File.ReadAllLines qpsmPath
            let test2 =
                let qpsmPath = relToDirectory "../../../data/AddDeducedPeptides/out/minimalQuant2.prot"
                File.ReadAllLines qpsmPath

            let compare =
                reference1 = test1 && reference2 = test2
            // cleanup
            File.Delete (relToDirectory "../../../data/AddDeducedPeptides/out/minimalQuant1.prot")
            File.Delete (relToDirectory "../../../data/AddDeducedPeptides/out/minimalQuant2.prot")
            File.Delete (relToDirectory "../../../data/AddDeducedPeptides/out/AddDeducedPeptides_log.txt")
            Expect.isTrue compare "Deduced Proteins are different"
    ]
