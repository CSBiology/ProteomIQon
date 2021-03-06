namespace ProteomIQon

open System
open MzIO
open FSharp.Stats
open FSharpAux
open ProteomIQon.Domain
open MzTABAux

module MzTABFormater =

    let protHeader path (mzTABParams: Domain.MzTABParams) =
        let expCount = mzTABParams.ExperimentNames.Length
        let stVarCount = mzTABParams.StudyVariables.Length
        let searchEnginecount = mzTABParams.SearchEngineNamesProt.Length
        let bestSearchEngineScore =
            formatOne searchEnginecount (sprintf "best_search_engine_score[%i]")
        let searchEngineScoreMS =
            formatTwo searchEnginecount expCount (sprintf "search_engine_score[%i]_ms_run[%i]")
        let psmsMSRun =
            formatOne expCount (sprintf "num_psms_ms_run[%i]")
        let pepsDistMSRun =
            formatOne expCount (sprintf "num_peptides_distinct_ms_run[%i]")
        let pepsUniqueMSRun =
            formatOne expCount (sprintf "num_peptides_unique_ms_run[%i]")
        let protAbundanceAssay =
            if mzTABParams.Labeled then
                formatOne (expCount*2) (sprintf "protein_abundance_assay[%i]")
            else
                formatOne expCount (sprintf "protein_abundance_assay[%i]")
        let protAbundanceStudVarParts =
            formatOneForThree 
                stVarCount
                (sprintf "protein_abundance_study_variable[%i]")
                (sprintf "protein_abundance_stdev_study_variable[%i]")
                (sprintf "protein_abundance_std_error_study_variable[%i]")
        let sb = new Text.StringBuilder()
        sb.AppendFormat(
            "PRH\taccession\tdescription\ttaxid\tspecies\tdatabase\tdatabase_version\tsearch_engine\t{0}\t{1}\treliability\t{2}\t{3}\t{4}\tambiguity_members\tmodifications\turi\tgo_terms\tprotein_coverage\t{5}\t{6}\topt_global_protein_group",
            bestSearchEngineScore,
            searchEngineScoreMS,
            psmsMSRun,
            pepsDistMSRun,
            pepsUniqueMSRun,
            protAbundanceAssay,
            protAbundanceStudVarParts
        ) |> ignore
        sb.AppendLine() |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let pepHeader path (mzTABParams: Domain.MzTABParams) =
        let expCount = mzTABParams.ExperimentNames.Length
        let stVarCount = mzTABParams.StudyVariables.Length
        let searchEnginecount = mzTABParams.SearchEngineNamesPep.Length
        let bestSearchEngineScore =
            formatOne searchEnginecount (sprintf "best_search_engine_score[%i]")
        let searchEngineScoreMS =
            formatTwo searchEnginecount expCount (sprintf "search_engine_score[%i]_ms_run[%i]")
        let pepAbundanceAssay =
            if mzTABParams.Labeled then
                formatOne (expCount*2) (sprintf "peptide_abundance_assay[%i]")
            else
                formatOne expCount (sprintf "peptide_abundance_assay[%i]")
        let pepAbundanceStudVarParts =
            formatOneForThree
                stVarCount
                (sprintf "peptide_abundance_study_variable[%i]")
                (sprintf "peptide_abundance_stdev_study_variable[%i]")
                (sprintf "peptide_abundance_std_error_study_variable[%i]")
        let sb = new Text.StringBuilder()
        sb.AppendFormat(
            "PEH\tsequence\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\t{0}\t{1}\treliability\tmodifications\tretention_time\tretention_time_window\tcharge\tmass_to_charge\turi\tspectra_ref\t{2}\t{3}\topt_global_global_mod\topt_global_protein_group",
            bestSearchEngineScore,
            searchEngineScoreMS,
            pepAbundanceAssay,
            pepAbundanceStudVarParts
        ) |> ignore
        sb.AppendLine() |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let psmHeader path (mzTABParams: Domain.MzTABParams) =
        let searchEnginecount = mzTABParams.SearchEngineNamesPSM.Length
        let searchEngineScoreMS =
            formatOne searchEnginecount (sprintf "search_engine_score[%i]")
        let sb = new Text.StringBuilder()
        sb.AppendFormat(
            "PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\t{0}\treliability\tmodifications\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\turi\tspectra_ref\tpre\tpost\tstart\tend\topt_global_protein_group",
            searchEngineScoreMS
        ) |> ignore
        sb.AppendLine() |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let protBody path (proteinS: ProteinSection[]) (sameProteinAcc: Map<string[],string>) =
        let sb = new Text.StringBuilder()
        printfn "%A" sameProteinAcc
        proteinS
        |> Array.map (fun prot ->
            sb.AppendFormat(
                "PRT\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}",
                (
                    let number = 
                        sameProteinAcc
                        |> Map.tryFind prot.accession
                    if number.IsSome then
                        prot.accession.[0] + number.Value
                    else
                        prot.accession.[0]
                )
                ,
                prot.description,
                prot.taxid,
                prot.species,
                prot.database,
                prot.database_version,
                prot.search_engine,
                //needs to be adapted for different search engines
                if prot.best_search_engine_score.IsSome then
                    string prot.best_search_engine_score.Value
                else
                    "null"
                ,
                prot.search_engine_score_ms_run
                |> Array.map (fun score ->
                    score
                    |> concatRuns "null"
                )
                |> String.concat "\t",
                //prot.reliability,
                "null",
                prot.num_psms_ms_run
                |> concatRuns "null",
                prot.num_peptides_distinct_ms_run
                |> concatRuns "null",
                prot.num_peptides_unique_ms_run
                |> Array.map (fun x -> "null")
                |> String.concat "\t",
                prot.ambiguity_members
                |> String.concat ",",
                prot.modifications,
                prot.uri,
                prot.go_terms
                |> String.concat "|",
                //prot.protein_coverage,
                "null",
                prot.protein_abundance_assay
                |> concatRuns "null",
                [|
                    prot.protein_abundance_study_variable;
                    prot.protein_abundance_stdev_study_variable;
                    prot.protein_abundance_std_error_study_variable
                |]
                |> concatRunsForDifferent "null",
                prot.accession
                |> String.concat ";"
            ) |> ignore
            sb.AppendLine()
        ) |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let pepBody path (peptideS: PeptideSection[]) (sameProteinAcc: Map<string[],string>) =
        let sb = new Text.StringBuilder()
        peptideS
        |> Array.map (fun pep ->
            pep.accession
            |> Array.iter (fun prot ->
                sb.AppendFormat(
                    "PEP\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}",
                    pep.sequence,
                    (
                        let number = 
                            sameProteinAcc
                            |> Map.tryFind prot
                        if number.IsSome then
                            prot.[0] + number.Value
                        else
                            prot.[0]
                    )
                    ,
                    pep.unique,
                    pep.database,
                    pep.database_version,
                    pep.search_engine,
                    //needs to be adapted for different search engines
                    pep.best_search_engine_score
                    |> Array.map (string << snd)
                    |> String.concat "\t"
                    ,
                    pep.search_engine_score_ms_run
                    |> Array.map (fun score ->
                        score
                        |> concatRuns "null"
                    )
                    |> String.concat "\t",
                    //pep.reliability,
                    "null",
                    pep.modifications,
                    pep.retention_time
                    |> string,
                    //pep.retention_time_window,
                    "null",
                    pep.charge,
                    pep.mass_to_charge,
                    pep.uri,
                    pep.spectra_ref,
                    pep.peptide_abundance_assay
                    |> concatRuns "null",
                    [|
                        pep.peptide_abundance_study_variable
                        pep.peptide_abundance_stdev_study_variable
                        pep.peptide_abundance_std_error_study_variable
                    |]
                    |> concatRunsForDifferent "null",
                    pep.labeling.toParam,
                    prot
                    |> String.concat ";"
                ) |> ignore
                sb.AppendLine()
                |> ignore
            )
        ) |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let psmBody path (psmS: PSMSection[]) (sameProteinAcc: Map<string[],string>) =
        let sb = new Text.StringBuilder()
        psmS
        |> Array.iteri (fun i psm ->
            sb.AppendFormat(
                "PSM\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}",
                psm.sequence,
                psm.PSM_ID,
                (
                    let number = 
                        sameProteinAcc
                        |> Map.tryFind psm.accession
                    if number.IsSome then
                        psm.accession.[0] + number.Value
                    else
                        psm.accession.[0]
                )
                ,
                psm.unique,
                psm.database,
                psm.database_version,
                psm.search_engine,
                //needs to be adapted for different search engines
                psm.search_engine_score
                |> Array.map string
                |> String.concat "\t",
                //psm.reliability,
                "null",
                psm.modifications,
                psm.retention_time,
                psm.charge,
                psm.exp_mass_to_charge,
                psm.calc_mass_to_charge,
                psm.uri,
                psm.spectra_ref,
                psm.pre,
                psm.post,
                //psm.start,
                //psm.ending
                "null",
                "null",
                psm.accession
                |> String.concat ";"
            ) |> ignore
            sb.AppendLine()
            |> ignore
            if i%10000 = 0 then
                IO.File.AppendAllText(path, sb.ToString())
                sb.Clear() |> ignore
        ) |> ignore
        IO.File.AppendAllText(path, sb.ToString())