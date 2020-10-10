namespace ProteomIQon

open System
open System.Linq
open MzIO
open FSharp.Stats
open FSharpAux
open FSharpAux.IO
open FSharpAux.IO.SchemaReader.Attribute
open FSharp.Plotly
open ProteomIQon.Dto
open ProteomIQon.Domain

module MzTAB =

    type TableSort =
        {
            [<FieldAttribute("\"Key #1\"")>]
            Protein             : string
            [<FieldAttribute("\"Key #2\"")>]
            Experiment          : string
            DistinctPeptideCount: float
            Quant_Heavy         : float
            Quant_Light         : float
            QValue              : float
            Ratio               : float
            Ratio_CV            : float
            Ratio_SEM           : float
            Ratio_StDev         : float
        }

    type InferredProteinClassItemOut =
        {
            GroupOfProteinIDs: string
            PeptideSequence  : string
            Class            : string
            TargetScore      : float
            DecoyScore       : float
            QValue           : float
        }

    type MetaDataSection =
        {
            mzTab_version                    : string
            mzTab_mode                       : string
            mzTab_type                       : string
            mzTab_ID                         : (string)option
            title                            : (string)option
            description                      : string
            sample_processing                : ((string*int)[])option
            instrument_name                  : ((string*int)[])option
            instrument_source                : ((string*int)[])option
            instrument_analyzer              : ((string[]*int)[])option
            instrument_detector              : ((string*int)[])option
            software                         : ((string*int)[])option
            software_setting                 : ((string[]*int)[])option
            protein_search_engine_score      : ((string*string*int)[])option
            peptide_search_engine_score      : ((string*string*int)[])option
            psm_search_engine_score          : ((string*string*int)[])option
            false_discovery_rate             : (string[])option
            publication                      : (string[])option
            contact_name                     : ((string*int)[])option
            contact_affiliation              : ((string*int)[])option
            contact_email                    : ((string*int)[])option
            uri                              : (string)option
            fixed_mod                        : (string*int)[]
            fixed_mod_site                   : ((string*int)[])option
            fixed_mod_position               : ((string*int)[])option
            variable_mod                     : (string*int)[]
            variable_mod_site                : ((string*int)[])option
            variable_mod_position            : ((string*int)[])option
            quantification_method            : (string)option
            protein_quantification_unit      : (string)option
            peptide_quantification_unit      : (string)option
            ms_run_format                    : ((string*int)[])option
            ms_run_location                  : (string*int)[]
            ms_run_id_format                 : ((string*int)[])option
            ms_run_fragmentation_method      : ((string[]*int)[])option
            ms_run_hash                      : ((string*int)[])option
            ms_run_hash_method               : ((string*int)[])option
            custom                           : (string[])option
            sample_species                   : ((string[]*int)[])option
            sample_tissue                    : ((string[]*int)[])option
            sample_cell_type                 : ((string[]*int)[])option
            sample_disease                   : ((string[]*int)[])option
            sample_description               : ((string*int)[])option
            sample_custom                    : ((string[]*int)[])option
            assay_quantification_reagent     : ((string*int)[])option
            assay_quantification_mod         : (((string*int)[]*int)[])option
            assay_quantification_mod_site    : (((string*int)[]*int)[])option
            assay_quantification_mod_position: (((string*int)[]*int)[])option
            assay_sample_ref                 : ((string*int)[])option
            assay_ms_run_ref                 : ((string*int)[])option
            study_variable_assay_refs        : ((int[]*int)[])option
            study_variable_sample_refs       : ((int[]*int)[])option
            study_variable_description       : ((string*int)[])option
            cv_label                         : ((string*int)[])option
            cv_full_name                     : ((string*int)[])option
            cv_version                       : ((string*int)[])option
            cv_url                           : ((string*int)[])option
            colunit_protein                  : (string)option
            colunit_psm                      : (string)option
            colunit_small_molecule           : (string)option
        }

    type ProteinSection =
        {
            accession                                 : string
            description                               : string
            taxid                                     : int
            species                                   : string
            database                                  : string
            database_version                          : string
            search_engine                             : string
            best_search_engine_score                  : float
            search_engine_score_ms_run                : (int*float option)[][]
            reliability                               : int
            num_psms_ms_run                           : (int*float option)[]
            num_peptides_distinct_ms_run              : (int*float option)[]
            num_peptides_unique_ms_run                : (int*float option)[]
            ambiguity_members                         : string[]
            modifications                             : string
            uri                                       : string
            go_terms                                  : string[]
            protein_coverage                          : float
            protein_abundance_assay                   : (int*float option)[]
            protein_abundance_study_variable          : (int*float option)[]
            protein_abundance_stdev_study_variable    : (int*float option)[]
            protein_abundance_std_error_study_variable: (int*float option)[]
        }

    type PeptideSection =
        {
            sequence                                  : string
            accession                                 : string[]
            unique                                    : int
            database                                  : string
            database_version                          : string
            search_engine                             : string
            best_search_engine_score                  : float
            search_engine_score_ms_run                : (int*float option)[][]
            reliability                               : int
            modifications                             : string
            retention_time                            : (float option * float option)
            retention_time_window                     : float*float
            charge                                    : int
            mass_to_charge                            : float
            uri                                       : string
            spectra_ref                               : string
            peptide_abundance_assay                   : (int*float option)[]
            peptide_abundance_study_variable          : (int*float option)[]
            peptide_abundance_stdev_study_variable    : (int*float option)[]
            peptide_abundance_std_error_study_variable: (int*float option)[]
        }

    type PSMSection =
        {
            sequence           : string
            PSM_ID             : int
            accession          : string[]
            unique             : int
            database           : string
            database_version   : string
            search_engine      : string
            search_engine_score: float[]
            reliability        : int
            modifications      : string
            retention_time     : float
            charge             : int
            exp_mass_to_charge : float
            calc_mass_to_charge: float
            uri                : string
            spectra_ref        : string
            pre                : string
            post               : string
            start              : int
            ending             : int
        }

    type AlignedComplete =
        {
            ProtInf  : InferredProteinClassItemOut
            Qpsm     : PSMStatisticsResult
            Quant    : QuantificationResult
            TableSort: TableSort
        }

    let createAlignedComplete protInf qpsm quant tableSort =
        {
            ProtInf   = protInf
            Qpsm      = qpsm
            Quant     = quant
            TableSort = tableSort
        }

    let getFilePaths path identifier =
        IO.Directory.GetFiles (path, identifier)

    let readQpsm (path: string) =
        SeqIO.Seq.fromFileWithCsvSchema<PSMStatisticsResult>(path, '\t', true)
        |> Seq.toArray

    let readQuant (path: string) =
        SeqIO.Seq.fromFileWithCsvSchema<QuantificationResult>(path, '\t', true)
        |> Seq.toArray

    let readProt (path: string) =
        SeqIO.Seq.fromFileWithCsvSchema<InferredProteinClassItemOut>(path, '\t', true)
        |> Seq.toArray
        |> Array.map (fun x ->
            {
                x with
                    GroupOfProteinIDs =
                        x.GroupOfProteinIDs
                        |> String.split ';'
                        |> Seq.map (fun y ->
                            y
                            |> String.replace "\"" ""
                        )
                        |> Seq.sort
                        |> String.concat ";"
            }
        )

    let readTab (path: string) =
        SeqIO.Seq.fromFileWithCsvSchema<TableSort>(path, '\t', true)
        |> Seq.toArray
        |> Array.map (fun x ->
            {
                x with
                    Protein =
                        x.Protein
                        |> String.split ';'
                        |> Seq.map (fun y ->
                            y
                            |> String.replace "\"" ""
                        )
                        |> Seq.sort
                        |> String.concat ";"
                    Experiment =
                        x.Experiment
                        |> String.replace "\"" ""
            }
        )
        |> Array.groupBy (fun x -> x.Experiment)

    let allignAllFiles (tabFile:string) (protFiles: string[]) (quantFiles: string[]) (qpsmFiles: string[]) =
        let tab = readTab tabFile
        let prot =
            protFiles
            |> Array.map (fun x ->
                IO.Path.GetFileNameWithoutExtension x, readProt x
            )
        let quant =
            quantFiles
            |> Array.map (fun x ->
                IO.Path.GetFileNameWithoutExtension x, readQuant x
            )
        let qpsm =
            qpsmFiles
            |> Array.map (fun x ->
                IO.Path.GetFileNameWithoutExtension x, readQpsm x
            )
        let alignedProtTab =
            tab
            |> Array.map (fun (exp,tabs) ->
                // looks for the prot file from the same experiment as the tab group
                let _,correspondingProts =
                    prot
                    |> Array.find (fun (exp',_) -> exp'=exp)
                // maps over the tab entries of that experiment and finds the corresponding prot entries
                tabs
                |> Array.map (fun x ->
                    let corrProt =
                        correspondingProts
                        |> Array.find (fun y -> y.GroupOfProteinIDs = x.Protein)
                    {|TableSort = x; ProtInf = corrProt|}
                )
            )
        // assigns the peptides that were used for the identification of the proteins to the proteins, and in a subsequent step, the psms that point to the peptides
        let alignedProtTabQuantQpsm =
            alignedProtTab
            |> Array.map (fun exp ->
                // check which experiment the current group belongs to
                let experiment = exp.[0].TableSort.Experiment
                // looks for the quant file from the same experiment as the tab group
                let _,correspondingQuant =
                    quant
                    |> Array.find (fun (exp',_) -> exp'=experiment)
                let _,correspondingQpsm =
                    qpsm
                    |> Array.find (fun (exp',_) -> exp'=experiment)
                exp
                |> Array.map (fun alignedInfo ->
                    let peptides =
                        alignedInfo.ProtInf.PeptideSequence
                        |> String.split ';'
                    // filter for all peptides that were used in the inference of the protein group
                    let corrQuant =
                        correspondingQuant
                        |> Array.filter (fun x -> peptides.Contains x.StringSequence)
                    // filter for all psms that match the quantified peptide
                    let corrQuantWithQpsm =
                        corrQuant
                        |> Array.map (fun quantItem ->
                            let psmItem =
                                correspondingQpsm
                                |> Array.filter (fun x -> x.ModSequenceID = quantItem.ModSequenceID && x.Charge = quantItem.Charge)
                            // create an entry for each psm with the peptide quantification it belongs to
                            psmItem
                            |> Array.map (fun psmI ->
                                {|Quant = quantItem; Qpsm = psmI|}
                            )
                        )
                        |> Array.concat
                    // map over all psms with their corresponding peptide and add the protein info that it points to
                    // results in multiple entries per protein/peptide, distinguishable by peptide/psm
                    corrQuantWithQpsm
                    |> Array.map (fun corrQuantQpsm ->
                        createAlignedComplete alignedInfo.ProtInf corrQuantQpsm.Qpsm corrQuantQpsm.Quant alignedInfo.TableSort
                    )
                )
                |> Array.concat
            )
        alignedProtTabQuantQpsm
        |> Array.concat

    let findValueNumberedProt (expNames: (string*int)[]) (proteins: TableSort []) (fieldName: string) =
        let fieldFunc (tableSort: TableSort) =
            let res = ReflectionHelper.tryGetPropertyValue tableSort fieldName
            match res with
            |None -> failwith (sprintf "Field %s doesn't exist" fieldName)
            |Some x -> x
        expNames
        |> Array.map (fun (experiment,number) ->
            let value =
                proteins
                |> Array.tryFind (fun prot -> prot.Experiment = experiment)
                |> fun x ->
                    match x with
                    |None -> None
                    |Some y ->
                        Some ((fieldFunc y) :?> float)
            number, value
        )

    let findValueNumberedPep (expNames: (string*int)[]) (peptides: (QuantificationResult*string) []) (fieldName: string) =
        let fieldFunc (tableSort: QuantificationResult) =
            let res = ReflectionHelper.tryGetPropertyValue tableSort fieldName
            match res with
            |None -> failwith (sprintf "Field %s doesn't exist" fieldName)
            |Some x -> x
        expNames
        |> Array.map (fun (experiment,number) ->
            let value =
                peptides
                |> Array.tryFind (fun (pep,exp) -> exp = experiment)
                |> fun x ->
                    match x with
                    |None -> None
                    |Some y ->
                        Some ((fieldFunc (fst y)) :?> float)
            number, value
        )

    let fieldFuncPSM (psm: PSMStatisticsResult) (fieldName: string) =
        let res = ReflectionHelper.tryGetPropertyValue psm fieldName
        match res with
        |None -> failwith (sprintf "Field %s doesn't exist" fieldName)
        |Some x -> x :?> float

    let concatRuns (fillEmpty: string) (data: (int*float option)[]) =
        data
        |> Array.sortBy fst
        |> Array.map (fun (run,value) ->
            match value with
            | None -> fillEmpty
            | Some s ->
                if System.Double.IsNaN s then
                    "null"
                else
                    string s
        )
        |> String.concat "\t"

    let sortAndPick (sortBy: 'T -> 'Key) (pick: 'T -> 'C) (input: 'T[]) =
        input
        |> Array.sortBy sortBy
        |> Array.map pick

    let formatOneMD (strF: int -> string -> string) (input: string[]) =
        [|
            for i=0 to (input.Length - 1) do
                yield strF (i+1) input.[i]
        |]
        |> String.concat "\n"

    let formatTwoMD (strF: int -> int -> string -> string) (input: string[][])=
        [|
            for i=0 to (input.Length - 1) do
                for j=0 to (input.[i].Length - 1) do
                    yield strF (i+1) (j+1) input.[i].[j]
        |]
        |> String.concat "\n"

    let matchOption (format: 'a -> string) (sb: Text.StringBuilder) (item: 'a option) =
        match item with
        |Some x -> 
            x
            |> format
            |> fun x -> sb.AppendLine(x)
        |None -> sb

    let metaDataSection (md: MetaDataSection) =
        let sb = new Text.StringBuilder()
        let mzTabVersion =
            sb.AppendFormat("MTD\tmzTab-version\t{0}", md.mzTab_version)
            |> ignore
            sb.AppendLine()
        let mzTabMode =
            sb.AppendFormat("MTD\tmzTab-mode\t{0}", md.mzTab_mode)
            |> ignore
            sb.AppendLine()
        let mzTabType =
            sb.AppendFormat("MTD\tmzTab-type\t{0}", md.mzTab_type)
            |> ignore
            sb.AppendLine()
        let mzTabID =
            md.mzTab_ID
            |> matchOption (sprintf "MTD\tmzTab-ID\t%s") sb
        let title =
            md.title
            |> matchOption (sprintf "MTD\ttitle\t%s") sb
        let description =
            sb.AppendFormat("MTD\tdescription\t{0}", md.description)
            |> ignore
            sb.AppendLine()
        let sampleProcessing =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tsample_processing[%i]\t%s")
                pickF >> strF
            md.sample_processing
            |> matchOption f sb
        let instrumentName =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tinstrument[%i]-name\t%s")
                pickF >> strF
            md.instrument_name
            |> matchOption f sb
        let instrumentSource =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tinstrument[%i]-source\t%s")
                pickF >> strF
            md.instrument_source
            |> matchOption f sb
        let instrumentAnalyzer =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatTwoMD (sprintf "MTD\tinstrument[%i]-analyzer[%i]\t%s")
                pickF >> strF
            md.instrument_analyzer
            |> matchOption f sb
        let instrumentDetector =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tinstrument[%i]-detector\t%s")
                pickF >> strF
            md.instrument_detector
            |> matchOption f sb
        let software =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tsoftware[%i]\t%s")
                pickF >> strF
            md.software
            |> matchOption f sb
        let softwareSetting =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    // documentation isn't clear about the number behind setting
                    // sample shows it without number, but according to instructions there should be a number
                    formatTwoMD (sprintf "MTD\tsoftware[%i]-setting[%i]\t%s")
                pickF >> strF
            md.software_setting
            |> matchOption f sb
        let protSearchEngineScore =
            let f =
                let pickF =
                    sortAndPick (fun (_,_,x)->x) (fun (y,_,_)->y)
                let strF =
                    formatOneMD (sprintf "MTD\tprotein_search_engine_score[%i]\t%s")
                pickF >> strF
            md.protein_search_engine_score
            |> matchOption f sb
        let pepSearchEngineScore =
            let f =
                let pickF =
                    sortAndPick (fun (_,_,x)->x) (fun (y,_,_)->y)
                let strF =
                    formatOneMD (sprintf "MTD\tpeptide_search_engine_score[%i]\t%s")
                pickF >> strF
            md.peptide_search_engine_score
            |> matchOption f sb
        let psmSearchEngineScore =
            let f =
                let pickF =
                    sortAndPick (fun (_,_,x)->x) (fun (y,_,_)->y)
                let strF =
                    formatOneMD (sprintf "MTD\tpsm_search_engine_score[%i]\t%s")
                pickF >> strF
            md.psm_search_engine_score
            |> matchOption f sb
        let fdr =
            let f =
                (String.concat "|") >> sprintf "MTD\tfalse_discovery_rate\t%s"
            md.false_discovery_rate
            |> matchOption f sb
        let publication =
            let strF =
                formatOneMD (sprintf "MTD\tpsm_search_engine_score[%i]\t%s")
            md.publication
            |> matchOption strF sb
        printfn "%s" (sb.ToString())

    let proteinSection (allAligned: AlignedComplete[]) (mzTABParams: Domain.MzTABParams) =
        let experimentNames = mzTABParams.ExperimentNames
        let groupedTab =
            allAligned
            |> Array.groupBy (fun x -> x.TableSort.Protein)
            |> Array.map snd
            |> Array.map (fun x ->
                x
                |> Array.sortBy (fun y -> y.TableSort.Experiment)
            )
            |> Array.map (fun x ->
                x
                |> Array.groupBy (fun y -> y.TableSort)
            )
        groupedTab
        |> Array.map (fun protGroup ->
            let studyVars =
                let quant =
                    findValueNumberedProt experimentNames (protGroup |> Array.map fst) "Ratio"
                    |> Array.sortBy fst
                mzTABParams.StudyVariables
                |> Array.map (fun (name,assays,number) ->
                    let corrQuant =
                        quant
                        |> Array.filter (fun x -> assays.Contains (fst x))
                    number,
                    corrQuant
                    |> Array.choose snd
                )

            {
                accession                                 = (fst protGroup.[0]).Protein
                description                               = ""
                taxid                                     = 3055
                species                                   = "Chlamydomonas reinhardtii"
                database                                  = "Chlamy.db"
                database_version                          = "19-Apr-20 21:44"
                search_engine                             =
                    mzTABParams.SearchEngineNamesProt

                    |> Array.map (fun (x,_,_) -> x)
                    |> String.concat "|"
                best_search_engine_score                  =
                    protGroup
                    |> Array.maxBy (fun (prot,peps) -> prot.QValue)
                    |> fun (prot,peps) -> prot.QValue
                search_engine_score_ms_run                =
                    mzTABParams.SearchEngineNamesProt
                    |> Array.map (fun (searchengine,fieldName,number) ->
                        findValueNumberedProt experimentNames (protGroup |> Array.map fst) fieldName
                        |> Array.sortBy fst
                    )
                reliability                               = 3
                num_psms_ms_run                           =
                    experimentNames
                    |> Array.map (fun (experiment,number) ->
                        let prot =
                            protGroup
                            |> Array.tryFind (fun (prot,psms) ->
                                prot.Experiment = experiment
                            )
                        match prot with
                        |None -> number, None
                        |Some (prot,psms) -> number, Some (float psms.Length)
                    )
                num_peptides_distinct_ms_run              =
                    findValueNumberedProt experimentNames (protGroup |> Array.map fst) "DistinctPeptideCount"
                    |> Array.sortBy fst
                // TODO: check uniqueness for all peptides
                num_peptides_unique_ms_run                =
                    findValueNumberedProt experimentNames (protGroup |> Array.map fst) "DistinctPeptideCount"
                    |> Array.sortBy fst
                ambiguity_members                         = [|"No decision yet"|]
                modifications                             = null
                uri                                       = @"C:\Users\jonat\source\repos\mzTAB\Chlamy.db"
                go_terms                                  = [||]
                // remember sequence
                protein_coverage                          = 0.
                protein_abundance_assay                   =
                    findValueNumberedProt experimentNames (protGroup |> Array.map fst) "Ratio"
                    |> Array.sortBy fst
                protein_abundance_study_variable          =
                    studyVars
                    |> Array.map (fun (number,variables) ->
                        number,
                        match variables with
                        | [||] -> None
                        | _ ->
                            Some (
                                variables
                                |> Array.average
                            )
                    )
                protein_abundance_stdev_study_variable    =
                    studyVars
                    |> Array.map (fun (number,variables) ->
                        number,
                        match variables with
                        | [||] -> None
                        | _ ->
                            Some (
                                variables
                                |> Seq.stDev
                            )
                    )
                protein_abundance_std_error_study_variable=
                    studyVars
                    |> Array.map (fun (number,variables) ->
                        number,
                        match variables with
                        | [||] -> None
                        | _ ->
                            Some (
                                variables
                                |> fun x ->
                                    let stDev = Seq.stDev x
                                    stDev / (sqrt (float x.Length))
                            )
                    )
            }
        )

    let peptideSection (allAligned: AlignedComplete[]) (mzTABParams: Domain.MzTABParams): PeptideSection[] =
        let experimentNames = mzTABParams.ExperimentNames
        let groupedTab =
            allAligned
            |> Array.groupBy (fun pep -> pep.Quant.StringSequence, pep.Quant.Charge)
            |> Array.map snd
            |> Array.map (fun x ->
                x
                |> Array.sortBy (fun y -> y.TableSort.Experiment)
            )
            |> Array.map (fun x ->
                x
                |> Array.groupBy (fun y -> y.Quant)
            )
        groupedTab
        |> Array.map (fun pepGroup ->
            let forF =
                pepGroup
                |> Array.map (fun x ->
                    fst x,
                    // this is only here as a check if my sorting isn't messed up
                    (snd x)
                    |> Array.map (fun y -> y.TableSort.Experiment)
                    |> Array.distinct
                    |> fun z ->
                        if z.Length <> 1 then
                            failwith "unexpected experiment number"
                        else
                            z.[0]
                )
            let corrProt =
                pepGroup
                |> Array.collect (fun (peptide,rest) ->
                    rest
                    |> Array.map (fun x -> x.TableSort.Protein)
                )
                |> Array.distinct
            let abundanceAssay =
                let light =
                    findValueNumberedPep experimentNames forF "Quant_Light"
                    |> Array.sortBy fst
                let heavy =
                    findValueNumberedPep experimentNames forF "Quant_Heavy"
                    |> Array.sortBy fst
                Array.map2 (fun (i,l) (j,h) ->
                    match l,h with
                    | Some x, Some y ->
                        printfn "%f,%f" x y
                        i, Some (x/y)
                    | _, None -> i, None
                    | None, _ -> i, None
                ) light heavy
            let studyVars =
                mzTABParams.StudyVariables
                |> Array.map (fun (name,assays,number) ->
                    let corrQuant =
                        abundanceAssay
                        |> Array.filter (fun x -> assays.Contains (fst x))
                    number,
                    corrQuant
                    |> Array.choose snd
                )
            {
                sequence                                  =
                    pepGroup
                    |> Array.map (fun x -> (fst x).StringSequence)
                    |> Array.distinct
                    |> fun z ->
                        if z.Length <> 1 then
                            failwith "unexpected experiment number"
                        else
                            z.[0]
                accession                                 =
                    corrProt
                unique                                    =
                    if corrProt.Length > 1 then
                        0
                    else
                        1
                database                                  =
                    "Chlamy.db"
                database_version                          =
                    "19-Apr-20 21:44"
                search_engine                             =
                    mzTABParams.SearchEngineNamesPep
                    |> Array.map (fun (x,_,_) -> x)
                    |> String.concat "|"
                best_search_engine_score                  =
                    pepGroup
                    |> Array.maxBy (fun (peptide,rest) -> peptide.MeanPercolatorScore)
                    |> fun (peptide,rest) -> peptide.MeanPercolatorScore
                search_engine_score_ms_run                =
                    mzTABParams.SearchEngineNamesPep
                    |> Array.map (fun (searchengine,fieldName,number) ->
                        findValueNumberedPep experimentNames forF fieldName
                        |> Array.sortBy fst
                    )
                reliability                               =
                    3
                modifications                             =
                    ""
                retention_time                            =
                    pepGroup
                    |> Array.sortBy (fun x -> (fst x).MeanPercolatorScore)
                    |> Array.head
                    |> fun (peptide,rest) ->
                        let labeled =
                            rest
                            |> Array.filter (fun x ->
                                x.Qpsm.GlobalMod = 1
                            )
                            |> fun x ->
                                if x.Length > 0 then
                                    Some (x |> Array.averageBy (fun y -> y.Qpsm.ScanTime))
                                else
                                    None
                        let unlabeled =
                            rest
                            |> Array.filter (fun x ->
                                x.Qpsm.GlobalMod = 0
                            )
                            |> fun x ->
                                if x.Length > 0 then
                                    Some (x |> Array.averageBy (fun y -> y.Qpsm.ScanTime))
                                else
                                    None
                        labeled,unlabeled
                retention_time_window                     =
                    1.,2.
                charge                                    =
                    (fst pepGroup.[0]).Charge
                mass_to_charge                            =
                    pepGroup
                    |> Array.sortBy (fun x -> (fst x).MeanPercolatorScore)
                    |> Array.head
                    |> fun (peptide,rest) -> peptide.PrecursorMZ
                uri                                       =
                    ""
                spectra_ref                               =
                    "null"
                peptide_abundance_assay                   =
                    abundanceAssay
                peptide_abundance_study_variable          =
                    studyVars
                    |> Array.map (fun (number,variables) ->
                        number,
                        match variables with
                        | [||] -> None
                        | _ ->
                            Some (
                                variables
                                |> Array.average
                            )
                    )
                peptide_abundance_stdev_study_variable    =
                    studyVars
                    |> Array.map (fun (number,variables) ->
                        number,
                        match variables with
                        | [||] -> None
                        | _ ->
                            Some (
                                variables
                                |> Seq.stDev
                            )
                    )
                peptide_abundance_std_error_study_variable=
                    studyVars
                    |> Array.map (fun (number,variables) ->
                        number,
                        match variables with
                        | [||] -> None
                        | _ ->
                            Some (
                                variables
                                |> fun x ->
                                    let stDev = Seq.stDev x
                                    stDev / (sqrt (float x.Length))
                            )
                    )
            }
        )

    let psmSection (allAligned: AlignedComplete[]) (mzTABParams: Domain.MzTABParams): PSMSection[] =
        let groupedTab =
            allAligned
            |> Array.groupBy (fun x -> x.Qpsm)
        groupedTab
            |> Array.map (fun (psm,rest) ->
            let corrProt =
                rest
                |> Array.map (fun x -> x.TableSort.Protein)
                |> Array.distinct
            {
                sequence                                  =
                    psm.StringSequence
                PSM_ID = psm.ScanNr
                accession                                 =
                    corrProt
                unique                                    =
                    if corrProt.Length <> 1 then
                        0
                    else
                        1
                database                                  =
                    "Chlamy.db"
                database_version                          =
                    "19-Apr-20 21:44"
                search_engine                             =
                    mzTABParams.SearchEngineNamesPSM
                    |> Array.map (fun (x,_,_) -> x)
                    |> String.concat "|"
                search_engine_score                       =
                    mzTABParams.SearchEngineNamesPSM
                    |> Array.map (fun (searchengine,fieldName,number) ->
                        fieldFuncPSM psm fieldName
                    )
                reliability                               =
                    3
                modifications                             =
                    null
                retention_time                            =
                    psm.ScanTime
                charge                                    =
                    psm.Charge
                exp_mass_to_charge                            =
                    psm.PrecursorMZ
                calc_mass_to_charge                            =
                    BioFSharp.Mass.toMZ psm.TheoMass (float psm.Charge)
                uri                                       =
                    "null"
                spectra_ref                               =
                    psm.PSMId
                pre="null"
                post="null"
                start=1
                ending=1
            }
        )

    let formatOne (n1: int) (strF: int -> string) =
        [|
            for i=1 to n1 do
                yield strF i
        |]
        |> String.concat "\t"

    let formatTwo (n1: int) (n2: int) (strF: int -> int -> string) =
        [|
            for i=1 to n1 do
                for j=1 to n2 do
                    yield strF i j
        |]
        |> String.concat "\t"

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
            formatOne expCount (sprintf "protein_abundance_assay[%i]")
        let protAbundanceStudVar =
            formatOne stVarCount (sprintf "protein_abundance_study_variable[%i]")
        let protAbundanceStDevStudVar =
            formatOne stVarCount (sprintf "protein_abundance_stdev_study_variable[%i]")
        let protAbundanceStdErrStudVar =
            formatOne stVarCount (sprintf "protein_abundance_std_error_study_variable[%i]")
        let sb = new Text.StringBuilder()
        sb.AppendFormat(
            "PRH\taccession\tdescription\ttaxid\tspecies\tdatabase\tdatabase_version\tsearch_engine\t{0}\t{1}\treliability\t{2}\t{3}\t{4}\tambiguity_members\tmodifications\turi\tgo_terms\tprotein_coverage\t{5}\t{6}\t{7}\t{8}",
            bestSearchEngineScore,
            searchEngineScoreMS,
            psmsMSRun,
            pepsDistMSRun,
            pepsUniqueMSRun,
            protAbundanceAssay,
            protAbundanceStudVar,
            protAbundanceStDevStudVar,
            protAbundanceStdErrStudVar
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
        let protAbundanceAssay =
            formatOne expCount (sprintf "peptide_abundance_assay[%i]")
        let protAbundanceStudVar =
            formatOne stVarCount (sprintf "peptide_abundance_study_variable[%i]")
        let protAbundanceStDevStudVar =
            formatOne stVarCount (sprintf "peptide_abundance_stdev_study_variable[%i]")
        let protAbundanceStdErrStudVar =
            formatOne stVarCount (sprintf "peptide_abundance_std_error_study_variable[%i]")
        let sb = new Text.StringBuilder()
        sb.AppendFormat(
            "PEH\tsequence\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\t{0}\t{1}\treliability\tmodifications\tretention_time\tretention_time_window\tcharge\tmass_to_charge\turi\tspectra_ref\t\t{2}\t{3}\t{4}\t{5}",
            bestSearchEngineScore,
            searchEngineScoreMS,
            protAbundanceAssay,
            protAbundanceStudVar,
            protAbundanceStDevStudVar,
            protAbundanceStdErrStudVar
        ) |> ignore
        sb.AppendLine() |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let psmHeader path (mzTABParams: Domain.MzTABParams) =
        let searchEnginecount = mzTABParams.SearchEngineNamesPep.Length
        let searchEngineScoreMS =
            formatOne searchEnginecount (sprintf "search_engine_score[%i]")
        let sb = new Text.StringBuilder()
        sb.AppendFormat(
            "PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\t{0}\treliability\tmodifications\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\turi\tspectra_ref\tpre\tpost\tstart\tend",
            searchEngineScoreMS
        ) |> ignore
        sb.AppendLine() |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let protBody path (proteinS: ProteinSection[]) =
        let sb = new Text.StringBuilder()
        proteinS
        |> Array.map (fun prot ->
            sb.AppendFormat(
                "PRT\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\t{21}",
                prot.accession,
                prot.description,
                prot.taxid,
                prot.species,
                prot.database,
                prot.database_version,
                prot.search_engine,
                //needs to be adapted for different search engines
                prot.best_search_engine_score,
                prot.search_engine_score_ms_run
                |> Array.map (fun score ->
                    score
                    |> concatRuns "null"
                )
                |> String.concat "\t",
                prot.reliability,
                prot.num_psms_ms_run
                |> concatRuns "0",
                prot.num_peptides_distinct_ms_run
                |> concatRuns "0",
                prot.num_peptides_unique_ms_run
                |> concatRuns "0",
                prot.ambiguity_members
                |> String.concat ",",
                prot.modifications,
                prot.uri,
                prot.go_terms,
                prot.protein_coverage,
                prot.protein_abundance_assay
                |> concatRuns "null",
                prot.protein_abundance_study_variable
                |> concatRuns "null",
                prot.protein_abundance_stdev_study_variable
                |> concatRuns "null",
                prot.protein_abundance_std_error_study_variable
                |> concatRuns "null"
            ) |> ignore
            sb.AppendLine()
        ) |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let pepBody path (peptideS: PeptideSection[]) =
        let sb = new Text.StringBuilder()
        peptideS
        |> Array.map (fun pep ->
            printfn "%A" pep.accession
            pep.accession
            |> Array.iter (fun prot ->
                sb.AppendFormat(
                    "PEP\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}",
                    pep.sequence,
                    prot,
                    pep.unique,
                    pep.database,
                    pep.database_version,
                    pep.search_engine,
                    //needs to be adapted for different search engines
                    pep.best_search_engine_score,
                    pep.search_engine_score_ms_run
                    |> Array.map (fun score ->
                        score
                        |> concatRuns "null"
                    )
                    |> String.concat "\t",
                    pep.reliability,
                    pep.modifications,
                    pep.retention_time
                    |> fun (x,y) ->
                        match (x,y) with
                        |None, Some b -> sprintf "null|%f" b
                        |Some a, None -> sprintf "%f|null" a
                        |Some a, Some b -> sprintf "%f|%f" a b
                    ,
                    pep.retention_time_window,
                    pep.charge,
                    pep.mass_to_charge,
                    pep.uri,
                    pep.spectra_ref,
                    pep.peptide_abundance_assay
                    |> concatRuns "null",
                    pep.peptide_abundance_study_variable
                    |> concatRuns "null",
                    pep.peptide_abundance_stdev_study_variable
                    |> concatRuns "null",
                    pep.peptide_abundance_std_error_study_variable
                    |> concatRuns "null"
                ) |> ignore
                sb.AppendLine()
                |> ignore
            )
        ) |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let psmBody path (psmS: PSMSection[]) =
        let sb = new Text.StringBuilder()
        psmS
        |> Array.map (fun psm ->
            psm.accession
            |> Array.iter (fun prot ->
                sb.AppendFormat(
                    "PSM\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}",
                    psm.sequence,
                    psm.PSM_ID,
                    prot,
                    psm.unique,
                    psm.database,
                    psm.database_version,
                    psm.search_engine,
                    //needs to be adapted for different search engines
                    psm.search_engine_score
                    |> Array.map string
                    |> String.concat "\t",
                    psm.reliability,
                    psm.modifications,
                    psm.retention_time,
                    psm.charge,
                    psm.exp_mass_to_charge,
                    psm.calc_mass_to_charge,
                    psm.uri,
                    psm.spectra_ref,
                    psm.pre,
                    psm.post,
                    psm.start,
                    psm.ending
                ) |> ignore
                sb.AppendLine()
                |> ignore
            )
        ) |> ignore
        IO.File.AppendAllText(path, sb.ToString())