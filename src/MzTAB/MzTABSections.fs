namespace ProteomIQon

open System
open System.Linq
open MzIO
open FSharp.Stats
open FSharpAux
open FSharp.Plotly
open ProteomIQon.Domain
open MzTABAux

module MzTABSections =

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

    let getSameAccessions (allAligned: AlignedComplete[]) =
        let proteinGroups =
            allAligned |> Array.map (fun x ->
                x.TableSort.Protein
                |> String.split ';'
                |> Array.sort
            )
        groupsWithSameProteinAccession proteinGroups

    let metaDataSection path (md: MetaDataSection) =
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
                let paramF (p: Ontologies.SampleProcessing[][]) =
                    p
                    |> Array.map (fun x ->
                        x
                        |> Array.map (fun y ->
                            y.toParam
                        )
                        |> String.concat "|"
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tsample_processing[%i]\t%s")
                pickF >> paramF >> strF
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
                let paramF (p: Ontologies.SearchEngineScore[]) =
                    p
                    |> Array.map (fun x ->
                        x.toParam
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tprotein_search_engine_score[%i]\t%s")
                pickF >> paramF >> strF
            md.protein_search_engine_score
            |> matchOption f sb
        let pepSearchEngineScore =
            let f =
                let pickF =
                    sortAndPick (fun (_,_,x)->x) (fun (y,_,_)->y)
                let paramF (p: Ontologies.SearchEngineScore[]) =
                    p
                    |> Array.map (fun x ->
                        x.toParam
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tpeptide_search_engine_score[%i]\t%s")
                pickF >> paramF >> strF
            md.peptide_search_engine_score
            |> matchOption f sb
        let psmSearchEngineScore =
            let f =
                let pickF =
                    sortAndPick (fun (_,_,x)->x) (fun (y,_,_)->y)
                let paramF (p: Ontologies.SearchEngineScore[]) =
                    p
                    |> Array.map (fun x ->
                        x.toParam
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tpsm_search_engine_score[%i]\t%s")
                pickF >> paramF >> strF
            md.psm_search_engine_score
            |> matchOption f sb
        let fdr =
            let f =
                (String.concat "|") >> sprintf "MTD\tfalse_discovery_rate\t%s"
            md.false_discovery_rate
            |> matchOption f sb
        let publication =
            let strF =
                formatOneMD (sprintf "MTD\tpublication[%i]\t%s")
            md.publication
            |> matchOption strF sb
        let contactName =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tcontact[%i]-name\t%s")
                pickF >> strF
            md.contact_name
            |> matchOption f sb
        let contactAffiliation =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tcontact[%i]-affiliation\t%s")
                pickF >> strF
            md.contact_affiliation
            |> matchOption f sb
        let contactEmail =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tcontact[%i]-email\t%s")
                pickF >> strF
            md.contact_email
            |> matchOption f sb
        let uri =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\turi[%i]\t%s")
                pickF >> strF
            md.uri
            |> matchOption f sb
        let fixedMod =
            let f =
                let pickF =
                    sortAndPick snd fst
                let paramF (p: Ontologies.Modification[]) =
                    p
                    |> Array.map (fun x ->
                        x.toParam
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tfixed_mod[%i]\t%s")
                pickF >> paramF >> strF
            md.fixed_mod
            |> fun x -> sb.AppendLine (f x)
        let fixedModSite =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tfixed_mod[%i]-site\t%s")
                pickF >> strF
            md.fixed_mod_site
            |> matchOption f sb
        let fixedModPosition =
            let f =
                let pickF =
                    sortAndPick snd fst
                let paramF (p: Ontologies.ModificationPosition[]) =
                    p
                    |> Array.map (fun x ->
                        x.toParam
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tfixed_mod[%i]-position\t%s")
                pickF >> paramF >> strF
            md.fixed_mod_position
            |> matchOption f sb
        let variableMod =
            let f =
                let pickF =
                    sortAndPick snd fst
                let paramF (p: Ontologies.Modification[]) =
                    p
                    |> Array.map (fun x ->
                        x.toParam
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tvariable_mod[%i]\t%s")
                pickF >> paramF >> strF
            md.variable_mod
            |> fun x -> sb.AppendLine (f x)
        let variableModSite =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tvariable_mod[%i]-site\t%s")
                pickF >> strF
            md.variable_mod_site
            |> matchOption f sb
        let variableModPosition =
            let f =
                let pickF =
                    sortAndPick snd fst
                let paramF (p: Ontologies.ModificationPosition[]) =
                    p
                    |> Array.map (fun x ->
                        x.toParam
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tvariable_mod[%i]-position\t%s")
                pickF >> paramF >> strF
            md.variable_mod_position
            |> matchOption f sb
        let quantificationMethod =
            md.quantification_method
            |> matchOption (sprintf "MTD\tquantification_method\t%s") sb
        let proteinQuantificationUnit =
            md.protein_quantification_unit
            |> matchOption (sprintf "MTD\tprotein_quantification_unit\t%s") sb
        let peptideQuantificationUnit =
            md.peptide_quantification_unit
            |> matchOption (sprintf "MTD\tpeptide_quantification_unit\t%s") sb
        let msRunFormat =
            let f =
                let pickF =
                    sortAndPick snd fst
                let paramF (p: Ontologies.FileFormats[]) =
                    p
                    |> Array.map (fun x ->
                        x.toParam
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tms_run[%i]-format\t%s")
                pickF >> paramF >> strF
            md.ms_run_format
            |> matchOption f sb
        let msRunLocation =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tms_run[%i]-location\t%s")
                pickF >> strF
            md.ms_run_location
            |> fun x -> sb.AppendLine (f x)
        let msRunIdFormat =
            let f =
                let pickF =
                    sortAndPick snd fst
                let paramF (p: Ontologies.IDFormats[]) =
                    p
                    |> Array.map (fun x ->
                        x.toParam
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tms_run[%i]-id_format\t%s")
                pickF >> paramF >> strF
            md.ms_run_id_format
            |> matchOption f sb
        let msRunFragmentationMethod =
            let f =
                let pickF =
                    sortAndPick snd fst
                let concatF (input: string[][])=
                    input
                    |> Array.map (fun x ->
                        x
                        |> String.concat "|"
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tms_run[%i]-fragmentation_method\t%s")
                pickF >> concatF >> strF
            md.ms_run_fragmentation_method
            |> matchOption f sb
        let msRunHash =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tms_run[%i]-hash\t%s")
                pickF >> strF
            md.ms_run_hash
            |> matchOption f sb
        let msRunHashMethod =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tms_run[%i]-hash_method\t%s")
                pickF >> strF
            md.ms_run_hash_method
            |> matchOption f sb
        let custom =
            let strF =
                    formatOneMD (sprintf "MTD\tcustom[%i]\t%s")
            md.custom
            |> matchOption strF sb
        let sampleSpecies =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatTwoMD (sprintf "MTD\tsample[%i]-species[%i]\t%s")
                pickF >> strF
            md.sample_species
            |> matchOption f sb
        let sampleTissue =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatTwoMD (sprintf "MTD\tsample[%i]-tissue[%i]\t%s")
                pickF >> strF
            md.sample_tissue
            |> matchOption f sb
        let sampleCellType =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatTwoMD (sprintf "MTD\tsample[%i]-cell_type[%i]\t%s")
                pickF >> strF
            md.sample_cell_type
            |> matchOption f sb
        let sampleDisease =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatTwoMD (sprintf "MTD\tsample[%i]-disease[%i]\t%s")
                pickF >> strF
            md.sample_disease
            |> matchOption f sb
        let sampleDescription =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tsample[%i]-description\t%s")
                pickF >> strF
            md.sample_description
            |> matchOption f sb
        let sampleCustom =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatTwoMD (sprintf "MTD\tsample[%i]-custom[%i]\t%s")
                pickF >> strF
            md.sample_custom
            |> matchOption f sb
        let assayQuantificationReagent =
            let f =
                let pickF =
                    sortAndPick snd fst
                let paramF (p: Ontologies.Labeling[]) =
                    p
                    |> Array.map (fun x ->
                        x.toParam
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tassay[%i]-quantification_reagent\t%s")
                pickF >> paramF >> strF
            md.assay_quantification_reagent
            |> f
            |> sb.AppendLine
        let assayQuantificationMod =
            let f =
                let pickF =
                    sortAndPick snd fst
                let pickF2 i =
                    i
                    |> Array.map (sortAndPick snd fst)
                let paramF (p: Ontologies.Modification[][]) =
                    p
                    |> Array.map (fun y -> 
                        y
                        |> Array.map (fun x ->
                            x.toParam
                        )
                    )
                let strF =
                    formatTwoMD (sprintf "MTD\tassay[%i]-quantification_mod[%i]\t%s")
                pickF >> pickF2 >> paramF >> strF
            md.assay_quantification_mod
            |> matchOption f sb
        let assayQuantificationModSite =
            let f =
                let pickF =
                    sortAndPick snd fst
                let pickF2 i =
                    i
                    |> Array.map (sortAndPick snd fst)
                let strF =
                    formatTwoMD (sprintf "MTD\tassay[%i]-quantification_mod[%i]-site\t%s")
                pickF >> pickF2 >> strF
            md.assay_quantification_mod_site
            |> matchOption f sb
        let assayQuantificationModPosition =
            let f =
                let pickF =
                    sortAndPick snd fst
                let pickF2 i =
                    i
                    |> Array.map (sortAndPick snd fst)
                let paramF (p: Ontologies.ModificationPosition[][]) =
                    p
                    |> Array.map (fun y -> 
                        y
                        |> Array.map (fun x ->
                            x.toParam
                        )
                    )
                let strF =
                    formatTwoMD (sprintf "MTD\tassay[%i]-quantification_mod[%i]-position\t%s")
                pickF >> pickF2 >> paramF >> strF
            md.assay_quantification_mod_position
            |> matchOption f sb
        let assaySampleRef =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tassay[%i]-sample_ref\t%s")
                pickF >> strF
            md.assay_sample_ref
            |> matchOption f sb
        let assayMsRunRef =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tassay[%i]-ms_run_ref\t%s")
                pickF >> strF
            md.assay_ms_run_ref
            |> matchOption f sb
        let studyVariableAssayRefs =
            let f =
                let pickF =
                    sortAndPick snd fst
                let concatF (input: int[][])=
                    input
                    |> Array.map (fun x ->
                        x
                        |> Array.map (fun y ->
                            sprintf "assay[%i]" y
                        )
                        |> String.concat ","
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tstudy_variable[%i]-assay_refs\t%s")
                pickF >> concatF >> strF
            md.study_variable_assay_refs
            |> matchOption f sb
        let studyVariableSampleRefs =
            let f =
                let pickF =
                    sortAndPick snd fst
                let concatF (input: int[][])=
                    input
                    |> Array.map (fun x ->
                        x
                        |> Array.map (fun y ->
                            sprintf "assay[%i]" y
                        )
                        |> String.concat ","
                    )
                let strF =
                    formatOneMD (sprintf "MTD\tstudy_variable[%i]-sample_refs\t%s")
                pickF >> concatF >> strF
            md.study_variable_sample_refs
            |> matchOption f sb
        let studyVariableDescription =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tstudy_variable[%i]-description\t%s")
                pickF >> strF
            md.study_variable_description
            |> matchOption f sb
        let cvLabel =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tcv[%i]-label\t%s")
                pickF >> strF
            md.cv_label
            |> matchOption f sb
        let cvFullName =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tcv[%i]-full_name\t%s")
                pickF >> strF
            md.cv_full_name
            |> matchOption f sb
        let cvVersion =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tcv[%i]-version\t%s")
                pickF >> strF
            md.cv_version
            |> matchOption f sb
        let cvUrl =
            let f =
                let pickF =
                    sortAndPick snd fst
                let strF =
                    formatOneMD (sprintf "MTD\tcv[%i]-url\t%s")
                pickF >> strF
            md.cv_url
            |> matchOption f sb
        let colunitProtein =
            md.colunit_protein
            |> matchOption (sprintf "MTD\tcolunit-protein\t%s") sb
        let colunitPeptide =
            md.colunit_peptide
            |> matchOption (sprintf "MTD\tcolunit-peptide\t%s") sb
        let colunitPSM =
            md.colunit_psm
            |> matchOption (sprintf "MTD\tcolunit-psm\t%s") sb
        IO.File.AppendAllText(path, sb.ToString())

    let proteinSection (allAligned: AlignedComplete[]) (mzTABParams: Domain.MzTABParams) =
        let experimentNames = mzTABParams.ExperimentNames
        let groupedTab =
            allAligned
            |> Array.groupBy (fun x -> x.TableSort.Protein)
            |> Array.map snd
            //|> Array.map (fun x ->
            //    x
            //    |> Array.sortBy (fun y -> y.TableSort.Experiment)
            //)
            |> Array.map (fun x ->
                x
                |> Array.groupBy (fun y -> y.TableSort)
            )
        let initRatio names tablesort =
            findValueNumberedProt names tablesort "Ratio"
        let initStDev names tablesort =
            findValueNumberedProt names tablesort "Ratio_StDev"
        let initSEM names tablesort =
            findValueNumberedProt names tablesort "Ratio_SEM"
        groupedTab
        |> Array.map (fun protGroup ->
            let studyVars =
                let quant =
                    initRatio experimentNames (protGroup |> Array.map fst)
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
            let studyVarsStDev =
                let quant =
                    initStDev experimentNames (protGroup |> Array.map fst)
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
            let studyVarsStdErr =
                let quant =
                    initSEM experimentNames (protGroup |> Array.map fst)
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
            let proteinGroup =
                (fst protGroup.[0]).Protein
                |> String.split ';'
                |> Array.sort
            {
                accession                                 = proteinGroup
                description                               = "null"
                taxid                                     = 3055
                species                                   =
                    let experimentsIdentifeidNumber =
                        protGroup
                        |> Array.map (fun (ts,_) ->
                            experimentNames
                            |> Array.find (fun (name,_) -> name = ts.Experiment)
                            |> fun (name,number) -> number
                        )
                    match mzTABParams.MetaData.sample_species with
                    |None -> "null"
                    |Some species -> 
                        let valSpecies =
                            species
                            |> Array.filter (fun (_,expN) -> experimentsIdentifeidNumber |> Array.contains expN)
                        valSpecies
                        |> Array.collect fst
                        |> Array.distinct
                        |> String.concat "|"
                database                                  = "Chlamy_JGI5_5(Cp_Mp).fasta"
                database_version                          = "2016_09"
                search_engine                             =
                    mzTABParams.MetaData.protein_search_engine_score.Value
                    |> Array.map (fun (x,_,_) -> x.toParam)
                    |> String.concat "|"
                best_search_engine_score                  =
                    protGroup
                    |> Array.minBy (fun (prot,peps) -> prot.QValue)
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
                ambiguity_members                         = 
                    if proteinGroup.Length = 1 then
                        [|"null"|]
                    else
                        proteinGroup.[1..]
                modifications                             = "null"
                // have to look how to handle uri with protein identifications from different sources
                uri                                       = "null"
                go_terms                                  = [|"null"|]
                // remember sequence
                protein_coverage                          = 0.
                protein_abundance_assay                   =
                    if mzTABParams.Labeled then
                        let light =
                            findValueNumberedProt experimentNames (protGroup |> Array.map fst) "Quant_Light"
                            |> Array.sortBy fst
                        let heavy =
                            findValueNumberedProt experimentNames (protGroup |> Array.map fst) "Quant_Heavy"
                            |> Array.sortBy fst
                        Array.map2 (fun light heavy -> [|light; heavy|]) light heavy
                        |> Array.concat
                    
                    else
                        findValueNumberedProt experimentNames (protGroup |> Array.map fst) "Quant_Light"
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
                    studyVarsStDev
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
                protein_abundance_std_error_study_variable=
                    studyVarsStdErr
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
                    
            }
        )

    let peptideSection (allAligned: AlignedComplete[]) (mzTABParams: Domain.MzTABParams): PeptideSection[] =
        let experimentNames = mzTABParams.ExperimentNames
        let groupedTab =
            allAligned
            // all same peptide ions are put into an array
            |> Array.groupBy (fun pep -> pep.Quant.StringSequence, pep.Quant.Charge, pep.Quant.GlobalMod)
            |> Array.map snd
            // sort peptides by experiment
            //|> Array.map (fun x ->
            //    x
            //    |> Array.sortBy (fun y -> y.TableSort.Experiment)
            //)
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
                    |> Array.map (fun x -> 
                        x.TableSort.Protein
                        |> String.split ';'
                    )
                    |> Array.sort
                )
                |> Array.distinct
            let ratio =
                let light =
                    findValueNumberedPep experimentNames forF "Quant_Light"
                    |> Array.sortBy fst
                let heavy =
                    findValueNumberedPep experimentNames forF "Quant_Heavy"
                    |> Array.sortBy fst
                Array.map2 (fun (i,l) (j,h) ->
                    if i <> j then failwith "Experimental data for peptide ratios missing"
                    match l,h with
                    | Some x, Some y -> i, Some (x/y)
                    | _, None -> i, None
                    | None, _ -> i, None
                ) light heavy
            let studyVars =
                mzTABParams.StudyVariables
                |> Array.map (fun (name,assays,number) ->
                    let corrQuant =
                        ratio
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
                    "null"
                database_version                          =
                    "null"
                search_engine                             =
                    mzTABParams.MetaData.peptide_search_engine_score.Value
                    |> Array.map (fun (x,_,_) -> x.toParam)
                    |> String.concat "|"
                best_search_engine_score                  =
                    mzTABParams.SearchEngineNamesPep
                    |> Array.map (fun (searchengine,fieldName,number) ->
                        findValueNumberedPep experimentNames forF fieldName
                        |> Array.sortBy fst
                        |>Array.choose (fun (i,x) ->
                            match x with
                            | None -> None
                            | Some v -> Some (i, v)
                        )
                        // need to add an identifier for search engines whether a high or low score is better
                        |> fun x ->
                            if fieldName = "QValue" then
                                x |> Array.minBy snd
                            else
                                x |> Array.maxBy snd
                            
                    )
                search_engine_score_ms_run                =
                    mzTABParams.SearchEngineNamesPep
                    |> Array.map (fun (searchengine,fieldName,number) ->
                        findValueNumberedPep experimentNames forF fieldName
                        |> Array.sortBy fst
                    )
                reliability                               =
                    3
                modifications                             =
                    "null"
                retention_time                            =
                    pepGroup
                    |> Array.sortBy (fun x -> (fst x).MeanPercolatorScore)
                    |> Array.head
                    |> fun (quant,rest) -> MzTABAux.getTargetScanTime quant
                    |> fun x -> x * 60.
                retention_time_window                     =
                    1.,2.
                charge                                    =
                    (fst pepGroup.[0]).Charge
                mass_to_charge                            =
                    pepGroup
                    |> Array.maxBy (fun x -> (fst x).MeanPercolatorScore)
                    |> fun (peptide,rest) -> peptide.PrecursorMZ
                uri                                       =
                    "null"
                spectra_ref                               =
                    "null"
                peptide_abundance_assay                   =
                    if mzTABParams.Labeled then
                        let light =
                            findValueNumberedPep experimentNames forF "Quant_Light"
                            |> Array.sortBy fst
                        let heavy =
                            findValueNumberedPep experimentNames forF "Quant_Heavy"
                            |> Array.sortBy fst
                        Array.map2 (fun light heavy -> [|light; heavy|]) light heavy
                        |> Array.concat
                    else
                        findValueNumberedPep experimentNames forF "Quant_Light"
                        |> Array.sortBy fst
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
                labeling =
                    let gMod =
                        pepGroup
                        |> Array.map (fun (x,_) ->
                            x.GlobalMod
                        )
                        |> Array.distinct
                    if gMod.Length <> 1 then failwith "gmod"
                    match gMod.[0] with
                    | 0 -> Ontologies.Labeling.N14
                    | 1 -> Ontologies.Labeling.N15
                    | _ -> failwith "Unexpected GlobalMod. GlobalMod must be either 0 or 1"
            }
        )

    let psmSection (allAligned: AlignedComplete[]) (mzTABParams: Domain.MzTABParams): PSMSection[] =
        let groupedTab =
            allAligned
            |> Array.groupBy (fun x -> x.Qpsm)
        groupedTab
            |> Array.map (fun (psm,rest) ->
                if rest.Length <> 1 then failwith "Unexpected number of entries for a single psm"
                let corrProt =
                    rest.[0].TableSort.Protein
                    |> String.split ';'
                    |> Array.sort
                {
                    sequence                                  =
                        psm.StringSequence
                    PSM_ID = psm.ScanNr
                    accession                                 =
                        corrProt
                    unique                                    =
                        // with the current check for a rest length of 1, only unique psms are possible
                        // as far as i know, the proteininference assigns peptides to one protein only
                        1
                        //if corrProt.Length <> 1 then
                        //    0
                        //else
                        //    1
                    database                                  =
                        "null"
                    database_version                          =
                        "null"
                    search_engine                             =
                        mzTABParams.MetaData.psm_search_engine_score.Value
                        |> Array.map (fun (x,_,_) -> x.toParam)
                        |> String.concat "|"
                    search_engine_score                       =
                        mzTABParams.SearchEngineNamesPSM
                        |> Array.map (fun (searchengine,fieldName,number) ->
                            fieldFuncPSM psm fieldName
                        )
                    reliability                               =
                        3
                    modifications                             =
                        match psm.GlobalMod with
                        | 0 -> Ontologies.Labeling.N14.toParam
                        | 1 -> Ontologies.Labeling.N15.toParam
                        | _ -> failwith "Unexpected GlobalMod. GlobalMod must be either 0 or 1"
                    retention_time                            =
                        psm.ScanTime
                        |> fun x -> x * 60.
                    charge                                    =
                        psm.Charge
                    exp_mass_to_charge                            =
                        psm.PrecursorMZ
                    calc_mass_to_charge                            =
                        BioFSharp.Mass.toMZ psm.TheoMass (float psm.Charge)
                    uri                                       =
                        "null"
                    spectra_ref                               =
                        let exp =
                            rest 
                            |> Array.map (fun x ->
                                x.TableSort.Experiment
                            )
                            |> Array.distinct
                        if exp.Length <> 1 then failwith "unexpected number of runs for psm"
                        let expNum =
                            mzTABParams.ExperimentNames
                            |> Array.find (fun (name,number) -> name=exp.[0])
                            |> snd

                        sprintf ("ms_run[%i]:%s") expNum psm.PSMId
                    pre="null"
                    post="null"
                    start=1
                    ending=1
                }
        )