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
            Protein: string
            [<FieldAttribute("\"Key #2\"")>]
            Experiment: string
            DistinctPeptideCount: float
            Quant_Heavy: float
            Quant_Light: float
            QValue: float
            Ratio: float
            Ratio_CV: float
            Ratio_SEM: float
            Ratio_StDev: float
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

    type PSMStatisticsResult = {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_
        PSMId                        : string
        GlobalMod                    : int
        PepSequenceID                : int
        ModSequenceID                : int
        Label                        : int
        // ascending ms2 id (file specific)
        ScanNr                       : int
        ScanTime                     : float
        Charge                       : int
        PrecursorMZ                  : float
        TheoMass                     : float
        AbsDeltaMass                 : float
        PeptideLength                : int
        MissCleavages                : int
        SequestScore                 : float
        SequestNormDeltaBestToRest   : float
        SequestNormDeltaNext         : float
        AndroScore                   : float
        AndroNormDeltaBestToRest     : float
        AndroNormDeltaNext           : float
        XtandemScore                 : float
        XtandemNormDeltaBestToRest   : float
        XtandemNormDeltaNext         : float
        PercolatorScore              : float
        QValue                       : float
        PEPValue                     : float
        StringSequence               : string
        ProteinNames                 : string
        }

    type QuantificationResult = {
        StringSequence                              : string
        GlobalMod                                   : int
        Charge                                      : int
        PepSequenceID                               : int
        ModSequenceID                               : int
        PrecursorMZ                                 : float
        MeasuredMass                                : float
        TheoMass                                    : float
        AbsDeltaMass                                : float
        MeanPercolatorScore                         : float
        QValue                                      : float
        PEPValue                                    : float
        ProteinNames                                : string
        QuantMz_Light                               : float
        Quant_Light                                 : float
        MeasuredApex_Light                          : float
        Seo_Light                                   : float
        Params_Light                                : string
        Difference_SearchRT_FittedRT_Light          : float
        KLDiv_Observed_Theoretical_Light            : float
        KLDiv_CorrectedObserved_Theoretical_Light   : float
        QuantMz_Heavy                               : float
        Quant_Heavy                                 : float
        MeasuredApex_Heavy                          : float
        Seo_Heavy                                   : float
        Params_Heavy                                : string
        Difference_SearchRT_FittedRT_Heavy          : float
        KLDiv_Observed_Theoretical_Heavy            : float
        KLDiv_CorrectedObserved_Theoretical_Heavy   : float
        Correlation_Light_Heavy                     : float
        QuantificationSource                        : string
        IsotopicPatternMz_Light                     : string
        IsotopicPatternIntensity_Observed_Light     : string
        IsotopicPatternIntensity_Corrected_Light    : string
        RtTrace_Light                               : string
        IntensityTrace_Observed_Light               : string
        IntensityTrace_Corrected_Light              : string
        IsotopicPatternMz_Heavy                     : string
        IsotopicPatternIntensity_Observed_Heavy     : string
        IsotopicPatternIntensity_Corrected_Heavy    : string
        RtTrace_Heavy                               : string
        IntensityTrace_Observed_Heavy               : string
        IntensityTrace_Corrected_Heavy              : string
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
            best_search_engine_score                  : (int*float) []
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
            labeling                                  : Ontologies.Labeling
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
        SeqIO.Seq.fromFileWithCsvSchema<PSMStatisticsResult>(path, '\t', true, schemaMode=SchemaReader.Csv.SchemaModes.Fill)
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
            let proteinGroup =
                (fst protGroup.[0]).Protein
                |> String.split ';'
            {
                accession                                 = proteinGroup |> Array.head
                description                               = ""
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
                database                                  = "Chlamy.db"
                database_version                          = "19-Apr-20 21:44"
                search_engine                             =
                    mzTABParams.MetaData.protein_search_engine_score.Value
                    |> Array.map (fun (x,_,_) -> x.toParam)
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
                ambiguity_members                         = 
                    if proteinGroup.Length = 1 then
                        [|"null"|]
                    else
                        proteinGroup.[1..]
                modifications                             = "null"
                // have to look how to handle uri with protein identifications from different sources
                uri                                       = "null"
                go_terms                                  = [||]
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
            // all same peptide ions are put into an array
            |> Array.groupBy (fun pep -> pep.Quant.StringSequence, pep.Quant.Charge, pep.Quant.GlobalMod)
            |> Array.map snd
            // sort peptides by experiment
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
                // must be adapted to search engine entries in metadata
                    let percolator =
                        pepGroup
                        |> Array.maxBy (fun (peptide,rest) -> peptide.MeanPercolatorScore)
                        |> fun (peptide,rest) -> peptide.MeanPercolatorScore
                    let qVal =
                        pepGroup
                        |> Array.minBy (fun (peptide,rest) -> peptide.QValue)
                        |> fun (peptide,rest) -> peptide.QValue
                    [|1,percolator; 2,qVal|]
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
            if mzTABParams.Labeled then
                formatOne (expCount*2) (sprintf "protein_abundance_assay[%i]")
            else
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
        let pepAbundanceAssay =
            if mzTABParams.Labeled then
                formatOne (expCount*2) (sprintf "peptide_abundance_assay[%i]")
            else
                formatOne expCount (sprintf "peptide_abundance_assay[%i]")
        let pepAbundanceStudVar =
            formatOne stVarCount (sprintf "peptide_abundance_study_variable[%i]")
        let pepAbundanceStDevStudVar =
            formatOne stVarCount (sprintf "peptide_abundance_stdev_study_variable[%i]")
        let pepAbundanceStdErrStudVar =
            formatOne stVarCount (sprintf "peptide_abundance_std_error_study_variable[%i]")
        let sb = new Text.StringBuilder()
        sb.AppendFormat(
            "PEH\tsequence\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\t{0}\t{1}\treliability\tmodifications\tretention_time\tretention_time_window\tcharge\tmass_to_charge\turi\tspectra_ref\t{2}\t{3}\t{4}\t{5}\topt_global_global_mod",
            bestSearchEngineScore,
            searchEngineScoreMS,
            pepAbundanceAssay,
            pepAbundanceStudVar,
            pepAbundanceStDevStudVar,
            pepAbundanceStdErrStudVar
        ) |> ignore
        sb.AppendLine() |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let psmHeader path (mzTABParams: Domain.MzTABParams) =
        let searchEnginecount = mzTABParams.SearchEngineNamesPSM.Length
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
                //prot.reliability,
                "null",
                prot.num_psms_ms_run
                |> concatRuns "0",
                prot.num_peptides_distinct_ms_run
                |> concatRuns "0",
                //prot.num_peptides_unique_ms_run
                //|> concatRuns "0",
                "null",
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
            pep.accession
            |> Array.iter (fun prot ->
                sb.AppendFormat(
                    "PEP\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}",
                    pep.sequence,
                    prot,
                    pep.unique,
                    pep.database,
                    pep.database_version,
                    pep.search_engine,
                    //needs to be adapted for different search engines
                    pep.best_search_engine_score
                    |> Array.sortBy fst
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
                    |> fun (x,y) ->
                        match (x,y) with
                        |None, Some b -> sprintf "null|%f" b
                        |Some a, None -> sprintf "%f|null" a
                        |Some a, Some b -> sprintf "%f|%f" a b
                    ,
                    //pep.retention_time_window,
                    "null",
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
                    |> concatRuns "null",
                    pep.labeling.toParam
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
                    //psm.start,
                    //psm.ending
                    "null",
                    "null"
                ) |> ignore
                sb.AppendLine()
                |> ignore
            )
        ) |> ignore
        IO.File.AppendAllText(path, sb.ToString())

    let createMzTab path tab prot quant qpsm (param: Domain.MzTABParams) =
        let protFiles =
            IO.Directory.GetFiles(prot,("*.prot"))
        let quantFiles =
            IO.Directory.GetFiles(quant,("*.quant"))
        let qpsmFiles =
            IO.Directory.GetFiles(qpsm,("*.qpsm"))
        let allAligned =
            allignAllFiles tab protFiles quantFiles qpsmFiles
        let protSection =
            proteinSection allAligned param
        let pepSection =
            peptideSection allAligned param
        let psmSection' =
            psmSection allAligned param
        metaDataSection path param.MetaData
        protHeader path param
        protBody path protSection
        pepHeader path param
        pepBody path pepSection
        psmHeader path param
        psmBody path psmSection'