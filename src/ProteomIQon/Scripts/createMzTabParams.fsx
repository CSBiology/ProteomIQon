#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"

open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain
open ProteomIQon.Ontologies

let defaultMzTabParams :Dto.MzTABParams = 
    let metaData =
        {
            mzTab_version                     = "1.0.0"
            mzTab_mode                        = "Complete"
            mzTab_type                        = "Quantification"
            mzTab_ID                          = None
            title                             = None
            description                       = "Your experiment description belongs here"
            sample_processing                 =
                Some
                    [|
                        ([|SampleProcessing.SDSPage; SampleProcessing.CoomassieStaining; SampleProcessing.ExcisedBand|],1);
                        ([|SampleProcessing.EnzymeDigestion; SampleProcessing.Trypsin; SampleProcessing.LysC|],2);
                        ([|SampleProcessing.SampleDesalting; SampleProcessing.DesaltingMicrocolumn|],3)
                        ([|SampleProcessing.PeakPicking|],4)
                    |]
            instrument_name                   = None
            instrument_source                 = None
            instrument_analyzer               = None
            instrument_detector               = None
            software                          = 
                Some
                    [|("ProteomIQon",1)|]
            software_setting                  = None
            protein_search_engine_score       = 
                Some
                    [|(SearchEngineScore.ProteinQValue, "QValue", 1)|]
            peptide_search_engine_score       = 
                Some
                    [|
                        (SearchEngineScore.PeptideQValue, "QValue", 1);
                        (SearchEngineScore.Percolator, "MeanPercolatorScore", 2)
                    |]
            psm_search_engine_score           = 
                Some
                    [|
                        (SearchEngineScore.Percolator, "PercolatorScore", 1);
                        (SearchEngineScore.SequestConsensusScore, "SequestScore", 2);
                        (SearchEngineScore.Andromeda, "AndroScore", 3);
                        (SearchEngineScore.XTandem, "XtandemScore", 4)
                    |]
            false_discovery_rate              = None
            publication                       = None
            contact_name                      = None
            contact_affiliation               = None
            contact_email                     = None
            uri                               = None
            fixed_mod                         = 
                [|(Modification.NoFixedModSearched, 1)|]
            fixed_mod_site                    = None
            fixed_mod_position                = None
            variable_mod                      =
                [|(Modification.NoVariableModsSearched, 1)|]
            variable_mod_site                 = None
            variable_mod_position             = None
            quantification_method             = None
            protein_quantification_unit       = 
                Some
                    "[PRIDE, PRIDE:0000395, Ratio, ]"
            peptide_quantification_unit       = 
                Some
                    "[PRIDE, PRIDE:0000395, Ratio, ]"
            ms_run_format                     =
                Some
                    [|(FileFormats.WIFF, 1)|]
            ms_run_location                   =
                [|("File1", 1)|]
            ms_run_id_format                  =
                Some
                    [|(IDFormats.WIFF, 1)|]
            ms_run_fragmentation_method       = None
            ms_run_hash                       = None
            ms_run_hash_method                = None
            custom                            = None
            sample_species                    = None
            sample_tissue                     = None
            sample_cell_type                  = None
            sample_disease                    = None
            sample_description                = None
            sample_custom                     = None
            assay_quantification_reagent      = None
            assay_quantification_mod          = None
            assay_quantification_mod_site     = None
            assay_quantification_mod_position = None
            assay_sample_ref                  = None
            assay_ms_run_ref                  = None
            study_variable_assay_refs         = None
            study_variable_sample_refs        = None
            study_variable_description        = None
            cv_label                          = None
            cv_full_name                      = None
            cv_version                        = None
            cv_url                            = None
            colunit_protein                   = None
            colunit_peptide                   = None
            colunit_psm                       = None
            colunit_small_molecule            = None
        }
    {
        ExperimentNames =
            [|"20170517 TM FScon3001",1; "20170517 TM FScon3003",2; "20170517 TM FScon3005",3;
            "20170517 TM FScon3007",4; "20170517 TM FScon3009",5; "20170517 TM FScon3011",6;
            "20170517 TM FScon3013",7; "20170517 TM FScon3015",8; "20170517 TM FScon3017",9|]
        StudyVariables = [|"Ratios",[|1 .. 9|],1|]
        SearchEngineNamesProt = metaData.protein_search_engine_score.Value
        SearchEngineNamesPep = metaData.peptide_search_engine_score.Value
        SearchEngineNamesPSM = metaData.psm_search_engine_score.Value
        MetaData = metaData
    }
    


let serialized = 
    defaultMzTabParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\MzTabParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\MzTabParams.json")
    |> Json.deserialize<Dto.MzTABParams>
    |> MzTABParams.toDomain
