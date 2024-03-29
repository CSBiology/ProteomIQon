namespace ProteomIQon

open BioFSharp
open BioFSharp.Mz
open MzIO.Binary

module Domain = 

    open BioFSharp.Mz.SearchDB

    type PaddingParams =
         {
            MaximumPaddingPoints    : int option
            Padding_MzTolerance     : float
            WindowSize              : int
            SpacingPerc             : float 
         }

    type YThreshold = 
        | Fixed of float
        | MinSpectrumIntensity

    type WaveletPeakPickingParams = 
        {
            /// Centroidization
            NumberOfScales          : int
            YThreshold              : YThreshold
            Centroid_MzTolerance    : float
            SNRS_Percentile         : float
            MinSNR                  : float
            RefineMZ                : bool
            SumIntensities          : bool
            PaddingParams           : PaddingParams option 
        } 

    type CentroidizationMode =
        | Manufacturer
        | Wavelet of WaveletPeakPickingParams
        //| Binning of BinningParams
        
    type PeakPicking = 
        | ProfilePeaks
        | Centroid of CentroidizationMode 

    type PreprocessingParams =
        {
            Compress                    : BinaryDataCompressionType
            StartRetentionTime          : float option 
            EndRetentionTime            : float option 
            MS1PeakPicking              : PeakPicking
            MS2PeakPicking              : PeakPicking
        }

    type PeptideDBParams = 
        {
        Name                : string
        FastaHeaderToName   : string -> string
        Protease            : Digestion.Protease
        MinMissedCleavages  : int
        MaxMissedCleavages  : int
        MaxMass             : float
        MinPepLength        : int
        MaxPepLength        : int
        IsotopicMod         : SearchInfoIsotopic list 
        MassMode            : MassMode
        MassFunction        : IBioItem -> float  
        FixedMods           : SearchModification list            
        VariableMods        : SearchModification list
        VarModThreshold     : int
        }

    type NTerminalSeries = ((IBioItem -> float) -> AminoAcids.AminoAcid list -> PeakFamily<TaggedMass.TaggedMass> list)
    
    type CTerminalSeries = ((IBioItem -> float) -> AminoAcids.AminoAcid list -> PeakFamily<TaggedMass.TaggedMass> list)

    type AndromedaParams = {
        /// selects the minimum and maximum amount of peaks retained in a 100 Da window, all combinations are tested and the best result is kept.
        PMinPMax                : int*int
        MatchingIonTolerancePPM : float       
        }

    type PeptideSpectrumMatchingParams = 
        {
            // Charge Determination Params
            ChargeStateDeterminationParams  : ChargeState.ChargeDetermParams             
            // +/- ppm of ion m/z to obtain target peptides from SearchDB. 
            LookUpPPM                       : float
            nTerminalSeries                 : NTerminalSeries
            cTerminalSeries                 : CTerminalSeries
            AndromedaParams                 : AndromedaParams
            ///
        }

    type PepValueFittingMethod = 
        //| LinearLogit
        | IRLS


    type EstimationParams = {
        QValueThreshold                 : float
        PepValueThreshold               : float 
        MaxIterations                   : int
        MinimumIncreaseBetweenIterations: float
        PepValueFittingMethod           : PepValueFittingMethod
    }

    type ScoreCutoff = {
        SequestLike : float
        Andromeda   : float 
        }

    type Threshold = 
        | Estimate of EstimationParams
        | Fixed of ScoreCutoff

    type PSMStatisticsParams = 
        {
            Threshold                   : Threshold 
            FastaHeaderToName           : string -> string
            KeepTemporaryFiles          : bool
        }

    type WindowSize = 
        | Fixed of int
        | EstimateUsingAutoCorrelation of float
    
    type SecondDerivativeParams = 
        {
            MinSNR                       : float  
            PolynomOrder                 : int
            WindowSize                   : WindowSize
        }
    
    type WaveletParameters = FSharpStats'.Wavelet.Parameters 

    type XicProcessing = 
        | SecondDerivative of SecondDerivativeParams
        | Wavelet of WaveletParameters
        
    type Window = 
        | Fixed of float
        | Estimate 

    type XicExtraction = 
        {
            ScanTimeWindow               : float 
            MzWindow_Da                  : Window 
            XicProcessing                : XicProcessing
            TopKPSMs                     : int option
        }
       
    type BaseLineCorrection = 
        {
            MaxIterations                : int 
            Lambda                       : int 
            P                            : float 
        }

    
    type Labeling =
        | Labelshift
        | Unlabeled
        | N15Labeling
        | N15LabelingOnly

    type QuantificationParams = 
        {
            PerformLabeledQuantification : Labeling
            XicExtraction                : XicExtraction
            BaseLineCorrection           : BaseLineCorrection option
        }

    type LabelEfficiencyCalculationParams =
        {
            LowerBound: float
            UpperBound: float
        }

    type AlignmentBasedQuantificationParams = 
        {
            PerformLabeledQuantification : bool
            PerformLocalWarp             : bool
            XicExtraction                : XicExtraction
            BaseLineCorrection           : BaseLineCorrection option
        }

    type FDRMethod =
        | Conservative
        | MAYU
        | DecoyTargetRatio

    type QValueMethod =
        | Storey
        | LogisticRegression of FDRMethod
        | NoQValue

    type ProteinInferenceParams = 
        {
            TryGetProteinIdentifier : string -> string option
            Protein                : ProteinInference.IntegrationStrictness
            Peptide                : ProteinInference.PeptideUsageForQuantification
            GroupFiles             : bool
            GetQValue              : QValueMethod
        }

    type SpectralLibraryParams =
        {
            MatchingTolerancePPM: float
        }

    type FilterOnField =
        {
            FieldName  : string
            UpperBound : float option
            LowerBound : float option
        }

    module FilterOnField =

        let create fieldName upperBound lowerBound =
            {
                FieldName  = fieldName
                UpperBound = upperBound
                LowerBound = lowerBound
            }

    type EssentialFields =
        {
            Light       : string
            Heavy       : string option
            ProteinIDs  : string
            PepSequence : string
            PepSequences: string
        }

    module EssentialFields =

        let create light heavy proteinIDs pepSequence pepSequences =
            {
                Light       = light
                Heavy       = heavy
                ProteinIDs  = proteinIDs
                PepSequence = pepSequence
                PepSequences= pepSequences
            }

    type AggregationMethod =
        | Sum
        | Mean
        | Median

    type Transform =
        | Log10
        | Log2
        | Ln
        | NoTransform

    type StatisticalMeasurement =
        | SEM
        | StDev
        | CV

    type TableSortParams =
        {
            SeparatorIn                 : string
            SeparatorOut                : char
            EssentialFields             : EssentialFields
            QuantFieldsToFilterOn       : FilterOnField[]
            ProtFieldsToFilterOn        : FilterOnField[]
            QuantColumnsOfInterest      : string[]
            ProtColumnsOfInterest       : (string*string)[]
            DistinctPeptideCount        : bool
            StatisticalMeasurements     : (string*StatisticalMeasurement)[]
            AggregatorFunction          : AggregationMethod
            AggregatorFunctionIntensity : AggregationMethod
            AggregatorPepToProt         : AggregationMethod
            Tukey                       : (string*float*Transform) []
        }

    module LabeledProteinQuantification = 
        
        type LabeledTransforms = 
            {
                Light: float -> float
                Heavy: float -> float
                Ratio: float -> float
            }   
        
        type LabeledAggregations = 
            {
                Light: seq<float> -> float
                Heavy: seq<float> -> float
                Ratio: seq<float> -> float
            }   
        
        type LabeledSingleFilters = 
            {
                Light: seq<(float -> bool)>
                Heavy: seq<(float -> bool)>
                Ratio: seq<(float -> bool)>
            }   

        type LabeledGroupFilters = 
            {
                Light: seq<(seq<float> -> float -> bool)>
                Heavy: seq<(seq<float> -> float -> bool)>
                Ratio: seq<(seq<float> -> float -> bool)>
            }   

        type AggregationParams = 
            {
                LabeledTransform        : LabeledTransforms 
                LabeledSingleFilters    : LabeledSingleFilters 
                LabeledGroupFilters     : LabeledGroupFilters 
                LabeledAggregation      : LabeledAggregations
            }

    type LabeledQuantificationParams = 
        {
            Correlation_Light_Heavy_Threshold: float option
            Alignment_QValue: float option
            ModificationFilter : string -> bool 
            AggregateGlobalModificationsParams   : LabeledProteinQuantification.AggregationParams
            AggregatePeptideChargeStatesParams   : LabeledProteinQuantification.AggregationParams option
            AggregateModifiedPeptidesParams      : LabeledProteinQuantification.AggregationParams option
            AggregateToProteinGroupsParams       : LabeledProteinQuantification.AggregationParams
        }

        
    module LabelFreeQuantification = 
        
        type Transforms = 
            {
                Light: float -> float
            }   
        
        type Aggregations = 
            {
                Light: seq<float> -> float
            }   
        
        type SingleFilters = 
            {
                Light: seq<(float -> bool)>
            }   

        type GroupFilters = 
            {
                Light: seq<(seq<float> -> float -> bool)>
            }   

        type AggregationParams = 
            {
                Transform        : Transforms 
                SingleFilters    : SingleFilters 
                GroupFilters     : GroupFilters 
                Aggregation      : Aggregations
            }

    type LabelFreeQuantificationParams = 
        {
            Alignment_QValue: float option
            ModificationFilter : string -> bool 
            AggregatePeptideChargeStatesParams   : LabelFreeQuantification.AggregationParams option
            AggregateModifiedPeptidesParams      : LabelFreeQuantification.AggregationParams option
            AggregateToProteinGroupsParams       : LabelFreeQuantification.AggregationParams
        } 

    type AlignmentBasedQuantStatisticsParams =
        {
            PositiveQuantMzCutoff: float
            NegativeQuantMzCutoff: float
            PositiveQuantCutoff: float
            NegativeQuantCutoff: float
        } 
               
    type ConsensusAlignmentAlgorithm = 
        | FastTree
        | Spline 

    ///
    type ConsensusSpectralLibraryParams = {
        // InitialPeptideSelection
        BinningWindowWidth                          : float
        FractionOfMostAbundandIonsPerBin            : float
        MinFragmentCount                            : int
        MinFragmentLadderIdx                        : int
        MinPeptideLength                            : int
        // XicExtraction
        RtWindowWidth                               : float
        // Matching
        FragMatchingBinWidth                        : float
        FragMatchingBinOffset                       : float
        MS2ScanRange                                : float*float
        // Filtering
        MaxRatioMostAbundandVsSecondAbundandPeak    : float
        ConsensusAlignmentAlgorithm :ConsensusAlignmentAlgorithm
        }

    type SpectrumSelection =
        |First
        |All

    ///
    type SWATHAnalysisParams = {
        // InitialPeptideSelection
        BinningWindowWidth                          : float
        FractionOfMostAbundandIonsPerBin            : float
        MinFragmentCount                            : int
        MinFragmentLadderIdx                        : int
        MinPeptideLength                            : int
        // XicExtraction
        RtWindowWidth                               : float
        // Matching
        FragMatchingBinWidth                        : float
        FragMatchingBinOffset                       : float
        MS2ScanRange                                : float*float
        // Filtering
        MaxRatioMostAbundandVsSecondAbundandPeak    : float
        }


    type MetaDataSection =
        {
            mzTab_version                    : string
            mzTab_mode                       : string
            mzTab_type                       : string
            mzTab_ID                         : (string)option
            title                            : (string)option
            description                      : string
            sample_processing                : ((Ontologies.SampleProcessing[]*int)[])option
            instrument_name                  : ((string*int)[])option
            instrument_source                : ((string*int)[])option
            instrument_analyzer              : ((string[]*int)[])option
            instrument_detector              : ((string*int)[])option
            software                         : ((string*int)[])option
            software_setting                 : ((string[]*int)[])option
            protein_search_engine_score      : ((Ontologies.SearchEngineScore*string*int)[])option
            peptide_search_engine_score      : ((Ontologies.SearchEngineScore*string*int)[])option
            psm_search_engine_score          : ((Ontologies.SearchEngineScore*string*int)[])option
            false_discovery_rate             : (string[])option
            publication                      : (string[])option
            contact_name                     : ((string*int)[])option
            contact_affiliation              : ((string*int)[])option
            contact_email                    : ((string*int)[])option
            uri                              : ((string*int)[])option
            fixed_mod                        : (Ontologies.Modification*int)[]
            fixed_mod_site                   : ((string*int)[])option
            fixed_mod_position               : ((Ontologies.ModificationPosition*int)[])option
            variable_mod                     : (Ontologies.Modification*int)[]
            variable_mod_site                : ((string*int)[])option
            variable_mod_position            : ((Ontologies.ModificationPosition*int)[])option
            quantification_method            : (string)option
            protein_quantification_unit      : (string)option
            peptide_quantification_unit      : (string)option
            ms_run_format                    : ((Ontologies.FileFormats*int)[])option
            // if unknown, "null" must be reported
            ms_run_location                  : (string*int)[]
            ms_run_id_format                 : ((Ontologies.IDFormats*int)[])option
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
            assay_quantification_reagent     : ((Ontologies.Labeling*int)[])
            assay_quantification_mod         : (((Ontologies.Modification*int)[]*int)[])option
            assay_quantification_mod_site    : (((string*int)[]*int)[])option
            assay_quantification_mod_position: (((Ontologies.ModificationPosition*int)[]*int)[])option
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
            colunit_peptide                  : (string)option
            colunit_psm                      : (string)option
            colunit_small_molecule           : (string)option
        }

    type TableSortFieldNames =
        {
            Proteingroup: string
            Experiment: string
            DistinctPeptideCount: string
            Quant_Heavy: string option
            Quant_Light: string
            QValue: string
            StudySubject: string
            Subject_SEM: string
            Subject_StDev: string
        }

    type MzTABParams =
        {
            ExperimentNames      : (string*int)[]
            StudyVariables       : (string*int[]*int)[]
            SearchEngineNamesProt: (Ontologies.SearchEngineScore*string*int)[]
            SearchEngineNamesPep : (Ontologies.SearchEngineScore*string*int)[]
            SearchEngineNamesPSM : (Ontologies.SearchEngineScore*string*int)[]
            Labeled              : bool
            MetaData             : MetaDataSection
            FieldNames           : TableSortFieldNames
        }

    type MzliteToMzMLParams =
        {
            Compress                    : BinaryDataCompressionType
            StartRetentionTime          : float option
            EndRetentionTime            : float option
        }

    type MzMLtoMzLiteParams = PreprocessingParams