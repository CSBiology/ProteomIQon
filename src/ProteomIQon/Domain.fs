namespace ProteomIQon

open BioFSharp
open BioFSharp.Mz

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
            PaddingParams           : PaddingParams option 
        } 

    type CentroidizationMode =
        | Manufacturer
        | Wavelet of WaveletPeakPickingParams
        
    type PeakPicking = 
        | ProfilePeaks
        | Centroid of CentroidizationMode 

    type PreprocessingParams =
        {
            Compress                    : bool
            StartRetentionTime          : float option 
            EndRetentionTime            : float option 
            MS1PeakPicking              : PeakPicking
            MS2PeakPicking              : PeakPicking
        }

    type PeptideDBParams = 
        {
        Name                : string
        FastaPath           : string
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
            // lowest m/z, highest m/z
            MS2ScanRange                    : float*float
            nTerminalSeries                 : NTerminalSeries
            cTerminalSeries                 : CTerminalSeries
            AndromedaParams                 : AndromedaParams
            ///
        }

    type PSMStatisticsParams = 
        {
            QValueThreshold             : float
            PepValueThreshold           : float
            FastaHeaderToName           : string -> string
            KeepTemporaryFiles          : bool
        }

    type WindowSize = 
        | Fixed of int
        | Estimate 

    type XicExtraction = 
        {
            ScanTimeWindow               : float 
            MzWindow_Da                  : float 
            MinSNR                       : float  
            PolynomOrder                 : int
            WindowSize                   : WindowSize
        }
       
    type BaseLineCorrection = 
        {
            MaxIterations                : int 
            Lambda                       : int 
            P                            : float 
        }

    type QuantificationParams = 
        {
            PerformLabeledQuantification : bool
            XicExtraction                : XicExtraction
            BaseLineCorrection           : BaseLineCorrection option
        }


    type ProteinInferenceParams = 
        {
            ProteinIdentifierRegex : string
            Protein                : ProteinInference.IntegrationStrictness
            Peptide                : ProteinInference.PeptideUsageForQuantification
            GroupFiles             : bool
        }
   