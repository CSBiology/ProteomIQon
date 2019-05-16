namespace ProteomIQon

open BioFSharp
open BioFSharp.Mz

module Domain = 

    open BioFSharp.Mz
    open System.Linq
    open BioFSharp.Mz.SignalDetection
                
    type PaddingParams =
         {
            /// Padding
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

    type SearchDbParams = BioFSharp.Mz.SearchDB.SearchDbParams

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

            ExpectedMinimalCharge   : int ///TODO: learn from Data
            ExpectedMaximumCharge   : int ///TODO: learn from Data
            Width                   : float
            /// RelativeToStartPeak
            MinIntensity            : float
            /// RelativeToPriorPeak
            DeltaMinIntensity       : float
            NrOfRndSpectra          : int
            
            // SearchParams
            Protease                : Digestion.Protease
            MinMissedCleavages      : int
            MaxMissedCleavages      : int
            MaxMass                 : float
            MinPepLength            : int
            MaxPepLength            : int
            // valid symbol name of isotopic label in label table i.e. #N15
            IsotopicMod             : SearchDB.SearchInfoIsotopic list 
            MassMode                : SearchDB.MassMode
            MassFunction            : IBioItem -> float  
            FixedMods               : SearchDB.SearchModification list            
            VariableMods            : SearchDB.SearchModification list
            VarModThreshold         : int  
            // +/- ppm of ion m/z to obtain target peptides from SearchDB. 
            LookUpPPM               : float
            // lowest m/z, highest m/z
            MS2ScanRange            : float*float
            nTerminalSeries         : NTerminalSeries
            cTerminalSeries         : CTerminalSeries
            Andromeda               : AndromedaParams
            ///
        }

    type PEPEParams = 
        {
            QValueThreshold             : float
            PepValueThreshold           : float
            ParseProteinID              : string -> string
        }

    type XicExtraction = 
        {
            ScanTimeWindow               : float 
            MzWindow_Da                  : float 
            MinSNR                       : float  
            PolynomOrder                 : int
            WindowSize                   : int
        }
       
    type BaseLineCorrection = 
        {
            MaxIterations                : int 
            Lambda                       : float 
            P                            : float 
        }

    type QuantificationParams = 
        {
            PerformLabeledQuantification : bool
            XicExtraction                : XicExtraction
            BaseLineCorrection           : BaseLineCorrection option
        }

    // Add by HLWeil
    //type ProteinInferenceParams = 
    //      ...
   