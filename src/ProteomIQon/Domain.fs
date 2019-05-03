namespace ProteomIQon

open BioFSharp
open BioFSharp.Mz

module Domain = 
            
    type MS1CentroidizationParams = 
        {
            NumberOfScales        : int
            YThreshold            : float
            MzTolerance           : float
            SNRS_Percentile       : float
            MinSNR                : float         
        }
    
    type MS2CentroidizationParams = 
        {
            /// Centroidization
            NumberOfScales          : int
            Centroid_MzTolerance    : float
            SNRS_Percentile         : float
            MinSNR                  : float 
            /// Padding
            MaximumPaddingPoints    : int option
            Padding_MzTolerance     : float
            WindowSize              : int
            SpacingPerc             : float 
        } 

    type CentroidizationParams =
        {
            Centroid                    : bool
            UseManufacturerCentroids    : bool
            StartRetentionTime          : float
            EndRetentionTime            : float
            MS1Centroidization          : MS1CentroidizationParams
            MS2CentroidizationParams    : MS2CentroidizationParams
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
