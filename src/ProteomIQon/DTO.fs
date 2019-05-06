namespace ProteomIQon

open BioFSharp
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open Domain

[<AutoOpen>]
module Common =
       
    module MassMode =
   
        let toDomain (massMode: MassMode) =
            match massMode with
            | MassMode.Monoisotopic -> BioItem.initMonoisoMassWithMemP
            | MassMode.Average      -> BioItem.initAverageMassWithMemP


    type Protease =
        | Trypsin

    module Protease = 

        let toDomain protease = 
            match protease with
            | Trypsin -> Digestion.Table.getProteaseBy "Trypsin"

 
    type Modification =
        | Acetylation'ProtNTerm'
        | Carbamidomethyl'Cys'
        | Oxidation'Met'
        | Phosphorylation'Ser'Thr'Tyr'
    
    module Modification  =

        let toDomain modification = 
            match modification with
            | Acetylation'ProtNTerm'        -> SearchDB.Table.acetylation'ProtNTerm'
            | Carbamidomethyl'Cys'          -> SearchDB.Table.carbamidomethyl'Cys'
            | Oxidation'Met'                -> SearchDB.Table.oxidation'Met'
            | Phosphorylation'Ser'Thr'Tyr'  -> SearchDB.Table.phosphorylation'Ser'Thr'Tyr'


    type IsotopicMod =
        | N15

    module IsotopicMod =
        let toDomain isoMod = 
            match isoMod with
            | N15 -> (SearchDB.createSearchInfoIsotopic "N15" Elements.Table.N Elements.Table.Heavy.N15)
 
    type NTerminalSeries = 
        | A
        | B
        | C
        | AB
        | AC
        | BC
        | ABC
        
    module NTerminalSeries = 
        let toDomain nTermSeries =
            match nTermSeries with
            | A   -> Fragmentation.Series.aOfBioList
            | B   -> Fragmentation.Series.bOfBioList
            | C   -> Fragmentation.Series.cOfBioList
            | AB  -> Fragmentation.Series.abOfBioList
            | AC  -> Fragmentation.Series.acOfBioList
            | BC  -> Fragmentation.Series.bcOfBioList
            | ABC -> Fragmentation.Series.abcOfBioList


    type CTerminalSeries = 
        | X
        | Y
        | Z
        | XY
        | XZ
        | YZ
        | XYZ
        
    module CTerminalSeries = 
        let toDomain nTermSeries =
            match nTermSeries with
            | X   -> Fragmentation.Series.xOfBioList
            | Y   -> Fragmentation.Series.yOfBioList
            | Z   -> Fragmentation.Series.zOfBioList
            | XY  -> Fragmentation.Series.xyOfBioList
            | XZ  -> Fragmentation.Series.xzOfBioList
            | YZ  -> Fragmentation.Series.yzOfBioList
            | XYZ -> Fragmentation.Series.xyzOfBioList

    let parseProteinIdUsing regex =
        match regex with
        | "ID" | "id" | "Id" | "" -> 
            id
        | pattern ->        
            (fun (inp : string)  -> System.Text.RegularExpressions.Regex.Match(inp,pattern).Value)


///
module Dto = 

    type CentroidizationParams =
        {
            Centroid                    : bool
            UseManufacturerCentroids    : bool
            StartRetentionTime          : float
            EndRetentionTime            : float
            MS1Centroidization          : Domain.MS1CentroidizationParams
            MS2CentroidizationParams    : Domain.MS2CentroidizationParams
        }


    module CentroidizationParams = 

        let toDomain (dtoCentroidizationParams: CentroidizationParams ) = 
                {
                    Centroid                    = dtoCentroidizationParams.Centroid                
                    UseManufacturerCentroids    = dtoCentroidizationParams.UseManufacturerCentroids
                    StartRetentionTime          = dtoCentroidizationParams.StartRetentionTime      
                    EndRetentionTime            = dtoCentroidizationParams.EndRetentionTime        
                    MS1Centroidization          = dtoCentroidizationParams.MS1Centroidization      
                    MS2CentroidizationParams    = dtoCentroidizationParams.MS2CentroidizationParams
                }

    type SearchDbParams = 
        {
        // name of database i.e. Creinhardtii_236_protein_full_labeled
        Name                        : string
        // path of db storage folder
        DbFolder                    : string
        FastaPath                   : string
        ParseProteinIDRegexPattern  : string 
        Protease                    : Protease
        MinMissedCleavages          : int
        MaxMissedCleavages          : int
        MaxMass                     : float
        MinPepLength                : int
        MaxPepLength                : int
        // valid symbol name of isotopic label in label table i.e. #N15
        IsotopicMod                 : IsotopicMod list 
        MassMode                    : MassMode
        FixedMods                   : Modification list            
        VariableMods                : Modification list
        VarModThreshold             : int
        }

    module SearchDbParams =

        let toDomain (dtoSearchDbParams: SearchDbParams ) = 
            {
            Name                = dtoSearchDbParams.Name
            DbFolder            = dtoSearchDbParams.DbFolder
            FastaPath           = dtoSearchDbParams.FastaPath
            FastaHeaderToName   = parseProteinIdUsing dtoSearchDbParams.ParseProteinIDRegexPattern
            Protease            = Protease.toDomain dtoSearchDbParams.Protease
            MinMissedCleavages  = dtoSearchDbParams.MinMissedCleavages
            MaxMissedCleavages  = dtoSearchDbParams.MaxMissedCleavages
            MaxMass             = dtoSearchDbParams.MaxMass
            MinPepLength        = dtoSearchDbParams.MinPepLength
            MaxPepLength        = dtoSearchDbParams.MaxPepLength
            // valid symbol name=of isotopic label in label table i.e. #N15
            IsotopicMod         =  List.map IsotopicMod.toDomain dtoSearchDbParams.IsotopicMod
            MassMode            = dtoSearchDbParams.MassMode
            MassFunction        = MassMode.toDomain dtoSearchDbParams.MassMode
            FixedMods           = List.map Modification.toDomain dtoSearchDbParams.FixedMods            
            VariableMods        = List.map Modification.toDomain dtoSearchDbParams.VariableMods
            VarModThreshold     = dtoSearchDbParams.VarModThreshold
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
            Protease                : Protease
            MinMissedCleavages      : int
            MaxMissedCleavages      : int
            MaxMass                 : float
            MinPepLength            : int
            MaxPepLength            : int
            // valid symbol name of isotopic label in label table i.e. #N15
            IsotopicMod             : IsotopicMod list 
            MassMode                : SearchDB.MassMode
            FixedMods               : Modification list            
            VariableMods            : Modification list
            VarModThreshold         : int  
            // +/- ppm of ion m/z to obtain target peptides from SearchDB. 
            LookUpPPM               : float
            // lowest m/z, highest m/z
            MS2ScanRange            : float*float
            nTerminalSeries         : NTerminalSeries
            cTerminalSeries         : CTerminalSeries
            Andromeda               : AndromedaParams
        }

    module PeptideSpectrumMatchingParams = 
        
        let toDomain (dtoPeptideSpectrumMatchingParams: PeptideSpectrumMatchingParams ) :Domain.PeptideSpectrumMatchingParams = 
            {
                ExpectedMinimalCharge  = dtoPeptideSpectrumMatchingParams.ExpectedMinimalCharge   
                ExpectedMaximumCharge  = dtoPeptideSpectrumMatchingParams.ExpectedMaximumCharge   
                Width                  = dtoPeptideSpectrumMatchingParams.Width                   
                MinIntensity           = dtoPeptideSpectrumMatchingParams.MinIntensity            
                DeltaMinIntensity      = dtoPeptideSpectrumMatchingParams.DeltaMinIntensity       
                NrOfRndSpectra         = dtoPeptideSpectrumMatchingParams.NrOfRndSpectra          
                Protease               = Protease.toDomain dtoPeptideSpectrumMatchingParams.Protease                
                MinMissedCleavages     = dtoPeptideSpectrumMatchingParams.MinMissedCleavages      
                MaxMissedCleavages     = dtoPeptideSpectrumMatchingParams.MaxMissedCleavages      
                MaxMass                = dtoPeptideSpectrumMatchingParams.MaxMass                 
                MinPepLength           = dtoPeptideSpectrumMatchingParams.MinPepLength            
                MaxPepLength           = dtoPeptideSpectrumMatchingParams.MaxPepLength            
                IsotopicMod            = List.map IsotopicMod.toDomain dtoPeptideSpectrumMatchingParams.IsotopicMod             
                MassMode               = dtoPeptideSpectrumMatchingParams.MassMode                
                MassFunction           = MassMode.toDomain dtoPeptideSpectrumMatchingParams.MassMode            
                FixedMods              = List.map Modification.toDomain dtoPeptideSpectrumMatchingParams.FixedMods                   
                VariableMods           = List.map Modification.toDomain dtoPeptideSpectrumMatchingParams.VariableMods            
                VarModThreshold        = dtoPeptideSpectrumMatchingParams.VarModThreshold         
                LookUpPPM              = dtoPeptideSpectrumMatchingParams.LookUpPPM               
                MS2ScanRange           = dtoPeptideSpectrumMatchingParams.MS2ScanRange            
                nTerminalSeries        = NTerminalSeries.toDomain dtoPeptideSpectrumMatchingParams.nTerminalSeries         
                cTerminalSeries        = CTerminalSeries.toDomain dtoPeptideSpectrumMatchingParams.cTerminalSeries         
                Andromeda              = dtoPeptideSpectrumMatchingParams.Andromeda               
            }                                   


    type PEPEParams = 
        {
            QValueThreshold             : float
            PepValueThreshold           : float
            ParseProteinIDRegexPattern  : string
        }

    module PEPEParams = 

        let toDomain (dtoPepeParams: PEPEParams ) = 
            {
                QValueThreshold                 = dtoPepeParams.QValueThreshold    
                PepValueThreshold               = dtoPepeParams.PepValueThreshold  
                ParseProteinID                  = parseProteinIdUsing dtoPepeParams.ParseProteinIDRegexPattern    
            }

    type QuantificationParams = 
        {
            PerformLabeledQuantification : bool
            XicExtraction                : XicExtraction
            BaseLineCorrection           : BaseLineCorrection option
        }
    
    module QuantificationParams = 

        let toDomain (dtoQuantificationParams: QuantificationParams ) = 
            {
                PerformLabeledQuantification = dtoQuantificationParams.PerformLabeledQuantification
                XicExtraction                = dtoQuantificationParams.XicExtraction
                BaseLineCorrection           = dtoQuantificationParams.BaseLineCorrection
            }
        


    // Add by HLWeil
    //type ProteinInferenceParams = 
    //      ...
    
    //module ProteinInferenceParams = 
        
    //    let toDomain (dtoProteinInferenceParams : ProteinInferenceParams) = 
    //      ...