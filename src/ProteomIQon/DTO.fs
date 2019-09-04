namespace ProteomIQon

open BioFSharp
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open Domain
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Csv
open FSharpAux.IO.SchemaReader.Attribute

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
        | Pyro_Glu'GluNterm'
        | Pyro_Glu'GlnNterm'

    module Modification  =

        open BioFSharp.AminoAcids
        open BioFSharp.ModificationInfo

        let pyro_Glu'GluNterm' =
            createSearchModification "Pyro_Glu'Glu'" "27" "Pyro-glu from E" true "H2O"
                [Specific(Glu,ModLocation.Nterm);] SearchModType.Minus "pe"

        let pyro_Glu'GlnNterm' =
            createSearchModification "Pyro_Glu'Gln'" "28" "	Pyro-glu from Q" true "H3N"
                [Specific(Gln,ModLocation.Nterm);] SearchModType.Minus "pq"
             
        let toDomain modification = 
            match modification with
            | Acetylation'ProtNTerm'        -> SearchDB.Table.acetylation'ProtNTerm'
            | Carbamidomethyl'Cys'          -> SearchDB.Table.carbamidomethyl'Cys'
            | Oxidation'Met'                -> SearchDB.Table.oxidation'Met'
            | Phosphorylation'Ser'Thr'Tyr'  -> SearchDB.Table.phosphorylation'Ser'Thr'Tyr'
            | Pyro_Glu'GluNterm'            -> pyro_Glu'GluNterm'
            | Pyro_Glu'GlnNterm'            -> pyro_Glu'GlnNterm'

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

    type PreprocessingParams =
        {
            Compress                    : bool
            StartRetentionTime          : float option
            EndRetentionTime            : float option
            MS1PeakPicking              : PeakPicking
            MS2PeakPicking              : PeakPicking
        }

    module PreprocessingParams =

        let toDomain (dtoCentroidizationParams: PreprocessingParams ) : Domain.PreprocessingParams =
                {
                    Compress                    = dtoCentroidizationParams.Compress
                    StartRetentionTime          = dtoCentroidizationParams.StartRetentionTime
                    EndRetentionTime            = dtoCentroidizationParams.EndRetentionTime
                    MS1PeakPicking              = dtoCentroidizationParams.MS1PeakPicking
                    MS2PeakPicking              = dtoCentroidizationParams.MS2PeakPicking
                }

    type PeptideDBParams =
        {
        // name of database i.e. Creinhardtii_236_protein_full_labeled
        Name                        : string
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

    module PeptideDBParams =

        let toDomain (dtoSearchDbParams: PeptideDBParams ) =
            {
            Name                = dtoSearchDbParams.Name
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
            ChargeStateDeterminationParams: ChargeState.ChargeDetermParams
            LookUpPPM                     : float
            MS2ScanRange                  : float*float
            nTerminalSeries               : NTerminalSeries
            cTerminalSeries               : CTerminalSeries
            Andromeda                     : AndromedaParams
        }

    module PeptideSpectrumMatchingParams =

        let toDomain (dtoPeptideSpectrumMatchingParams: PeptideSpectrumMatchingParams ) :Domain.PeptideSpectrumMatchingParams =
            {
                ChargeStateDeterminationParams  = dtoPeptideSpectrumMatchingParams.ChargeStateDeterminationParams
                LookUpPPM                       = dtoPeptideSpectrumMatchingParams.LookUpPPM
                MS2ScanRange                    = dtoPeptideSpectrumMatchingParams.MS2ScanRange
                nTerminalSeries                 = NTerminalSeries.toDomain dtoPeptideSpectrumMatchingParams.nTerminalSeries
                cTerminalSeries                 = CTerminalSeries.toDomain dtoPeptideSpectrumMatchingParams.cTerminalSeries
                AndromedaParams                 = dtoPeptideSpectrumMatchingParams.Andromeda
            }

    type PeptideSpectrumMatchingResult =
        //{
        //    // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
        //    [<FieldAttribute(0)>]
        //    PSMId                        : string
        //    [<FieldAttribute(1)>]
        //    GlobalMod                    : int
        //    [<FieldAttribute(2)>]
        //    PepSequenceID                : int
        //    [<FieldAttribute(3)>]
        //    ModSequenceID                : int
        //    [<FieldAttribute(4)>]
        //    Label                        : int
        //    // ascending ms2 id (file specific)
        //    [<FieldAttribute(5)>]
        //    ScanNr                       : int
        //    [<FieldAttribute(6)>]
        //    ScanTime                     : float
        //    [<FieldAttribute(7)>]
        //    Charge                       : int
        //    [<FieldAttribute(8)>]
        //    PrecursorMZ                  : float
        //    [<FieldAttribute(9)>]
        //    TheoMass                     : float
        //    [<FieldAttribute(10)>]
        //    AbsDeltaMass                 : float
        //    [<FieldAttribute(11)>]
        //    PeptideLength                : int
        //    [<FieldAttribute(12)>]
        //    MissCleavages                : int
        //    [<FieldAttribute(13)>]
        //    SequestScore                 : float
        //    [<FieldAttribute(14)>]
        //    SequestNormDeltaBestToRest   : float
        //    [<FieldAttribute(15)>]
        //    SequestNormDeltaNext         : float
        //    [<FieldAttribute(16)>]
        //    AndroScore                   : float
        //    [<FieldAttribute(17)>]
        //    AndroNormDeltaBestToRest     : float
        //    [<FieldAttribute(18)>]
        //    AndroNormDeltaNext           : float
        //    [<FieldAttribute(19)>]
        //    StringSequence               : string
        //}
        {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
        [<FieldAttribute(0)>]
        PSMId                        : string;
        [<FieldAttribute(1)>]
        GlobalMod                    : int
        [<FieldAttribute(2)>]
        PepSequenceID                : int
        [<FieldAttribute(3)>]
        ModSequenceID                : int
        [<FieldAttribute(4)>]
        Label                        : int
        // ascending ms2 id (file specific)
        [<FieldAttribute(5)>]
        ScanNr                       : int;
        [<FieldAttribute(6)>]
        Charge                       : int;
        [<FieldAttribute(7)>]
        PrecursorMZ                  : float
        [<FieldAttribute(8)>]
        TheoMass                     : float;
        [<FieldAttribute(9)>]
        AbsDeltaMass                 : float;
        [<FieldAttribute(10)>]
        PeptideLength                : int;
        [<FieldAttribute(11)>]
        MissCleavages                : int;
        [<FieldAttribute(12)>]
        SequestScore                 : float;
        [<FieldAttribute(13)>]
        SequestNormDeltaBestToRest   : float;
        [<FieldAttribute(14)>]
        SequestNormDeltaNext         : float;
        [<FieldAttribute(15)>]
        AndroScore                   : float;
        [<FieldAttribute(16)>]
        AndroNormDeltaBestToRest     : float;
        [<FieldAttribute(17)>]
        AndroNormDeltaNext           : float;
        [<FieldAttribute(18)>]
        StringSequence               : string;
        [<FieldAttribute(19)>]
        ProteinNames                 : string
        }

    type PSMStatisticsParams =
        {
            QValueThreshold             : float
            PepValueThreshold           : float
            ParseProteinIDRegexPattern  : string
            KeepTemporaryFiles          : bool
        }

    module PSMStatisticsParams =

        let toDomain (dtoPSMStatisticsParams: PSMStatisticsParams ): Domain.PSMStatisticsParams =
            {
                QValueThreshold                 = dtoPSMStatisticsParams.QValueThreshold
                PepValueThreshold               = dtoPSMStatisticsParams.PepValueThreshold
                FastaHeaderToName               = parseProteinIdUsing dtoPSMStatisticsParams.ParseProteinIDRegexPattern
                KeepTemporaryFiles              = dtoPSMStatisticsParams.KeepTemporaryFiles
            }

    type PSMStatisticsResult = {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
        [<FieldAttribute(0)>]
        PSMId                        : string;
        [<FieldAttribute(1)>]
        GlobalMod                    : int
        [<FieldAttribute(2)>]
        PepSequenceID                : int
        [<FieldAttribute(3)>]
        ModSequenceID                : int
        [<FieldAttribute(4)>]
        Label                        : int
        // ascending ms2 id (file specific)
        [<FieldAttribute(5)>]
        ScanNr                       : int;
        [<FieldAttribute(6)>]
        Charge                       : int;
        [<FieldAttribute(7)>]
        PrecursorMZ                  : float
        [<FieldAttribute(8)>]
        TheoMass                     : float;
        [<FieldAttribute(9)>]
        AbsDeltaMass                 : float;
        [<FieldAttribute(10)>]
        PeptideLength                : int;
        [<FieldAttribute(11)>]
        MissCleavages                : int;
        [<FieldAttribute(12)>]
        SequestScore                 : float;
        [<FieldAttribute(13)>]
        SequestNormDeltaBestToRest   : float;
        [<FieldAttribute(14)>]
        SequestNormDeltaNext         : float;
        [<FieldAttribute(15)>]
        AndroScore                   : float;
        [<FieldAttribute(16)>]
        AndroNormDeltaBestToRest     : float;
        [<FieldAttribute(17)>]
        AndroNormDeltaNext           : float;
        [<FieldAttribute(18)>]
        PercolatorScore              : float;
        [<FieldAttribute(19)>]
        QValue                       : float;
        [<FieldAttribute(20)>]
        PEPValue                     : float;
        [<FieldAttribute(21)>]
        StringSequence               : string;
        [<FieldAttribute(22)>]
        ProteinNames                 : string;
        }

    type QuantificationParams =
        {
            PerformLabeledQuantification : bool
            XicExtraction                : XicExtraction
            //10 6 0.05
            BaseLineCorrection           : BaseLineCorrection option
        }

    module QuantificationParams =

        let toDomain (dtoQuantificationParams: QuantificationParams ): Domain.QuantificationParams =
            {
                PerformLabeledQuantification = dtoQuantificationParams.PerformLabeledQuantification
                XicExtraction                = dtoQuantificationParams.XicExtraction
                BaseLineCorrection           = dtoQuantificationParams.BaseLineCorrection
            }

    type ProteinInferenceParams =
          {
              ProteinIdentifierRegex : string
              Protein                : ProteinInference.IntegrationStrictness
              Peptide                : ProteinInference.PeptideUsageForQuantification
              GroupFiles             : bool
          }

    module ProteinInferenceParams =

        let toDomain (dtoProteinInferenceParams : ProteinInferenceParams): Domain.ProteinInferenceParams =
            {
                ProteinIdentifierRegex = dtoProteinInferenceParams.ProteinIdentifierRegex
                Protein                = dtoProteinInferenceParams.Protein
                Peptide                = dtoProteinInferenceParams.Peptide
                GroupFiles             = dtoProteinInferenceParams.GroupFiles
            }