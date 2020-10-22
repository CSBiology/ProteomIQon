namespace ProteomIQon

open BioFSharp
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open Domain
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute
open MzIO.Binary
open MzIO.Processing

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
            Compress                    : BinaryDataCompressionType
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
        {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
        [<FieldAttribute(0)>]
        PSMId                        : string
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
        ScanNr                       : int
        [<FieldAttribute(6)>]
        ScanTime                     : float
        [<FieldAttribute(7)>]
        Charge                       : int
        [<FieldAttribute(8)>]
        PrecursorMZ                  : float
        [<FieldAttribute(9)>]
        TheoMass                     : float
        [<FieldAttribute(10)>]
        AbsDeltaMass                 : float
        [<FieldAttribute(11)>]
        PeptideLength                : int
        [<FieldAttribute(12)>]
        MissCleavages                : int
        [<FieldAttribute(13)>]
        SequestScore                 : float
        [<FieldAttribute(14)>]
        SequestNormDeltaBestToRest   : float
        [<FieldAttribute(15)>]
        SequestNormDeltaNext         : float
        [<FieldAttribute(16)>]
        AndroScore                   : float
        [<FieldAttribute(17)>]
        AndroNormDeltaBestToRest     : float
        [<FieldAttribute(18)>]
        AndroNormDeltaNext           : float
        [<FieldAttribute(19)>]
        XtandemScore                 : float
        [<FieldAttribute(20)>]
        XtandemNormDeltaBestToRest   : float
        [<FieldAttribute(21)>]
        XtandemNormDeltaNext         : float
        [<FieldAttribute(22)>]
        StringSequence               : string
        [<FieldAttribute(23)>]
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
        PSMId                        : string
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
        ScanNr                       : int
        [<FieldAttribute(6)>]
        ScanTime                     : float
        [<FieldAttribute(7)>]
        Charge                       : int
        [<FieldAttribute(8)>]
        PrecursorMZ                  : float
        [<FieldAttribute(9)>]
        TheoMass                     : float
        [<FieldAttribute(10)>]
        AbsDeltaMass                 : float
        [<FieldAttribute(11)>]
        PeptideLength                : int
        [<FieldAttribute(12)>]
        MissCleavages                : int
        [<FieldAttribute(13)>]
        SequestScore                 : float
        [<FieldAttribute(14)>]
        SequestNormDeltaBestToRest   : float
        [<FieldAttribute(15)>]
        SequestNormDeltaNext         : float
        [<FieldAttribute(16)>]
        AndroScore                   : float
        [<FieldAttribute(17)>]
        AndroNormDeltaBestToRest     : float
        [<FieldAttribute(18)>]
        AndroNormDeltaNext           : float
        [<FieldAttribute(19)>]
        XtandemScore                 : float
        [<FieldAttribute(20)>]
        XtandemNormDeltaBestToRest   : float
        [<FieldAttribute(21)>]
        XtandemNormDeltaNext         : float
        [<FieldAttribute(22)>]
        PercolatorScore              : float
        [<FieldAttribute(23)>]
        QValue                       : float
        [<FieldAttribute(24)>]
        PEPValue                     : float
        [<FieldAttribute(25)>]
        StringSequence               : string
        [<FieldAttribute(26)>]
        ProteinNames                 : string
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

    ///
    type QuantificationSource = 
        | PSM
        | Alignment

    type QuantSourceConverter() = 
        inherit ConverterAttribute()
        override this.convertToObj = 
            Converter.Single(fun (str : string) -> 
                if (str) = "PSM" then QuantificationSource.PSM else QuantificationSource.Alignment 
                |> box
                )

    type TraceConverter() = 
        inherit ConverterAttribute()
        override this.convertToObj = 
            Converter.Single(fun (str : string) -> 
                let tmp = (str |> String.filter (fun x -> x <> '|' && x <> '[' && x <> ']' )).Trim() 
                if tmp = "" then 
                    [||]
                else
                    tmp.Split(';')
                    |> Array.map float
                |> box 
                )
                
    ///
    type QuantificationResult = {
        [<FieldAttribute(0)>]
        StringSequence                              : string
        [<FieldAttribute(1)>]
        GlobalMod                                   : int
        [<FieldAttribute(2)>]
        Charge                                      : int
        [<FieldAttribute(3)>]
        PepSequenceID                               : int
        [<FieldAttribute(4)>]
        ModSequenceID                               : int
        [<FieldAttribute(5)>]
        PrecursorMZ                                 : float
        [<FieldAttribute(6)>]
        MeasuredMass                                : float 
        [<FieldAttribute(7)>]
        TheoMass                                    : float
        [<FieldAttribute(8)>]
        AbsDeltaMass                                : float
        [<FieldAttribute(9)>]
        MeanPercolatorScore                         : float
        [<FieldAttribute(10)>]
        QValue                                      : float
        [<FieldAttribute(11)>]
        PEPValue                                    : float
        [<FieldAttribute(12)>]
        ProteinNames                                : string
        [<FieldAttribute(13)>]
        QuantMz_Light                               : float
        [<FieldAttribute(14)>]
        Quant_Light                                 : float
        [<FieldAttribute(15)>]
        MeasuredApex_Light                          : float 
        [<FieldAttribute(16)>]
        Seo_Light                                   : float
        [<FieldAttribute(17)>][<TraceConverter>]
        Params_Light                                : float []
        [<FieldAttribute(18)>]
        Difference_SearchRT_FittedRT_Light          : float
        [<FieldAttribute(19)>]
        KLDiv_Observed_Theoretical_Light            : float
        [<FieldAttribute(20)>]
        KLDiv_CorrectedObserved_Theoretical_Light   : float

        [<FieldAttribute(21)>]
        QuantMz_Heavy                               : float
        [<FieldAttribute(22)>]
        Quant_Heavy                                 : float
        [<FieldAttribute(23)>]
        MeasuredApex_Heavy                          : float
        [<FieldAttribute(24)>]
        Seo_Heavy                                   : float
        [<FieldAttribute(25)>][<TraceConverter>]
        Params_Heavy                                : float []        
        [<FieldAttribute(26)>]
        Difference_SearchRT_FittedRT_Heavy          : float
        [<FieldAttribute(27)>]
        KLDiv_Observed_Theoretical_Heavy            : float
        [<FieldAttribute(28)>]
        KLDiv_CorrectedObserved_Theoretical_Heavy   : float


        [<FieldAttribute(29)>]
        Correlation_Light_Heavy                     : float
        [<FieldAttribute(30)>][<QuantSourceConverter>]
        QuantificationSource                        : QuantificationSource

        [<FieldAttribute(31)>][<TraceConverter>]
        IsotopicPatternMz_Light                     : float []
        [<FieldAttribute(32)>][<TraceConverter>]
        IsotopicPatternIntensity_Observed_Light     : float []
        [<FieldAttribute(33)>][<TraceConverter>]
        IsotopicPatternIntensity_Corrected_Light    : float []
        [<FieldAttribute(34)>][<TraceConverter>]
        RtTrace_Light                               : float []
        [<FieldAttribute(35)>][<TraceConverter>]
        IntensityTrace_Observed_Light               : float []
        [<FieldAttribute(36)>][<TraceConverter>]
        IntensityTrace_Corrected_Light              : float []
        [<FieldAttribute(37)>][<TraceConverter>]
        IsotopicPatternMz_Heavy                     : float []
        [<FieldAttribute(38)>][<TraceConverter>]
        IsotopicPatternIntensity_Observed_Heavy     : float []
        [<FieldAttribute(39)>][<TraceConverter>]
        IsotopicPatternIntensity_Corrected_Heavy    : float []
        [<FieldAttribute(40)>][<TraceConverter>]
        RtTrace_Heavy                               : float []
        [<FieldAttribute(41)>][<TraceConverter>]
        IntensityTrace_Observed_Heavy               : float []
        [<FieldAttribute(42)>][<TraceConverter>]
        IntensityTrace_Corrected_Heavy              : float []
        }

    module QuantificationResult = 
        /// Retrieves the scan time based on the fitted parameter values (HULQ output).
        let getTargetIntensity (qp:QuantificationResult) = 
            try
            if qp.GlobalMod = 0 then
                qp.Quant_Light
            else
                qp.Quant_Heavy
            with
            | _ -> nan

        /// Retrieves the scan time based on the fitted parameter values (HULQ output).
        let getTargetScanTime (qp:QuantificationResult) = 
            try
            if qp.GlobalMod = 0 then
                qp.Params_Light.[1] 
            else
                qp.Params_Heavy.[1] 
            with
            | _ -> nan

        /// Retrieves the scan time based on the fitted parameter values (HULQ output).
        let getTargetStabw (qp:QuantificationResult) = 
            try
            if qp.GlobalMod = 0 then
                qp.Params_Light.[2] 
            else
                qp.Params_Heavy.[2] 
            with
            | _ -> nan

        /// Retrieves the scan time based on the fitted parameter values (HULQ output).
        let getTargetScanTimeDifference (qp:QuantificationResult) = 
            try
            if qp.GlobalMod = 0 then
                qp.Difference_SearchRT_FittedRT_Light
            else
                qp.Difference_SearchRT_FittedRT_Heavy
            with
            | _ -> nan

        /// Retrieves the scan time based on the fitted parameter values (HULQ output).
        let tryTargetGetScanTime (qp:QuantificationResult) = 
            try
            if qp.GlobalMod = 0 then
                qp.Params_Light.[1] 
                |> float
                |> Some
            else
                qp.Params_Heavy.[1] 
                |> float
                |> Some
            with
            | _ -> None

        /// 
        let getTargetRtTrace (qp:QuantificationResult) = 
            try
            if qp.GlobalMod = 0 then
                qp.RtTrace_Light
            else
                qp.RtTrace_Heavy 
            with
            | _ -> [||]

        /// 
        let getTargetIntensityTrace (qp:QuantificationResult) = 
            try
            if qp.GlobalMod = 0 then
                qp.IntensityTrace_Corrected_Light
            else
                qp.IntensityTrace_Corrected_Heavy
            with
            | _ -> [||]

        /// 
        let getIsotopicPatternMz (qp:QuantificationResult) = 
            try
            if qp.GlobalMod = 0 then
                qp.IsotopicPatternMz_Light
            else
                qp.IsotopicPatternMz_Heavy
            with
            | _ -> [||]

        /// 
        let getIsotopicPatternIntensity_Observed (qp:QuantificationResult) = 
            try
            if qp.GlobalMod = 0 then
                qp.IsotopicPatternIntensity_Corrected_Light
            else
                qp.IsotopicPatternIntensity_Corrected_Heavy
            with
            | _ -> [||]

    ///
    type AlignmentParams = {
        Placeholder : bool 
        }

    ///
    type AlignmentResult = 
        {
            [<FieldAttribute(0)>]
            StringSequence               : string
            [<FieldAttribute(1)>]
            GlobalMod                    : int
            [<FieldAttribute(2)>]
            Charge                       : int
            [<FieldAttribute(3)>]
            PepSequenceID                : int
            [<FieldAttribute(4)>]
            ModSequenceID                : int
            [<FieldAttribute(5)>]
            Mz                           : float
            [<FieldAttribute(6)>]
            ProteinNames                 : string
            [<FieldAttribute(7)>]
            PredictedScanTime            : float
            [<FieldAttribute(8)>]
            ScanTime_SourceFile          : float
            [<FieldAttribute(9)>][<TraceConverter>]
            RtTrace_SourceFile           : float []
            [<FieldAttribute(10)>][<TraceConverter>]
            IntensityTrace_SourceFile    : float []
            [<FieldAttribute(11)>][<TraceConverter>]
            IsotopicPatternMz_SourceFile                    : float []            
            [<FieldAttribute(12)>][<TraceConverter>]
            IsotopicPatternIntensity_Observed_SourceFile    : float []       
        } 

    /////
    //type AlignmentModelMetrics<'Metrics> = 
    //    {
    //    Metrics                             : 'Metrics
    //    Sequence                            : string []
    //    GlobalMod                           : int []
    //    Charge                              : int []
    //    PepSequenceID                       : int []
    //    ModSequenceID                       : int []
    //    X_Intensities                       : float []
    //    X_Stabw                             : float []        
    //    X_Test                              : float []
    //    X_IsotopicPatternMz                 : float [][]
    //    X_IsotopicPatternIntensity_Observed : float [][]    
    //    Y_Test                              : float []
    //    YHat_Test                           : float []
    //    YHat_Refined_Test                   : float []
    //    Y_IsotopicPatternMz                 : float [][]
    //    Y_IsotopicPatternIntensity_Observed : float [][]       
    //    DtwDistanceBefore                   : float []
    //    DtwDistanceAfter                    : float []
    //    }

    ///
    type AlignmentMetricsDTO = 
       {
           [<FieldAttribute(0)>]
           Sequence                             : string
           [<FieldAttribute(1)>]
           GlobalMod                            : int
           [<FieldAttribute(2)>]
           Charge                               : int
           [<FieldAttribute(3)>]
           PepSequenceID                        : int
           [<FieldAttribute(4)>]
           ModSequenceID                        : int
           [<FieldAttribute(5)>]
           X_FileName                           : string 
           [<FieldAttribute(6)>]
           X_Intensities                        : float 
           [<FieldAttribute(7)>]
           X_Stabw                              : float
           [<FieldAttribute(8)>]
           X_Test                               : float 
           [<FieldAttribute(9)>] [<TraceConverter>]
           X_IsotopicPatternMz                  : float []
           [<FieldAttribute(10)>] [<TraceConverter>]
           X_IsotopicPatternIntensity_Observed  : float []
           [<FieldAttribute(11)>] [<TraceConverter>]
           X_RtTrace                            : float []
           [<FieldAttribute(12)>] [<TraceConverter>]
           X_IntensityTrace                     : float []   
           [<FieldAttribute(13)>]
           Y_Test                               : float 
           [<FieldAttribute(14)>]
           YHat_Test                            : float 
           [<FieldAttribute(15)>]
           YHat_Refined_Test                    : float
           [<FieldAttribute(16)>]
           Y_Intensity                          : float
           [<FieldAttribute(17)>] [<TraceConverter>]
           Y_IsotopicPatternMz                  : float []
           [<FieldAttribute(18)>] [<TraceConverter>]
           Y_IsotopicPatternIntensity_Observed  : float [] 
           [<FieldAttribute(19)>] [<TraceConverter>]
           Y_RtTrace                            : float []
           [<FieldAttribute(20)>] [<TraceConverter>]
           Y_IntensityTrace                     : float []
           [<FieldAttribute(21)>]
           DtwDistanceBefore                    : float 
       }

    type AlignmentBasedQuantificationParams =
        {
            PerformLabeledQuantification : bool
            PerformLocalWarp             : bool
            XicExtraction                : XicExtraction
            //10 6 0.05
            BaseLineCorrection           : BaseLineCorrection option
        }

    module AlignmentBasedQuantificationParams =

        let toDomain (dtoQuantificationParams: AlignmentBasedQuantificationParams ): Domain.AlignmentBasedQuantificationParams =
            {
                PerformLabeledQuantification = dtoQuantificationParams.PerformLabeledQuantification
                PerformLocalWarp             = dtoQuantificationParams.PerformLocalWarp
                XicExtraction                = dtoQuantificationParams.XicExtraction
                BaseLineCorrection           = dtoQuantificationParams.BaseLineCorrection
            }

    type ProteinInferenceParams =
          {
              ProteinIdentifierRegex : string
              Protein                : ProteinInference.IntegrationStrictness
              Peptide                : ProteinInference.PeptideUsageForQuantification
              GroupFiles             : bool
              GetQValue              : QValueMethod
          }

    module ProteinInferenceParams =

        let toDomain (dtoProteinInferenceParams: ProteinInferenceParams): Domain.ProteinInferenceParams =
            {
                ProteinIdentifierRegex = dtoProteinInferenceParams.ProteinIdentifierRegex
                Protein                = dtoProteinInferenceParams.Protein
                Peptide                = dtoProteinInferenceParams.Peptide
                GroupFiles             = dtoProteinInferenceParams.GroupFiles
                GetQValue              = dtoProteinInferenceParams.GetQValue
            }

    type SpectralLibraryParams =
        {
            ChargeList          : float list
            MatchingTolerancePPM: float
        }

    module SpectralLibraryParams =

        let toDomain (dtoSpectralLibraryParams: SpectralLibraryParams): Domain.SpectralLibraryParams =
            {
                ChargeList           = dtoSpectralLibraryParams.ChargeList
                MatchingTolerancePPM = dtoSpectralLibraryParams.MatchingTolerancePPM
            }


    type TableSortParams =
        {
            SeparatorIn                 : string
            SeparatorOut                : char
            EssentialFields             : EssentialFields
            QuantFieldsToFilterOn       : FilterOnField[]
            ProtFieldsToFilterOn        : FilterOnField[]
            QuantColumnsOfInterest      : string[]
            ProtColumnsOfInterest       : string[]
            StatisticalMeasurements     : (string*StatisticalMeasurement)[]
            AggregatorFunction          : AggregationMethod
            AggregatorFunctionIntensity : AggregationMethod
            AggregatorPepToProt         : AggregationMethod
            Tukey                       : (string*float*Transform) []
        }

    module TableSortParams =
        
        let inline toDomain (dtoTableSortParams: TableSortParams): Domain.TableSortParams =
            {
                SeparatorIn                 = dtoTableSortParams.SeparatorIn
                SeparatorOut                = dtoTableSortParams.SeparatorOut
                EssentialFields             = dtoTableSortParams.EssentialFields
                QuantFieldsToFilterOn       = dtoTableSortParams.QuantFieldsToFilterOn
                ProtFieldsToFilterOn        = dtoTableSortParams.ProtFieldsToFilterOn
                QuantColumnsOfInterest      = dtoTableSortParams.QuantColumnsOfInterest
                ProtColumnsOfInterest       = dtoTableSortParams.ProtColumnsOfInterest
                StatisticalMeasurements     = dtoTableSortParams.StatisticalMeasurements
                AggregatorFunction          = dtoTableSortParams.AggregatorFunction
                AggregatorFunctionIntensity = dtoTableSortParams.AggregatorFunctionIntensity
                AggregatorPepToProt         = dtoTableSortParams.AggregatorPepToProt
                Tukey                       = dtoTableSortParams.Tukey
            }
    type ConsensusSpectralLibraryParams =
        {
            RTTolerance: float
            iRTPeptides: string list
        }

    module ConsensusSpectralLibraryParams =

        let toDomain (dtoConsensusSpectralLibraryParams: ConsensusSpectralLibraryParams): Domain.ConsensusSpectralLibraryParams =
            {
                RTTolerance = dtoConsensusSpectralLibraryParams.RTTolerance
                iRTPeptides = dtoConsensusSpectralLibraryParams.iRTPeptides
            }

    type SWATHAnalysisParams =
        {
            PeptideList         : string [] option
            MatchingTolerancePPM: float
            QueryOffsetRange    : float
            SpectrumSelectionF  : SpectrumSelection
            AggregationF        : AggregationMethod
            XicProcessing       : XicProcessing
        }

    module SWATHAnalysisParams =

        let toDomain (dtoSWATHAnalysisParams: SWATHAnalysisParams): Domain.SWATHAnalysisParams =
            {
                PeptideList          = dtoSWATHAnalysisParams.PeptideList
                MatchingTolerancePPM = dtoSWATHAnalysisParams.MatchingTolerancePPM
                QueryOffsetRange     = dtoSWATHAnalysisParams.QueryOffsetRange
                SpectrumSelectionF   = dtoSWATHAnalysisParams.SpectrumSelectionF
                AggregationF         = dtoSWATHAnalysisParams.AggregationF
                XicProcessing        = dtoSWATHAnalysisParams.XicProcessing
            }

    type MzTABParams =
        {
            ExperimentNames      : (string*int)[]
            StudyVariables       : (string*int[]*int)[]
            SearchEngineNamesProt: (Ontologies.SearchEngineScore*string*int)[]
            SearchEngineNamesPep : (Ontologies.SearchEngineScore*string*int)[]
            SearchEngineNamesPSM : (Ontologies.SearchEngineScore*string*int)[]
            MetaData             : MetaDataSection
        }

    module MzTABParams =

        let toDomain (dtoMzTABParams: MzTABParams): Domain.MzTABParams =
            {
                ExperimentNames       = dtoMzTABParams.ExperimentNames
                StudyVariables        = dtoMzTABParams.StudyVariables
                SearchEngineNamesProt = dtoMzTABParams.SearchEngineNamesProt
                SearchEngineNamesPep  = dtoMzTABParams.SearchEngineNamesPep
                SearchEngineNamesPSM  = dtoMzTABParams.SearchEngineNamesPSM
                MetaData              = dtoMzTABParams.MetaData
            }

    type FragmentIon = {
        Charge                            : float
        Iontype                           : string
        Number                            : int
        // Brauch man glaub ich nicht weil das in der Quantinfo steckt
        GlobalMod                         : int
        CountAbsolute                     : int
        CountFraction                     : float
        MeanFragMz                        : float 
        CvMeanFragMz                      : float
        MaxIntensity                      :float
        MinIntensity                      :float
        MeanRelativeIntensity_Total       : float
        CVRelativeIntensity_Total         : float
        MeanRelativeIntensity_Frags       : float
        CVRelativeIntensity_Frags         : float
        }
    
    ///
    type PeptideIon = {
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
        Quant                                       : float
        RelativeQuant                               : float
        MeasuredApex                                : float
        RelativeMeasuredApex                        : float
        Seo                                         : float
        ScanTime                                    : float
        ElutionWidth                                : float
        Fragments                                   : FragmentIon list
        }
 