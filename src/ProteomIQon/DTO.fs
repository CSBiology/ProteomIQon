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
        | Trypsin_P
        | LysC
        | LysC_P
        | Chymotrypsin
        | PepsinA


    module Protease =
        
        open BioFSharp.AminoAcids

        let toDomain protease =
            match protease with
            | Trypsin   -> Digestion.Table.getProteaseBy "Trypsin"
            | Trypsin_P -> 
                BioFSharp.Digestion.createProtease "Trypsin/P" 
                    (let _p1 = [AminoAcid.Lys;AminoAcid.Arg] |> Set.ofList 
                    fun p4 p3 p2 p1 p1' p2' -> 
                    match p1,p1' with
                    | Some a1,Some a1' -> _p1.Contains(a1)
                    | _   -> false                     
                    )    
            | LysC      -> 
                BioFSharp.Digestion.createProtease "LysC" 
                    (let _p1 = [AminoAcid.Lys] |> Set.ofList 
                    fun p4 p3 p2 p1 p1' p2' -> 
                    match p1,p1' with
                    | Some a1,Some a1' -> _p1.Contains(a1) && not (a1' = AminoAcid.Pro)
                    | _   -> false                     
                    )
            | LysC_P      -> 
                BioFSharp.Digestion.createProtease "LysC/P" 
                    (let _p1 = [AminoAcid.Lys] |> Set.ofList 
                    fun p4 p3 p2 p1 p1' p2' -> 
                    match p1,p1' with
                    | Some a1,Some a1' -> _p1.Contains(a1)
                    | _   -> false                     
                    )
            | Chymotrypsin      -> 
                BioFSharp.Digestion.createProtease "Chymotrypsin" 
                    (let _p1 = [AminoAcid.Phe;AminoAcid.Tyr;AminoAcid.Trp;AminoAcid.Leu] |> Set.ofList 
                    fun p4 p3 p2 p1 p1' p2' -> 
                    match p1,p1' with
                    | Some a1,Some a1' -> _p1.Contains(a1) && not (a1' = AminoAcid.Pro)
                    | _   -> false                     
                    )
            | PepsinA      -> 
                BioFSharp.Digestion.createProtease "PepsinA" 
                    (let _p1 = [AminoAcid.Phe;AminoAcid.Leu] |> Set.ofList 
                    fun p4 p3 p2 p1 p1' p2' -> 
                    match p1,p1' with
                    | Some a1,Some a1' -> _p1.Contains(a1)
                    | _   -> false                     
                    )
    
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
        | C13
        | O17
        | O18
        | D

    module IsotopicMod =
       
        open BioFSharp.Elements

        module Elements' = 
            let C13 = Di (createDi "C13" (Isotopes.Table.C13,Isotopes.Table.C12.NatAbundance) (Isotopes.Table.C12,Isotopes.Table.C13.NatAbundance) )
            let O17 = Tri (createTri "O17" (Isotopes.Table.O17,Isotopes.Table.O16.NatAbundance) (Isotopes.Table.O16,Isotopes.Table.O17.NatAbundance) (Isotopes.Table.O18,Isotopes.Table.O18.NatAbundance) )
            let O18 = Tri (createTri "O18" (Isotopes.Table.O18,Isotopes.Table.O16.NatAbundance) (Isotopes.Table.O16,Isotopes.Table.O18.NatAbundance) (Isotopes.Table.O17,Isotopes.Table.O17.NatAbundance) )
            let D = Di (createDi "D" (Isotopes.Table.H2,Isotopes.Table.H1.NatAbundance) (Isotopes.Table.H1,Isotopes.Table.H2.NatAbundance) )
        
        let toDomain isoMod =
            match isoMod with
            | N15 -> (SearchDB.createSearchInfoIsotopic "N15" Elements.Table.N Elements.Table.Heavy.N15)
            | C13 -> (SearchDB.createSearchInfoIsotopic "C13" Elements.Table.C Elements'.C13)
            | O17 -> (SearchDB.createSearchInfoIsotopic "O17" Elements.Table.O Elements'.O17)
            | O18 -> (SearchDB.createSearchInfoIsotopic "O18" Elements.Table.O Elements'.O18)
            | D   -> (SearchDB.createSearchInfoIsotopic "D" Elements.Table.H Elements'.D)

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

    let tryParseProteinIdUsing regex =
        match regex with
        | "ID" | "id" | "Id" | "" ->
            id >> Some
        | pattern ->
            fun (inp : string)  -> 
                let tmp = System.Text.RegularExpressions.Regex.Match(inp,pattern)
                if tmp.Success then 
                    Some tmp.Value
                else None

    type Compression =
        | NoCompression = 0
        | ZLib = 1
        | NumPress = 2
        | NumPressZLib = 3

    module Compression =
        let toDomain (compression: Compression) =
            match compression with
            | Compression.NoCompression -> BinaryDataCompressionType.NoCompression
            | Compression.ZLib          -> BinaryDataCompressionType.ZLib
            | Compression.NumPress      -> BinaryDataCompressionType.NumPress
            | Compression.NumPressZLib  -> BinaryDataCompressionType.NumPressZLib

    type UseModifiedPeptides =
        | All 
        | No 
        | UseOnly of seq<Modification>

    module UseModifiedPeptides =

        let toDomain (useModifiedPeptides:UseModifiedPeptides)  = 
            match useModifiedPeptides with 
            | All -> 
                fun (x:string) -> true
            | No -> 
                fun (x:string) -> x.Contains "[" |> not
            | UseOnly mods -> 
                let searchMods = 
                    mods 
                    |> Seq.map (Modification.toDomain >> (fun x -> x.XModCode))
                    |> Set.ofSeq
                fun (x:string) -> 
                    let pattern = @"(?<=\[)(\w*)(?=\])"
                    let matches = System.Text.RegularExpressions.Regex.Matches(x,pattern)
                    // using this in favor of the out commented code since there seems to be a netstandard related issue that 
                    // hinders resolving that matches implements IEnumerable.
                    matches.Count = 0 || [for m in matches do searchMods.Contains m.Value] |> Seq.forall id
                    //Seq.isEmpty matches || (matches |> Seq.forall (fun m -> searchMods.Contains m.Value)) 
    
    open FSharp.Stats

    type NumericTransform =
        | Log2
        | Substract of float 
        | Add of float
        | DivideBy of float
        | MultiplyBy of float 

    module NumericTransform =     
        
        let toDomain transform =
            match transform with 
            | Log2          -> log2 
            | Substract v   -> fun x -> x - v 
            | Add v         -> fun x -> v + x
            | DivideBy v    -> fun x -> x / v
            | MultiplyBy v  -> fun x -> x * v 
            
    type NumericAggregation = 
        | Mean 
        | Median 
        | Sum
    
    module NumericAggregation =     
        
        let toDomain aggregation  = 
            match aggregation with 
            | Mean      -> Seq.mean 
            | Median    -> Seq.median 
            | Sum       -> Seq.sum 

    type DispersionMeasure =
        | SEM
        | StDev
        | CV
    
    module DispersionMeasure =     
        let toDomain dipsMeasure  = 
            match dipsMeasure with 
            | SEM   -> 
                fun (x:seq<float>) ->
                    let stDev = Seq.stDev x
                    stDev / (sqrt (x |> Seq.length |> float))
            | StDev -> Seq.stDev 
            | CV    -> Seq.cv 

    type NumericFilter = 
        | IsBiggerThan of float 
        | IsSmallerThan of float

    module NumericFilter =     
        let toDomain filter = 
            match filter with 
            | IsBiggerThan v -> fun x -> x > v
            | IsSmallerThan v -> fun x -> x < v

    type GroupFilter = 
        | Tukey of float  
        | Stdev of float
        | TopX  of int


    module GroupFilter =
        let toDomain gf =
            match gf with 
            | Tukey threshold -> 
                fun (values:seq<float>) -> 
                    let v = Array.ofSeq values
                    let tukey = FSharp.Stats.Signal.Outliers.tukey threshold
                    match tukey v with 
                    | Intervals.Interval.ClosedInterval (lower, upper) -> 
                        (fun v -> v <= upper && v >= lower)
                    | _ -> fun v -> false
            | Stdev threshold ->  
                fun values -> 
                    let mean = Seq.mean values
                    let stdev = Seq.stDev values
                    (fun v -> v <= (mean+stdev*threshold) && v >= (mean-stdev*threshold)
                    )
            | TopX count ->
                fun (values:seq<float>) -> 
                    let topX =
                        values
                        |> Seq.sortDescending
                        |> Seq.take count
                    fun v -> topX |> Seq.contains v


    module LabeledProteinQuantification = 
        
        type LabeledTransforms = 
                {
                    Light: NumericTransform option 
                    Heavy: NumericTransform option
                    Ratio: NumericTransform option
                }   
            
        type LabeledAggregations = 
            {
                Light: NumericAggregation
                Heavy: NumericAggregation
                Ratio: NumericAggregation
            }   
        
        type LabeledSingleFilters = 
            {
                Light: seq<NumericFilter> option
                Heavy: seq<NumericFilter> option
                Ratio: seq<NumericFilter> option
            }   

        type LabeledGroupFilters = 
            {
                Light: seq<GroupFilter> option
                Heavy: seq<GroupFilter> option
                Ratio: seq<GroupFilter> option
            }   

        type AggregationParams = 
            {
                LabeledTransform        : LabeledTransforms option
                LabeledSingleFilters    : LabeledSingleFilters option
                LabeledGroupFilters     : LabeledGroupFilters option
                LabeledAggregation      : LabeledAggregations 
            }
        
        module AggregationParams = 

            let toDomain (aggP:AggregationParams) :Domain.LabeledProteinQuantification.AggregationParams = 
                let labeledTransforms :Domain.LabeledProteinQuantification.LabeledTransforms = 
                    match aggP.LabeledTransform with
                    | Some t -> 
                        let unpackTransform t = 
                            match t with 
                            | Some t -> NumericTransform.toDomain t
                            | None -> id
                        {Light = unpackTransform t.Light;Heavy = unpackTransform t.Heavy;Ratio = unpackTransform t.Ratio}
                    | None -> 
                        {Light = id;Heavy = id;Ratio = id}
                let labeledSingleFilters :Domain.LabeledProteinQuantification.LabeledSingleFilters = 
                    match aggP.LabeledSingleFilters with
                    | Some f -> 
                        let unpackTransform (fs:seq<NumericFilter> option) = 
                            match fs with 
                            | Some t -> t |> Seq.map NumericFilter.toDomain 
                            | None -> Seq.empty
                        {Light = unpackTransform f.Light;Heavy = unpackTransform f.Heavy;Ratio = unpackTransform f.Ratio}
                    | None -> 
                        {Light = Seq.empty;Heavy = Seq.empty;Ratio = Seq.empty}
                let labeledGroupFilters :Domain.LabeledProteinQuantification.LabeledGroupFilters = 
                    match aggP.LabeledGroupFilters with
                    | Some f -> 
                        let unpackTransform (fs:seq<GroupFilter> option) = 
                            match fs with 
                            | Some t -> t |> Seq.map GroupFilter.toDomain 
                            | None -> Seq.empty
                        {Light = unpackTransform f.Light;Heavy = unpackTransform f.Heavy;Ratio = unpackTransform f.Ratio}
                    | None -> 
                        {Light = Seq.empty;Heavy = Seq.empty;Ratio = Seq.empty}
                let labeledAggregations :Domain.LabeledProteinQuantification.LabeledAggregations =                     
                    {
                        Light = aggP.LabeledAggregation.Light   |> NumericAggregation.toDomain 
                        Heavy = aggP.LabeledAggregation.Heavy   |> NumericAggregation.toDomain
                        Ratio = aggP.LabeledAggregation.Ratio   |> NumericAggregation.toDomain
                     }     
                {
                    LabeledTransform        = labeledTransforms 
                    LabeledSingleFilters    = labeledSingleFilters 
                    LabeledGroupFilters     = labeledGroupFilters 
                    LabeledAggregation      = labeledAggregations
                }

    module LabelFreeQuantification = 
        
        type Transforms = 
                {
                    Light: NumericTransform option 
                }   
            
        type Aggregations = 
            {
                Light: NumericAggregation
            }   
        
        type SingleFilters = 
            {
                Light: seq<NumericFilter> option
            }   

        type GroupFilters = 
            {
                Light: seq<GroupFilter> option
            }   

        type AggregationParams = 
            {
                Transform        : Transforms option
                SingleFilters    : SingleFilters option
                GroupFilters     : GroupFilters option
                Aggregation      : Aggregations 
            }
        
        module AggregationParams = 

            let toDomain (aggP:AggregationParams) :Domain.LabelFreeQuantification.AggregationParams = 
                let transforms :Domain.LabelFreeQuantification.Transforms = 
                    match aggP.Transform with
                    | Some t -> 
                        let unpackTransform t = 
                            match t with 
                            | Some t -> NumericTransform.toDomain t
                            | None -> id
                        {Light = unpackTransform t.Light}
                    | None -> 
                        {Light = id}
                let singleFilters :Domain.LabelFreeQuantification.SingleFilters = 
                    match aggP.SingleFilters with
                    | Some f -> 
                        let unpackTransform (fs:seq<NumericFilter> option) = 
                            match fs with 
                            | Some t -> t |> Seq.map NumericFilter.toDomain 
                            | None -> Seq.empty
                        {Light = unpackTransform f.Light}
                    | None -> 
                        {Light = Seq.empty}
                let groupFilters :Domain.LabelFreeQuantification.GroupFilters = 
                    match aggP.GroupFilters with
                    | Some f -> 
                        let unpackTransform (fs:seq<GroupFilter> option) = 
                            match fs with 
                            | Some t -> t |> Seq.map GroupFilter.toDomain 
                            | None -> Seq.empty
                        {Light = unpackTransform f.Light}
                    | None -> 
                        {Light = Seq.empty}
                let aggregations :Domain.LabelFreeQuantification.Aggregations =                     
                    {
                        Light = aggP.Aggregation.Light |> NumericAggregation.toDomain 
                    }     
                {
                    Transform        = transforms 
                    SingleFilters    = singleFilters 
                    GroupFilters     = groupFilters 
                    Aggregation      = aggregations
                }
/// 
module Dto =

    type PreprocessingParams =
        {
            Compress                    : Compression
            StartRetentionTime          : float option
            EndRetentionTime            : float option
            MS1PeakPicking              : PeakPicking
            MS2PeakPicking              : PeakPicking
        }

    module PreprocessingParams =

        let toDomain (dtoCentroidizationParams: PreprocessingParams ) : Domain.PreprocessingParams =
                {
                    Compress                    = Compression.toDomain dtoCentroidizationParams.Compress
                    StartRetentionTime          = dtoCentroidizationParams.StartRetentionTime
                    EndRetentionTime            = dtoCentroidizationParams.EndRetentionTime
                    MS1PeakPicking              = dtoCentroidizationParams.MS1PeakPicking
                    MS2PeakPicking              = dtoCentroidizationParams.MS2PeakPicking
                }

    
    type PeptideDBParams =
        {
        // name of database i.e. Creinhardtii_236_protein_full_labeled
        Name                        : string
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
            nTerminalSeries               : NTerminalSeries
            cTerminalSeries               : CTerminalSeries
            Andromeda                     : AndromedaParams
        }

    module PeptideSpectrumMatchingParams =

        let toDomain (dtoPeptideSpectrumMatchingParams: PeptideSpectrumMatchingParams ) :Domain.PeptideSpectrumMatchingParams =
            {
                ChargeStateDeterminationParams  = dtoPeptideSpectrumMatchingParams.ChargeStateDeterminationParams
                LookUpPPM                       = dtoPeptideSpectrumMatchingParams.LookUpPPM
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
        }


    type PSMStatisticsParams = 
        {
            Threshold                   : Domain.Threshold 
            ParseProteinIDRegexPattern  : string 
            KeepTemporaryFiles          : bool
        }


    module PSMStatisticsParams =

        let toDomain (dtoPSMStatisticsParams: PSMStatisticsParams ): Domain.PSMStatisticsParams =
            {
                Threshold           = dtoPSMStatisticsParams.Threshold
                FastaHeaderToName   = parseProteinIdUsing dtoPSMStatisticsParams.ParseProteinIDRegexPattern
                KeepTemporaryFiles  = dtoPSMStatisticsParams.KeepTemporaryFiles
            }

    type PSMStatisticsResultFragpipe = {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space separated by '_'
        [<FieldAttribute(0)>]
        PSMId                        : string
        [<FieldAttribute(1)>]
        PepSequenceID                : int
        [<FieldAttribute(2)>]
        ModSequenceID                : int
        [<FieldAttribute(3)>]
        ScanTime                     : float
        [<FieldAttribute(4)>]
        Charge                       : int
        [<FieldAttribute(5)>]
        PrecursorMZ                  : float
        [<FieldAttribute(6)>]
        TheoMass                     : float
        [<FieldAttribute(7)>]
        AbsDeltaMass                 : float
        [<FieldAttribute(8)>]
        IonMobility                  : float
        [<FieldAttribute(9)>]
        PeptideLength                : int
        [<FieldAttribute(10)>]
        MissCleavages                : int
        [<FieldAttribute(11)>]
        Expectscore                  : float
        [<FieldAttribute(12)>]
        Hyperscore                   : float
        [<FieldAttribute(13)>]
        Probability                  : float
        [<FieldAttribute(14)>]
        StringSequence               : string
        [<FieldAttribute(15)>]
        ProteinNames                 : string
        [<FieldAttribute(16)>]
        GlobalMod                    : int
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
        ModelScore                   : float
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
            PerformLabeledQuantification : Labeling
            FragPipe                     : bool
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
                let tmp = (str |> String.filter (fun x -> x <> '|' && x <> '[' && x <> ']' && x <> '"' )).Trim() 
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
        [<FieldAttribute(43)>]
        AlignmentScore                              : float
        [<FieldAttribute(44)>]
        AlignmentQValue                             : float
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

        /// Retrieves the scan time based on the fitted parameter values (HULQ output).
        let tryLightGetScanTime (qp:QuantificationResult) = 
            try
                qp.Params_Light.[1] 
                |> float
                |> Some
            with
            | _ -> None

        /// Retrieves the scan time based on the fitted parameter values (HULQ output).
        let tryHeavyGetScanTime (qp:QuantificationResult) = 
            try
                qp.Params_Heavy.[1] 
                |> float
                |> Some
            with
            | _ -> None

        /// Retrieves the scan time based on the fitted parameter values (HULQ output).
        let tryLightGetStabw (qp:QuantificationResult) = 
            try
                qp.Params_Light.[2] 
                |> float
                |> Some
            with
            | _ -> None

        /// Retrieves the scan time based on the fitted parameter values (HULQ output).
        let tryHeavyGetStabw (qp:QuantificationResult) = 
            try
                qp.Params_Heavy.[2] 
                |> float
                |> Some
            with
            | _ -> None
    
    ///
    type LabelEfficiencyCalculationParams =
        {
            LowerBound: float
            UpperBound: float
        }

    module LabelEfficiencyCalculationParams =
    
        let toDomain (dtoLabelEfficiencyParams: LabelEfficiencyCalculationParams ): Domain.LabelEfficiencyCalculationParams =
            {
                LowerBound = dtoLabelEfficiencyParams.LowerBound
                UpperBound = dtoLabelEfficiencyParams.UpperBound
            }
        
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
            [<FieldAttribute(9)>]
            ApexIntensity_SourceFile     : float
            [<FieldAttribute(10)>]
            Quant_SourceFile             : float
            [<FieldAttribute(11)>][<TraceConverter>]
            RtTrace_SourceFile           : float []
            [<FieldAttribute(12)>][<TraceConverter>]
            IntensityTrace_SourceFile    : float []
            [<FieldAttribute(13)>][<TraceConverter>]
            IsotopicPatternMz_SourceFile                    : float []            
            [<FieldAttribute(14)>][<TraceConverter>]
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

    ///
    type PeptideEvidenceClassConverter() = 
        inherit ConverterAttribute()
        override this.convertToObj = 
            Converter.Single(fun (str : string) -> 
                match str with 
                | "Unknown" -> PeptideClassification.PeptideEvidenceClass.Unknown 
                | "C1a"     -> PeptideClassification.PeptideEvidenceClass.C1a 
                | "C1b"     -> PeptideClassification.PeptideEvidenceClass.C1b 
                | "C2a"     -> PeptideClassification.PeptideEvidenceClass.C2a 
                | "C2b"     -> PeptideClassification.PeptideEvidenceClass.C2b 
                | "C3a"     -> PeptideClassification.PeptideEvidenceClass.C3a 
                | "C3b"     -> PeptideClassification.PeptideEvidenceClass.C3b 
                | _ -> 
                    printfn "%s: PeptideEvidenceClass could not be parsed" str
                    PeptideClassification.PeptideEvidenceClass.Unknown 
                |> box  
                )
          
    /// For a group of proteins, contains information about all peptides that are put into the output file.
    type ProteinInferenceResult =
        {
            [<FieldAttribute(0)>]
            ProteinGroup: string
            [<FieldAttribute(1)>]
            PeptideSequence  : string
            [<FieldAttribute(2)>] [<PeptideEvidenceClassConverter>]
            Class            : PeptideClassification.PeptideEvidenceClass
            [<FieldAttribute(3)>]
            TargetScore      : float
            [<FieldAttribute(4)>]
            DecoyScore       : float
            [<FieldAttribute(5)>]
            QValue           : float
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
                TryGetProteinIdentifier = tryParseProteinIdUsing dtoProteinInferenceParams.ProteinIdentifierRegex
                Protein                = dtoProteinInferenceParams.Protein
                Peptide                = dtoProteinInferenceParams.Peptide
                GroupFiles             = dtoProteinInferenceParams.GroupFiles
                GetQValue              = dtoProteinInferenceParams.GetQValue
            }

    type ProteinAssignedQuantifiedIon = { 
        [<FieldAttribute(0)>]
        FileName                                    : string
        [<FieldAttribute(1)>]
        ProteinGroup                                : string
        [<FieldAttribute(2)>]
        ProteinGroup_QValue                         : float
        [<FieldAttribute(3)>]
        StringSequence                              : string
        [<FieldAttribute(4)>]
        PepSequenceID                               : int
        [<FieldAttribute(5)>]
        ModSequenceID                               : int
        [<FieldAttribute(6)>]
        Charge                                      : int
        [<FieldAttribute(7)>]
        GlobalMod                                   : int                    
        [<FieldAttribute(8)>]
        PrecursorMZ                                 : float
        [<FieldAttribute(9)>]
        MeasuredMass                                : float 
        [<FieldAttribute(10)>]
        TheoMass                                    : float
        [<FieldAttribute(11)>]
        AbsDeltaMass                                : float
        [<FieldAttribute(12)>]
        MeanPercolatorScore                         : float
        [<FieldAttribute(13)>]
        QValue                                      : float
        [<FieldAttribute(14)>]
        PEPValue                                    : float
        [<FieldAttribute(15)>]
        ProteinNames                                : string
        [<FieldAttribute(16)>]
        QuantMz_Light                               : float
        [<FieldAttribute(17)>]
        Quant_Light                                 : float
        [<FieldAttribute(18)>]
        MeasuredApex_Light                          : float 
        [<FieldAttribute(19)>]
        Seo_Light                                   : float
        [<FieldAttribute(20)>][<TraceConverter>]
        Params_Light                                : float []
        [<FieldAttribute(21)>]
        Difference_SearchRT_FittedRT_Light          : float
        [<FieldAttribute(22)>]
        KLDiv_Observed_Theoretical_Light            : float
        [<FieldAttribute(23)>]
        KLDiv_CorrectedObserved_Theoretical_Light   : float
        [<FieldAttribute(24)>]
        QuantMz_Heavy                               : float
        [<FieldAttribute(25)>]
        Quant_Heavy                                 : float
        [<FieldAttribute(26)>]
        MeasuredApex_Heavy                          : float
        [<FieldAttribute(27)>]
        Seo_Heavy                                   : float
        [<FieldAttribute(28)>][<TraceConverter>]
        Params_Heavy                                : float []        
        [<FieldAttribute(29)>]
        Difference_SearchRT_FittedRT_Heavy          : float
        [<FieldAttribute(30)>]
        KLDiv_Observed_Theoretical_Heavy            : float
        [<FieldAttribute(31)>]
        KLDiv_CorrectedObserved_Theoretical_Heavy   : float
        [<FieldAttribute(32)>]
        Correlation_Light_Heavy                     : float
        [<FieldAttribute(33)>][<QuantSourceConverter>]
        QuantificationSource                        : QuantificationSource
        [<FieldAttribute(34)>][<TraceConverter>]
        IsotopicPatternMz_Light                     : float []
        [<FieldAttribute(35)>][<TraceConverter>]
        IsotopicPatternIntensity_Observed_Light     : float []
        [<FieldAttribute(36)>][<TraceConverter>]
        IsotopicPatternIntensity_Corrected_Light    : float []
        [<FieldAttribute(37)>][<TraceConverter>]
        RtTrace_Light                               : float []
        [<FieldAttribute(38)>][<TraceConverter>]
        IntensityTrace_Observed_Light               : float []
        [<FieldAttribute(39)>][<TraceConverter>]
        IntensityTrace_Corrected_Light              : float []
        [<FieldAttribute(40)>][<TraceConverter>]
        IsotopicPatternMz_Heavy                     : float []
        [<FieldAttribute(41)>][<TraceConverter>]
        IsotopicPatternIntensity_Observed_Heavy     : float []
        [<FieldAttribute(42)>][<TraceConverter>]
        IsotopicPatternIntensity_Corrected_Heavy    : float []
        [<FieldAttribute(43)>][<TraceConverter>]
        RtTrace_Heavy                               : float []
        [<FieldAttribute(44)>][<TraceConverter>]
        IntensityTrace_Observed_Heavy               : float []
        [<FieldAttribute(45)>][<TraceConverter>]
        IntensityTrace_Corrected_Heavy              : float []
        [<FieldAttribute(46)>]
        AlignmentScore                              : float
        [<FieldAttribute(47)>]
        AlignmentQValue                             : float
        }

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
            StatisticalMeasurements     : (string*Domain.StatisticalMeasurement)[]
            AggregatorFunction          : Domain.AggregationMethod
            AggregatorFunctionIntensity : Domain.AggregationMethod
            AggregatorPepToProt         : Domain.AggregationMethod
            Tukey                       : (string*float*Domain.Transform) []
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
                DistinctPeptideCount        = dtoTableSortParams.DistinctPeptideCount
                StatisticalMeasurements     = dtoTableSortParams.StatisticalMeasurements
                AggregatorFunction          = dtoTableSortParams.AggregatorFunction
                AggregatorFunctionIntensity = dtoTableSortParams.AggregatorFunctionIntensity
                AggregatorPepToProt         = dtoTableSortParams.AggregatorPepToProt
                Tukey                       = dtoTableSortParams.Tukey
            }
                       


    type LabeledQuantificationParams = 
        {
        Correlation_Light_Heavy_Threshold    : float option
        Alignment_QValue                     : float option
        ModificationFilter                   : UseModifiedPeptides 
        AggregateGlobalModificationsParams   : LabeledProteinQuantification.AggregationParams
        AggregatePeptideChargeStatesParams   : LabeledProteinQuantification.AggregationParams option
        AggregateModifiedPeptidesParams      : LabeledProteinQuantification.AggregationParams option
        AggregateToProteinGroupsParams       : LabeledProteinQuantification.AggregationParams
        }   
       
    module LabeledQuantificationParams =
    
        let inline toDomain (dtoTableSortParams: LabeledQuantificationParams): Domain.LabeledQuantificationParams = 
            {
                Correlation_Light_Heavy_Threshold       = dtoTableSortParams.Correlation_Light_Heavy_Threshold
                Alignment_QValue                        = dtoTableSortParams.Alignment_QValue
                ModificationFilter                      = dtoTableSortParams.ModificationFilter |> UseModifiedPeptides.toDomain              
                AggregateGlobalModificationsParams      = dtoTableSortParams.AggregateGlobalModificationsParams |> LabeledProteinQuantification.AggregationParams.toDomain
                AggregatePeptideChargeStatesParams      = dtoTableSortParams.AggregatePeptideChargeStatesParams |> Option.map LabeledProteinQuantification.AggregationParams.toDomain
                AggregateModifiedPeptidesParams         = dtoTableSortParams.AggregateModifiedPeptidesParams    |> Option.map LabeledProteinQuantification.AggregationParams.toDomain
                AggregateToProteinGroupsParams          = dtoTableSortParams.AggregateToProteinGroupsParams     |> LabeledProteinQuantification.AggregationParams.toDomain
            }
    
    type LabelFreeQuantificationParams = 
        {
        Alignment_QValue                     : float option
        ModificationFilter                   : UseModifiedPeptides 
        AggregatePeptideChargeStatesParams   : LabelFreeQuantification.AggregationParams option
        AggregateModifiedPeptidesParams      : LabelFreeQuantification.AggregationParams option
        AggregateToProteinGroupsParams       : LabelFreeQuantification.AggregationParams
        }   
       
    module LabelFreeQuantificationParams =
    
        let inline toDomain (dtoTableSortParams: LabelFreeQuantificationParams): Domain.LabelFreeQuantificationParams = 
            {
                Alignment_QValue                        = dtoTableSortParams.Alignment_QValue
                ModificationFilter                      = dtoTableSortParams.ModificationFilter |> UseModifiedPeptides.toDomain              
                AggregatePeptideChargeStatesParams      = dtoTableSortParams.AggregatePeptideChargeStatesParams |> Option.map Common.LabelFreeQuantification.AggregationParams.toDomain 
                AggregateModifiedPeptidesParams         = dtoTableSortParams.AggregateModifiedPeptidesParams    |> Option.map Common.LabelFreeQuantification.AggregationParams.toDomain
                AggregateToProteinGroupsParams          = dtoTableSortParams.AggregateToProteinGroupsParams     |> Common.LabelFreeQuantification.AggregationParams.toDomain 
            }

    type AlignmentBasedQuantStatisticsParams =
        {
            PositiveQuantMzCutoff: float
            NegativeQuantMzCutoff: float
            PositiveQuantCutoff: float
            NegativeQuantCutoff: float
        }
        
    module AlignmentBasedQuantStatisticsParams =
        
        let inline toDomain (dtoAlignmentBasedQuantStatisticsParams) : Domain.AlignmentBasedQuantStatisticsParams =
            {
                PositiveQuantMzCutoff = dtoAlignmentBasedQuantStatisticsParams.PositiveQuantMzCutoff
                NegativeQuantMzCutoff = dtoAlignmentBasedQuantStatisticsParams.NegativeQuantMzCutoff
                PositiveQuantCutoff = dtoAlignmentBasedQuantStatisticsParams.PositiveQuantCutoff
                NegativeQuantCutoff = dtoAlignmentBasedQuantStatisticsParams.NegativeQuantCutoff
            }
            
    type SpectralLibraryParams =
        {
            MatchingTolerancePPM: float
        }

    module SpectralLibraryParams =

        let toDomain (dtoSpectralLibraryParams: SpectralLibraryParams): Domain.SpectralLibraryParams =
            {
                MatchingTolerancePPM = dtoSpectralLibraryParams.MatchingTolerancePPM
            }

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
        ConsensusAlignmentAlgorithm                 : ConsensusAlignmentAlgorithm
        }

    module ConsensusSpectralLibraryParams =

        let toDomain (dtoConsensusSpectralLibraryParams: ConsensusSpectralLibraryParams): Domain.ConsensusSpectralLibraryParams =
            {
                // InitialPeptideSelection
                BinningWindowWidth                          = dtoConsensusSpectralLibraryParams.BinningWindowWidth                          
                FractionOfMostAbundandIonsPerBin            = dtoConsensusSpectralLibraryParams.FractionOfMostAbundandIonsPerBin            
                MinFragmentCount                            = dtoConsensusSpectralLibraryParams.MinFragmentCount                            
                MinFragmentLadderIdx                        = dtoConsensusSpectralLibraryParams.MinFragmentLadderIdx                        
                MinPeptideLength                            = dtoConsensusSpectralLibraryParams.MinPeptideLength                            
                // XicExtraction                            = // XicExtraction
                RtWindowWidth                               = dtoConsensusSpectralLibraryParams.RtWindowWidth                               
                // Matching                                 = // Matching
                FragMatchingBinWidth                        = dtoConsensusSpectralLibraryParams.FragMatchingBinWidth                        
                FragMatchingBinOffset                       = dtoConsensusSpectralLibraryParams.FragMatchingBinOffset                       
                MS2ScanRange                                = dtoConsensusSpectralLibraryParams.MS2ScanRange                                
                // Filtering                                = // Filtering
                MaxRatioMostAbundandVsSecondAbundandPeak    = dtoConsensusSpectralLibraryParams.MaxRatioMostAbundandVsSecondAbundandPeak    
                ConsensusAlignmentAlgorithm                 = dtoConsensusSpectralLibraryParams.ConsensusAlignmentAlgorithm
            }

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

    module SWATHAnalysisParams =

        let toDomain (dtoSWATHAnalysisParams: SWATHAnalysisParams): Domain.SWATHAnalysisParams =
            {
            
            BinningWindowWidth                         = dtoSWATHAnalysisParams.BinningWindowWidth                       
            FractionOfMostAbundandIonsPerBin           = dtoSWATHAnalysisParams.FractionOfMostAbundandIonsPerBin         
            MinFragmentCount                           = dtoSWATHAnalysisParams.MinFragmentCount                         
            MinFragmentLadderIdx                       = dtoSWATHAnalysisParams.MinFragmentLadderIdx                     
            MinPeptideLength                           = dtoSWATHAnalysisParams.MinPeptideLength                                                          
            RtWindowWidth                              = dtoSWATHAnalysisParams.RtWindowWidth                                                                 
            FragMatchingBinWidth                       = dtoSWATHAnalysisParams.FragMatchingBinWidth                     
            FragMatchingBinOffset                      = dtoSWATHAnalysisParams.FragMatchingBinOffset                    
            MS2ScanRange                               = dtoSWATHAnalysisParams.MS2ScanRange                                                                
            MaxRatioMostAbundandVsSecondAbundandPeak   = dtoSWATHAnalysisParams.MaxRatioMostAbundandVsSecondAbundandPeak 
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

    module MzTABParams =

        let toDomain (dtoMzTABParams: MzTABParams): Domain.MzTABParams =
            {
                ExperimentNames       = dtoMzTABParams.ExperimentNames
                StudyVariables        = dtoMzTABParams.StudyVariables
                SearchEngineNamesProt = dtoMzTABParams.SearchEngineNamesProt
                SearchEngineNamesPep  = dtoMzTABParams.SearchEngineNamesPep
                SearchEngineNamesPSM  = dtoMzTABParams.SearchEngineNamesPSM
                Labeled               = dtoMzTABParams.Labeled
                MetaData              = dtoMzTABParams.MetaData
                FieldNames            = dtoMzTABParams.FieldNames
            }

    type FragmentIon = {
        Charge                            : float
        Iontype                           : string
        Number                            : int
        // Brauch man glaub ich nicht weil das in der Quantinfo steckt
        GlobalMod                         : int
        CountAbsolute                     : int
        CountFraction                     : float
        CalculatedMz                      : float
        MeanFragMz                        : float 
        MeanAbsMzDelta                    : float
        MeanTheoMinusXMzDelta             : float
        CvMeanFragMz                      : float
        MaxIntensity                      : float
        MinIntensity                      : float
        MeanIntensity                     : float
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

    type MzliteToMzMLParams =
        {
            Compress                    : Compression
            StartRetentionTime          : float option
            EndRetentionTime            : float option
        }

    module MzliteToMzMLParams =

        let toDomain (dtoMzMLConverterParams: MzliteToMzMLParams ) : Domain.MzliteToMzMLParams =
                {
                    Compress                    = Compression.toDomain dtoMzMLConverterParams.Compress
                    StartRetentionTime          = dtoMzMLConverterParams.StartRetentionTime
                    EndRetentionTime            = dtoMzMLConverterParams.EndRetentionTime
                }

    type MzMLtoMzLiteParams = PreprocessingParams

    module MzMLtoMzLiteParams =

        let toDomain (dtoCentroidizationParams: MzMLtoMzLiteParams ) : Domain.MzMLtoMzLiteParams =
                {
                    Compress                    = Compression.toDomain dtoCentroidizationParams.Compress
                    StartRetentionTime          = dtoCentroidizationParams.StartRetentionTime
                    EndRetentionTime            = dtoCentroidizationParams.EndRetentionTime
                    MS1PeakPicking              = dtoCentroidizationParams.MS1PeakPicking
                    MS2PeakPicking              = dtoCentroidizationParams.MS2PeakPicking
                }
 