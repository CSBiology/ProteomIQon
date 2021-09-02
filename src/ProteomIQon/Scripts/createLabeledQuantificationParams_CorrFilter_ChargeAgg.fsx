// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\MzMLToMzLite\net5.0\Newtonsoft.Json.dll"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"


open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Dto.LabeledQuantificationParams
open ProteomIQon.Domain

    
let defaultLabeledQuantificationParams :Dto.LabeledQuantificationParams = 

    let globalModParams :Common.LabeledProteinQuantification.AggregationParams = 
        {
                LabeledTransform        = None
                LabeledSingleFilters    = None
                LabeledGroupFilters     = None
                LabeledAggregation      = {Light= NumericAggregation.Mean; Heavy=NumericAggregation.Mean; Ratio=NumericAggregation.Mean} 
        }
    let chargeParams :Common.LabeledProteinQuantification.AggregationParams = 
        {
                LabeledTransform        = None
                LabeledSingleFilters    = None
                LabeledGroupFilters     = None
                LabeledAggregation      = {Light= NumericAggregation.Mean; Heavy=NumericAggregation.Mean; Ratio=NumericAggregation.Mean} 
        }

    let protParams :Common.LabeledProteinQuantification.AggregationParams = 
        {
                LabeledTransform        = None
                LabeledSingleFilters    = None
                LabeledGroupFilters     = None
                LabeledAggregation      = {Light= NumericAggregation.Mean; Heavy=NumericAggregation.Mean; Ratio=NumericAggregation.Mean}
        }
    {
        Correlation_Light_Heavy_Threshold    = Some 0.0
        ModificationFilter                   = UseModifiedPeptides.All 
        AggregateGlobalModificationsParams   = globalModParams
        AggregatePeptideChargeStatesParams   = Some chargeParams
        AggregateModifiedPeptidesParams      = None
        AggregateToProteinGroupsParams       = protParams
    }   


let serialized = 
    defaultLabeledQuantificationParams
    |> Json.serialize


System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\LabeledQuantificationParams_CorrFilter_ChargeAgg.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\LabeledQuantificationParams_CorrFilter_ChargeAgg.json")
    |> Json.deserialize<Dto.LabeledQuantificationParams>
    |> LabeledQuantificationParams.toDomain
