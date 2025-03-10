// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\MzMLToMzLite\net5.0\Newtonsoft.Json.dll"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"


open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain

    
let defaultLabelFreeQuantificationParams :Dto.LabelFreeQuantificationParams = 
    let chargeParams :Common.LabelFreeQuantification.AggregationParams = 
        {
                Transform        = None
                SingleFilters    = None
                GroupFilters     = None
                Aggregation      = {Light= NumericAggregation.Mean}
        }
    let modParams :Common.LabelFreeQuantification.AggregationParams = 
        {
                Transform        = None
                SingleFilters    = None
                GroupFilters     = None
                Aggregation      = {Light= NumericAggregation.Mean}
        }
    let protParams :Common.LabelFreeQuantification.AggregationParams = 
        {
                Transform        = None
                SingleFilters    = None
                GroupFilters     = None
                Aggregation      = {Light= NumericAggregation.Mean}
        }
    {
        ModificationFilter                   = UseModifiedPeptides.All 
        AggregatePeptideChargeStatesParams   = Some modParams
        AggregateModifiedPeptidesParams      = Some modParams
        AggregateToProteinGroupsParams       = protParams
    }   


let serialized = 
    defaultLabelFreeQuantificationParams
    |> Json.serialize


System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\LabelFreeQuantificationParams_ChargeAgg_ModAgg.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\LabelFreeQuantificationParams_ChargeAgg_ModAgg.json")
    |> Json.deserialize<Dto.LabelFreeQuantificationParams>
    |> LabelFreeQuantificationParams.toDomain

    