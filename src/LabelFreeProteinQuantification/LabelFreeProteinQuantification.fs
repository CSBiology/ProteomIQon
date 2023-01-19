namespace ProteomIQon

open FSharpAux
open Deedle
open FSharpAux.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute
open Dto
open FSharp.Stats
open System.IO
open ProteomIQon.Drafo

module LabelFreeProteinQuantification = 
    
    open Domain.LabelFreeQuantification

    let tryGetQVal s =
        match s |> Seq.filter (nan.Equals >> not) |> Seq.tryItem 0 with 
        | Some x -> x 
        | None -> nan 

    let performChargeOrModAggregation (transformParams:Transforms) (singleFilters:SingleFilters) (groupFilters:GroupFilters) (aggregations:Aggregations) alignment_QValue modifyKeyColumns keyColumns peptidesAndProteinsIndexed =
            let light       = peptidesAndProteinsIndexed |> Core.getColumn<float>"Quant_Light"
            let proteinQVal = peptidesAndProteinsIndexed |> Core.getColumn<float>"ProteinGroup_QValue"
            let alignmentQVal = peptidesAndProteinsIndexed |> Core.getColumn<float>"AlignmentQValue"
            /// Transformed 
            let lightT       = light    |> Core.transform transformParams.Light
            /// ThresholdFilters
            let lightTF = singleFilters.Light |> Seq.map (fun f -> Core.createFilter f lightT)
            /// Alignment QValue Filter
            let alignmentQvalF = 
                match alignment_QValue with 
                | Some a -> 
                    Core.createFilter (fun x -> 
                        x < a || isNan x
                    ) alignmentQVal
                | None -> Core.createFilter (fun x -> true ) alignmentQVal
            /// GroupFilters = 
            let lightGF = groupFilters.Light |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns lightT)
            /// CombinedFilters 
            let cf = alignmentQvalF::([lightTF;lightGF] |> Seq.concat |> List.ofSeq) 
            /// Aggregated
            let light_agg :Series<_,float> = Core.aggregate aggregations.Light modifyKeyColumns keyColumns cf lightT 
            let proteinQVal :Series<_,float> = Core.aggregate tryGetQVal modifyKeyColumns keyColumns cf proteinQVal 
            /// Assembeled
            let assembeled = 
                Core.assemble 
                    [
                    "Quant_Light", light_agg :> ISeries<Core.Key>
                    "ProteinGroup_QValue", proteinQVal :> ISeries<Core.Key>
                    ]
            assembeled

    let performPeptideToProteinAggregation (transformParams:Transforms) (singleFilters:SingleFilters) (groupFilters:GroupFilters) (aggregations:Aggregations) alignment_QValue modifyKeyColumns keyColumns peptidesAndProteinsIndexed =
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let light       = peptidesAndProteinsIndexed |> Core.getColumn<float>"Quant_Light"
            let proteinQVal = peptidesAndProteinsIndexed |> Core.getColumn<float>"ProteinGroup_QValue"
            let alignmentQVal = peptidesAndProteinsIndexed |> Core.getColumn<float>"AlignmentQValue"
            /// Transformed 
            let lightT       = light    |> Core.transform transformParams.Light
            /// 
            /// ThresholdFilters
            let lightTF = singleFilters.Light |> Seq.map (fun f -> Core.createFilter f lightT)
            /// Alignment QValue Filter
            let alignmentQvalF = 
                match alignment_QValue with 
                | Some a -> 
                    Core.createFilter (fun x -> 
                        x < a || isNan x
                    ) alignmentQVal
                | None -> Core.createFilter (fun x -> true ) alignmentQVal
            /// GroupFilters = 
            let lightGF = groupFilters.Light |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns lightT)
            /// CombinedFilters 
            let cf = alignmentQvalF::([lightTF;lightGF] |> Seq.concat |> List.ofSeq) 
            /// Aggregated
            /// Light
            let light_agg :Series<_,float> = Core.aggregate aggregations.Light modifyKeyColumns keyColumns cf lightT 
            let nLight :Series<_,int> = Core.aggregate (Seq.length) modifyKeyColumns keyColumns cf lightT       
            let light_cv :Series<_,float> = Core.aggregate Seq.cv modifyKeyColumns keyColumns cf lightT 
            let light_stDev :Series<_,float> = Core.aggregate (Seq.stDevBy float) modifyKeyColumns keyColumns cf lightT 
            let light_SEM :Series<_,float> = Core.aggregate sem modifyKeyColumns keyColumns cf lightT                   
            let proteinQVal :Series<_,float> = Core.aggregate tryGetQVal modifyKeyColumns keyColumns cf proteinQVal 
            /// Assembeled
            let assembeled = 
                Core.assemble 
                    [
                    "ProteinGroup_QValue", proteinQVal :> ISeries<Core.Key>                 
                    "ItemsUsedForQuant_Light", nLight :> ISeries<Core.Key>
                    "Quant_Light", light_agg :> ISeries<Core.Key>
                    "CV_Quant_Light", light_cv :> ISeries<Core.Key>
                    "StDev_Quant_Light", light_stDev :> ISeries<Core.Key>
                    "SEM_Quant_Light", light_SEM :> ISeries<Core.Key>
                    ]
            assembeled
            
    let labelFreeQuantification (labelfreeQuantificationParams:Domain.LabelFreeQuantificationParams) (outputDir:string) (instrumentOutput:string[]) =
        let logger = Logging.createLogger "LabeledQuantification"
        logger.Trace (sprintf "Input files: %A" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        let modPepFilter (p:ProteinAssignedQuantifiedIon) = p.StringSequence |> labelfreeQuantificationParams.ModificationFilter
        let keyCols =         
            [|
                "FileName"      
                "ProteinGroup"  
                "StringSequence"
                "PepSequenceID" 
                "ModSequenceID" 
                "Charge"        
                "GlobalMod"     
            |]
        let peptidesAndProteinsIndexed = 
            instrumentOutput
            |> Array.map (fun fp ->
                    let peptidesAndProteinsIndexed' =
                        Csv.CsvReader<ProteinAssignedQuantifiedIon>(SchemaMode=Csv.Fill).ReadFile(fp,'\t',false,1)
                        |> Array.ofSeq
                        |> Array.filter (fun x -> x.GlobalMod = 0)
                        |> Array.filter modPepFilter
                        |> Frame.ofRecords
                        |> Core.indexWithColumnValues keyCols
                        |> Frame.sliceCols
                            [|
                                "Quant_Light"
                                "ProteinGroup_QValue"
                                "AlignmentQValue"
                            |]
                    logger.Trace (sprintf "QuantAndProt file with name:%s contributes %i quantifications" (System.IO.Path.GetFileNameWithoutExtension fp) (peptidesAndProteinsIndexed'.RowCount))
                    peptidesAndProteinsIndexed'
                )
            |> Frame.mergeAll
        logger.Trace (sprintf "Starting Aggregation with %i quantifications" (peptidesAndProteinsIndexed.RowCount))       
      /// Step 1: Aggregate Charges
        ///  
        let chargesAggregated = 
            match labelfreeQuantificationParams.AggregatePeptideChargeStatesParams with 
            | Some chParams ->
                let res = 
                    performChargeOrModAggregation 
                        chParams.Transform chParams.SingleFilters chParams.GroupFilters chParams.Aggregation labelfreeQuantificationParams.Alignment_QValue
                            Core.dropKeyColumns ["Charge";"GlobalMod"] peptidesAndProteinsIndexed
                logger.Trace (sprintf "Charge State to modified peptide sequence aggregation yields %i quantifications." (res.RowCount))
                res  
            | None ->  
                logger.Trace (sprintf "Charge state aggregation is skipped.") 
                peptidesAndProteinsIndexed       
        /// Step 2: Aggregate Modifications
        let modificationsAggregated = 
            match labelfreeQuantificationParams.AggregateModifiedPeptidesParams with
            | Some modParams ->
                let res = 
                    performChargeOrModAggregation 
                        modParams.Transform modParams.SingleFilters modParams.GroupFilters modParams.Aggregation labelfreeQuantificationParams.Alignment_QValue
                            Core.dropKeyColumns ["StringSequence";"ModSequenceID"] chargesAggregated
                logger.Trace (sprintf "Modified peptide sequence to peptide Sequence aggregation yields %i quantifications." (res.RowCount)) 
                res
            | None -> 
                logger.Trace (sprintf "Modified peptide sequence to peptide Sequence aggregation is skipped.")
                chargesAggregated
        /// Step 3: Aggregate To Proteins
        let proteins = 
            let pParams = labelfreeQuantificationParams.AggregateToProteinGroupsParams
            let res = 
                performPeptideToProteinAggregation
                    pParams.Transform pParams.SingleFilters pParams.GroupFilters pParams.Aggregation labelfreeQuantificationParams.Alignment_QValue
                                 Core.dropAllKeyColumnsBut ["FileName";"ProteinGroup"] modificationsAggregated
            logger.Trace (sprintf "Peptide to Protein Group aggregation yields %i quantifications." (res.RowCount))                
            res
        /// 
        let proteinsPivoted = 
           proteins
           |> Core.pivot "FileName"
           |> Core.rowKeyToColumns 
        logger.Trace (sprintf "Finally, %i Protein Groups are reported across all files. %i of all protein groups were quantified across all samples." (proteinsPivoted.RowCount) (proteinsPivoted.DropSparseRows().RowCount))           
        
        if labelfreeQuantificationParams.AggregatePeptideChargeStatesParams.IsSome then 
            let outFilePathC =
               let fileName = "ChargeAggregation.txt"
               Path.Combine [|outputDir;fileName|]
            (chargesAggregated |> Core.rowKeyToColumns).SaveCsv(outFilePathC,separator='\t',includeRowKeys=false)
        if labelfreeQuantificationParams.AggregateModifiedPeptidesParams.IsSome then 
            let outFilePathM =
               let fileName = "ModificationAggregation.txt"
               Path.Combine [|outputDir;fileName|]
            (modificationsAggregated |> Core.rowKeyToColumns).SaveCsv(outFilePathM,separator='\t',includeRowKeys=false)        
        let outFilePathP =
           let fileName = "ProteinAggregation.txt"
           Path.Combine [|outputDir;fileName|]
        (proteins |> Core.rowKeyToColumns).SaveCsv(outFilePathP,separator='\t',includeRowKeys=false)
        let outFilePath =
           let fileName = "LabelFreeQuant.txt"
           Path.Combine [|outputDir;fileName|]
        proteinsPivoted.SaveCsv(outFilePath,separator='\t',includeRowKeys=false)


































