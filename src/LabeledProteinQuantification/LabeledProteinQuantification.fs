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

module LabeledProteinQuantification = 
    
    open Domain.LabeledProteinQuantification

    let tryGetQVal s =
        match s |> Seq.filter (nan.Equals >> not) |> Seq.tryItem 0 with 
        | Some x -> x 
        | None -> nan 

    let performGlobalModAggregation (transformParams:LabeledTransforms) (singleFilters:LabeledSingleFilters) (groupFilters:LabeledGroupFilters) (aggregations:LabeledAggregations) correlation_Light_Heavy_Threshold modifyKeyColumns keyColumns peptidesAndProteinsIndexed =
            let light       = peptidesAndProteinsIndexed |> Core.getColumn<float>"Quant_Light"
            let heavy       = peptidesAndProteinsIndexed |> Core.getColumn<float>"Quant_Heavy"
            let correlation = peptidesAndProteinsIndexed |> Core.getColumn<float>"Correlation_Light_Heavy"
            let proteinQVal = peptidesAndProteinsIndexed |> Core.getColumn<float>"ProteinGroup_QValue"
            /// Zipped
            let n14ByN15 = Core.zip (fun x y -> x / y) light heavy 
            /// Transformed 
            let lightT       = light    |> Core.transform transformParams.Light
            let heavyT       = heavy    |> Core.transform transformParams.Heavy
            let n14ByN15T    = n14ByN15 |> Core.transform transformParams.Ratio
            /// 
            /// ThresholdFilters
            let lightTF = singleFilters.Light |> Seq.map (fun f -> Core.createFilter f lightT)
            let heavyTF = singleFilters.Heavy |> Seq.map (fun f -> Core.createFilter f heavyT)
            let n14ByN15TF = singleFilters.Ratio |> Seq.map (fun f -> Core.createFilter f n14ByN15T)
            let correlationF = 
                match correlation_Light_Heavy_Threshold with 
                | Some c -> Core.createFilter (fun x -> x > c ) correlation
                | None -> Core.createFilter (fun x -> true ) correlation
            /// GroupFilters = 
            let lightGF = groupFilters.Light |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns lightT)
            let heavyGF = groupFilters.Heavy |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns heavyT) 
            let n14ByN15GF = groupFilters.Ratio |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns n14ByN15T)       
            /// CombinedFilters 
            let cf = correlationF::([lightTF;heavyTF;n14ByN15TF;lightGF;heavyGF;n14ByN15GF] |> Seq.concat |> List.ofSeq) 
            /// Aggregated
            let light_agg :Series<_,float> = Core.aggregate aggregations.Light modifyKeyColumns keyColumns cf lightT 
            let heavy_agg :Series<_,float> = Core.aggregate aggregations.Heavy modifyKeyColumns keyColumns cf heavyT 
            let ratio_agg :Series<_,float> = Core.aggregate aggregations.Ratio modifyKeyColumns keyColumns cf n14ByN15T
            let proteinQVal :Series<_,float> = Core.aggregate tryGetQVal modifyKeyColumns keyColumns cf proteinQVal  
            /// Assembeled
            let assembeled = 
                Core.assemble 
                    [
                    "Quant_Light", light_agg :> ISeries<Core.Key>
                    "Quant_Heavy", heavy_agg :> ISeries<Core.Key>
                    "Ratio_LightByHeavy", ratio_agg :> ISeries<Core.Key>
                    "ProteinGroup_QValue", proteinQVal :> ISeries<Core.Key>
                    ]
            assembeled

    let performChargeOrModAggregation (transformParams:LabeledTransforms) (singleFilters:LabeledSingleFilters) (groupFilters:LabeledGroupFilters) (aggregations:LabeledAggregations) modifyKeyColumns keyColumns peptidesAndProteinsIndexed =
            let light       = peptidesAndProteinsIndexed |> Core.getColumn<float>"Quant_Light"
            let heavy       = peptidesAndProteinsIndexed |> Core.getColumn<float>"Quant_Heavy"
            let n14ByN15    = peptidesAndProteinsIndexed |> Core.getColumn<float>"Ratio_LightByHeavy"
            let proteinQVal = peptidesAndProteinsIndexed |> Core.getColumn<float>"ProteinGroup_QValue"
            /// Transformed 
            let lightT       = light    |> Core.transform transformParams.Light
            let heavyT       = heavy    |> Core.transform transformParams.Heavy
            /// Zipped
            let n14ByN15T    = n14ByN15 |> Core.transform transformParams.Ratio
            /// 
            /// ThresholdFilters
            let lightTF = singleFilters.Light |> Seq.map (fun f -> Core.createFilter f lightT)
            let heavyTF = singleFilters.Heavy |> Seq.map (fun f -> Core.createFilter f heavyT)
            let n14ByN15TF = singleFilters.Ratio |> Seq.map (fun f -> Core.createFilter f n14ByN15T)
            /// GroupFilters = 
            let lightGF = groupFilters.Light |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns lightT)
            let heavyGF = groupFilters.Heavy |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns heavyT) 
            let n14ByN15GF = groupFilters.Ratio |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns n14ByN15T)       
            /// CombinedFilters 
            let cf = ([lightTF;heavyTF;n14ByN15TF;lightGF;heavyGF;n14ByN15GF] |> Seq.concat |> List.ofSeq) 
            /// Aggregated
            let light_agg :Series<_,float> = Core.aggregate aggregations.Light modifyKeyColumns keyColumns cf lightT 
            let heavy_agg :Series<_,float> = Core.aggregate aggregations.Heavy modifyKeyColumns keyColumns cf heavyT 
            let ratio_agg :Series<_,float> = Core.aggregate aggregations.Ratio modifyKeyColumns keyColumns cf n14ByN15T 
            let proteinQVal :Series<_,float> = Core.aggregate tryGetQVal modifyKeyColumns keyColumns cf proteinQVal 
            /// Assembeled
            let assembeled = 
                Core.assemble 
                    [
                    "Quant_Light", light_agg :> ISeries<Core.Key>
                    "Quant_Heavy", heavy_agg :> ISeries<Core.Key>
                    "Ratio_LightByHeavy", ratio_agg :> ISeries<Core.Key>
                    "ProteinGroup_QValue", proteinQVal :> ISeries<Core.Key>
                    ]
            assembeled

    let performPeptideToProteinAggregation (transformParams:LabeledTransforms) (singleFilters:LabeledSingleFilters) (groupFilters:LabeledGroupFilters) (aggregations:LabeledAggregations) modifyKeyColumns keyColumns peptidesAndProteinsIndexed =
            let sem (x:seq<float>) =
                let stDev = Seq.stDev x
                stDev / (sqrt (x |> Seq.length |> float))
            let light       = peptidesAndProteinsIndexed |> Core.getColumn<float>"Quant_Light"
            let heavy       = peptidesAndProteinsIndexed |> Core.getColumn<float>"Quant_Heavy"
            let n14ByN15    = peptidesAndProteinsIndexed |> Core.getColumn<float>"Ratio_LightByHeavy"
            let proteinQVal = peptidesAndProteinsIndexed |> Core.getColumn<float>"ProteinGroup_QValue"
            /// Zipped
            /// Transformed 
            let lightT       = light    |> Core.transform transformParams.Light
            let heavyT       = heavy    |> Core.transform transformParams.Heavy
            let n14ByN15T    = n14ByN15 |> Core.transform transformParams.Ratio           
            /// 
            /// ThresholdFilters
            let lightTF = singleFilters.Light |> Seq.map (fun f -> Core.createFilter f lightT)
            let heavyTF = singleFilters.Heavy |> Seq.map (fun f -> Core.createFilter f heavyT)
            let n14ByN15TF = singleFilters.Ratio |> Seq.map (fun f -> Core.createFilter f n14ByN15T)
            /// GroupFilters = 
            let lightGF = groupFilters.Light |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns lightT)
            let heavyGF = groupFilters.Heavy |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns heavyT) 
            let n14ByN15GF = groupFilters.Ratio |> Seq.map (fun f -> Core.createGroupFilter f modifyKeyColumns keyColumns n14ByN15T)       
            /// CombinedFilters 
            let cf = ([lightTF;heavyTF;n14ByN15TF;lightGF;heavyGF;n14ByN15GF] |> Seq.concat |> List.ofSeq) 
            /// Aggregated
            /// Light
            let light_agg :Series<_,float> = Core.aggregate aggregations.Light modifyKeyColumns keyColumns cf lightT 
            let nLight :Series<_,int> = Core.aggregate (Seq.length) modifyKeyColumns keyColumns cf lightT       
            let light_cv :Series<_,float> = Core.aggregate Seq.cv modifyKeyColumns keyColumns cf lightT 
            let light_stDev :Series<_,float> = Core.aggregate (Seq.stDevBy float) modifyKeyColumns keyColumns cf lightT 
            let light_SEM :Series<_,float> = Core.aggregate sem modifyKeyColumns keyColumns cf lightT 
            
            /// Heavy
            let heavy_agg :Series<_,float> = Core.aggregate aggregations.Heavy modifyKeyColumns keyColumns cf heavyT 
            let nHeavy :Series<_,int> = Core.aggregate (Seq.length) modifyKeyColumns keyColumns cf heavyT     
            let heavy_cv :Series<_,float> = Core.aggregate Seq.cv modifyKeyColumns keyColumns cf heavyT 
            let heavy_stDev :Series<_,float> = Core.aggregate (Seq.stDevBy float) modifyKeyColumns keyColumns cf heavyT 
            let heavy_SEM :Series<_,float> = Core.aggregate sem modifyKeyColumns keyColumns cf heavyT 
            
            /// Ratio
            let ratio_agg :Series<_,float> = Core.aggregate aggregations.Ratio modifyKeyColumns keyColumns cf n14ByN15T 
            let nRatio :Series<_,int> = Core.aggregate (Seq.length) modifyKeyColumns keyColumns cf n14ByN15T
            let ratio_cv :Series<_,float> = Core.aggregate Seq.cv modifyKeyColumns keyColumns cf n14ByN15T 
            let ratio_stDev :Series<_,float> = Core.aggregate (Seq.stDevBy float) modifyKeyColumns keyColumns cf n14ByN15T 
            let ratio_SEM :Series<_,float> = Core.aggregate sem modifyKeyColumns keyColumns cf n14ByN15T  
                       
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
                    
                    "ItemsUsedForQuant_Heavy", nHeavy :> ISeries<Core.Key>
                    "Quant_Heavy", heavy_agg :> ISeries<Core.Key>
                    "CV_Quant_Heavy", heavy_cv :> ISeries<Core.Key>
                    "StDev_Quant_Heavy", heavy_stDev :> ISeries<Core.Key>
                    "SEM_Quant_Heavy", heavy_SEM :> ISeries<Core.Key>
                    
                    "ItemsUsedForQuant_LightByHeavy", nRatio :> ISeries<Core.Key>
                    "Ratio_LightByHeavy", ratio_agg :> ISeries<Core.Key>
                    "CV_Quant_Ratio", ratio_cv :> ISeries<Core.Key>
                    "StDev_Quant_Ratio", ratio_stDev :> ISeries<Core.Key>
                    "SEM_Quant_Ratio", ratio_SEM :> ISeries<Core.Key>
                    ]
            assembeled
            
    let labeledQuantification (labeledQuantificationParams:Domain.LabeledQuantificationParams) (outputDir:string) (instrumentOutput:string[]) =
        let logger = Logging.createLogger "LabeledQuantification"
        logger.Trace (sprintf "Input files: %A" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        let modPepFilter (p:ProteinAssignedQuantifiedIon) = p.StringSequence |> labeledQuantificationParams.ModificationFilter
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
                        |> Array.filter modPepFilter
                        |> Frame.ofRecords
                        |> Core.indexWithColumnValues keyCols 
                    logger.Trace (sprintf "QuantAndProt file with name:%s contributes %i quantifications" (System.IO.Path.GetFileNameWithoutExtension fp) (peptidesAndProteinsIndexed'.RowCount))
                    peptidesAndProteinsIndexed'
                )
            |> Frame.mergeAll
        logger.Trace (sprintf "Starting Aggregation with %i quantifications" (peptidesAndProteinsIndexed.RowCount))       
        /// Step 1: Aggregate GlobalModifications
        ///
        let globalModAggregated =
            let gParams = labeledQuantificationParams.AggregateGlobalModificationsParams
            performGlobalModAggregation 
                gParams.LabeledTransform gParams.LabeledSingleFilters gParams.LabeledGroupFilters gParams.LabeledAggregation labeledQuantificationParams.Correlation_Light_Heavy_Threshold 
                    Core.dropKeyColumns ["GlobalMod";"ModSequenceID"] peptidesAndProteinsIndexed
        logger.Trace (sprintf "Global Modification aggregation yields %i quantifications" (globalModAggregated.RowCount))   
        /// Step 2: Aggregate Charges
        ///  
        let chargesAggregated = 
            match labeledQuantificationParams.AggregatePeptideChargeStatesParams with 
            | Some chParams ->
                let res = 
                    performChargeOrModAggregation 
                        chParams.LabeledTransform chParams.LabeledSingleFilters chParams.LabeledGroupFilters chParams.LabeledAggregation 
                            Core.dropKeyColumns ["Charge"] globalModAggregated
                logger.Trace (sprintf "Charge State to modified peptide sequence aggregation yields %i quantifications." (res.RowCount))
                res  
            | None ->  
                logger.Trace (sprintf "Charge state aggregation is skipped.") 
                globalModAggregated       
        /// Step 3: Aggregate Modifications
        let modificationsAggregated = 
            match labeledQuantificationParams.AggregateModifiedPeptidesParams with
            | Some modParams ->
                let res = 
                    performChargeOrModAggregation 
                        modParams.LabeledTransform modParams.LabeledSingleFilters modParams.LabeledGroupFilters modParams.LabeledAggregation 
                            Core.dropKeyColumns ["StringSequence";] chargesAggregated
                logger.Trace (sprintf "Modified peptide sequence to peptide Sequence aggregation yields %i quantifications." (res.RowCount)) 
                res
            | None -> 
                logger.Trace (sprintf "Modified peptide sequence to peptide Sequence aggregation is skipped.")
                chargesAggregated
        /// Step 4: Aggregate To Proteins
        let proteins = 
            let pParams = labeledQuantificationParams.AggregateToProteinGroupsParams
            let res = 
                performPeptideToProteinAggregation
                    pParams.LabeledTransform pParams.LabeledSingleFilters pParams.LabeledGroupFilters pParams.LabeledAggregation 
                                 Core.dropAllKeyColumnsBut ["FileName";"ProteinGroup"] modificationsAggregated
            logger.Trace (sprintf "Peptide to Protein Group aggregation yields %i quantifications." (res.RowCount))                
            res
        /// 
        let proteinsPivoted = 
           proteins
           |> Core.pivot "FileName"
           |> Core.rowKeyToColumns 
        logger.Trace (sprintf "Finally, %i Protein Groups are reported across all files. %i of all protein groups were quantified across all samples." (proteinsPivoted.RowCount) (proteinsPivoted.DropSparseRows().RowCount))           
        
        /// Write files
        let outFilePathG =
           let fileName = "GlobModAggregation.txt"
           Path.Combine [|outputDir;fileName|]
        (globalModAggregated |> Core.rowKeyToColumns).SaveCsv(outFilePathG,separator='\t',includeRowKeys=false)
        if labeledQuantificationParams.AggregatePeptideChargeStatesParams.IsSome then 
            let outFilePathC =
               let fileName = "ChargeAggregation.txt"
               Path.Combine [|outputDir;fileName|]
            (chargesAggregated |> Core.rowKeyToColumns).SaveCsv(outFilePathC,separator='\t',includeRowKeys=false)
        if labeledQuantificationParams.AggregateModifiedPeptidesParams.IsSome then 
            let outFilePathM =
               let fileName = "ModificationAggregation.txt"
               Path.Combine [|outputDir;fileName|]
            (modificationsAggregated |> Core.rowKeyToColumns).SaveCsv(outFilePathM,separator='\t',includeRowKeys=false)        
        let outFilePathP =
           let fileName = "ProteinAggregation.txt"
           Path.Combine [|outputDir;fileName|]
        (proteins |> Core.rowKeyToColumns).SaveCsv(outFilePathP,separator='\t',includeRowKeys=false)
        let outFilePath =
           let fileName = "LabeledQuant.txt"
           Path.Combine [|outputDir;fileName|]
        proteinsPivoted.SaveCsv(outFilePath,separator='\t',includeRowKeys=false)


































