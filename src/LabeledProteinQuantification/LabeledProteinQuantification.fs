namespace ProteomIQon

open FSharpAux
open Deedle
open FSharpAux.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute
open Dto

module Core =
    open Deedle
    open DynamicObj
    open System.Collections.Generic

    type Key() = 
        inherit DynamicObj ()
       
        member this.addCol(columns:(string*'a)) =
            this.SetValue columns
            this  

        override this.ToString() = 
            let sb = new System.Text.StringBuilder()
            // sb.Append
            (this.GetProperties true)        
            |> Seq.iter (fun x -> 
                let value = x.Value.ToString() 
                sb.AppendLine(x.Key+": "+value)
                |> ignore
                )
            (sb.ToString())

        override this.Equals(b) =
            match b with
            | :? Key as p -> 
                let propA = (this.GetProperties true) |> Seq.map (fun x -> x.Key,x.Value) 
                let propB = (p.GetProperties true) |> Seq.map (fun x -> x.Key,x.Value) 
                Seq.map2 (fun (x1,x2) (y1,y2) -> x1 = unbox y1 && x2 = unbox y2) propA propB 
                |> Seq.contains false
                |> not
            | _ -> false

        override this.GetHashCode() = 
            let sb = new System.Text.StringBuilder()
            // sb.Append
            (this.GetProperties true)        
            |> Seq.iter (fun x -> 
                let value = x.Value.ToString() 
                sb.Append(value)
                |> ignore
                )
            (sb.ToString())
            |> hash   

    ///
    let indexWithColumnValues (keyCols:(seq<string>) ) (f:Frame<_,string>) :Frame<Key,_>= 
        f
        |> Frame.indexRowsUsing (fun s -> 
                keyCols
                |> Seq.fold (fun key x -> 
                    let value = s.GetAs<string>(x) 
                    key.addCol (x,value)
                    ) (Key())
            )  

    ///
    let readFrame fp = Frame.ReadCsv(fp,hasHeaders=true,inferTypes=false,separators="\t")
      
    ///
    let readAndIndexFrame keyCols fp = 
        readFrame fp
        |> indexWithColumnValues keyCols

    ///
    let inline getColumn<'a> column (f:Frame<Key, string> )  =
        f.GetColumn<'a>(column)

    ///
    let inline seriesToFrame (s: Series<Key, 'a>) =
        s
        |> Series.map (fun k s -> 
            (k.GetProperties true) 
            |> Seq.map (fun x -> x.Key,x.Value)
            |> Series.ofObservations
        )    
        |> Frame.ofRows
        |> Frame.addCol "Value" s
        |> Frame.indexRowsOrdinally

    ///
    let inline rowKeyToColumns (f: Frame<Key, string>) =
        let rowKeysAsColumns = 
            f
            |> Frame.mapRows (fun k s -> 
                (k.GetProperties true) 
                |> Seq.map (fun x -> x.Key,x.Value)
                |> Series.ofObservations
            )    
            |> Frame.ofRows
        Frame.join JoinKind.Inner rowKeysAsColumns f 
        |> Frame.indexRowsOrdinally

    ///
    let createFilter (op : 'a -> bool) (s: Series<'KeyType, 'a>) = 
        s
        |> Series.mapValues op

    ///
    let transform (op : 'a -> 'b) (s: Series<'KeyType, 'a>) = 
        s
        |> Series.mapValues op

    ///
    let zip (op : 'a -> 'a -> 'b) (s1: Series<'KeyType, 'a>) (s2: Series<'KeyType, 'a>) = 
        Series.zipInner s1 s2
        |> Series.mapValues (fun (s1,s2) -> 
            op s1 s2
            )

    ///
    let dropAllPropertiesBut (properties:seq<string>) (key:Key) = 
        let newK = Key()
        key.GetProperties true
        |> Seq.filter (fun x -> properties |> Seq.contains x.Key )
        |> Seq.fold (fun (k:Key) x -> k.addCol (x.Key,x.Value) ) newK

    ///
    let dropProperties (properties:seq<string>) (key:Key) = 
        let newK = Key()
        key.GetProperties true
        |> Seq.filter (fun x -> properties |> Seq.contains x.Key |> not)
        |> Seq.fold (fun (k:Key) x -> k.addCol (x.Key,x.Value) ) newK

    ///
    let groupTransform (op :'a [] -> 'a -> 'b) (newKeys:seq<string>) (s: Series<'KeyType, 'a>) =
        s
        |> Series.groupBy (fun k v -> dropAllPropertiesBut newKeys k )
        |> Series.mapValues (fun valueCol -> 
            let fInit = valueCol |> Series.values |> Array.ofSeq |> op
            let filterCol = 
                valueCol
                |> Series.mapValues fInit
            filterCol 
        )
        |> Series.values
        |> Series.mergeAll

    ///
    let createGroupFilter (op :'a [] -> 'a -> bool) (newKeys:seq<string>) (s: Series<'KeyType, 'a>) =
        groupTransform op newKeys s

    ///
    let aggregate (op : seq<'a> -> 'b) (newKeys:seq<string>) (filters:seq<Series<Key,bool>>) (col:Series<Key,'a>) = 
        let filtered = 
            filters
            |> Seq.map (fun s -> 
                System.Guid.NewGuid(),
                s |> Series.filterValues id)
            |> Series.ofObservations
            |> Frame.ofColumns
            |> Frame.dropSparseRows
        let colID = (System.Guid.NewGuid())
        filtered
        |> Frame.addCol colID col
        |> Frame.dropSparseRows
        |> Frame.applyLevel (dropAllPropertiesBut newKeys) (fun s -> s |> Series.values |> op)
        |> Frame.getCol colID

    ///
    let assemble (cols:seq<(string * #ISeries<Key>)>) =
        Frame.ofColumns cols

    ///
    let pivot (pivotCol:string) (assembeledFrame:Frame<Key,string>) =
        assembeledFrame
        |> Frame.nestBy (fun k -> 
            let value: string option = k.TryGetTypedValue pivotCol
            value.Value
            )
        |> Series.map (fun k f ->
                f
                |> Frame.mapColKeys (fun ck -> sprintf "%s.%s" k ck)
                |> Frame.mapRowKeys (dropProperties [pivotCol])
            )
        |> Series.values
        |> Frame.mergeAll

//module Library = 
//    open System.IO
  

//    let labeledQuantification (correlationFilter:float option) (useModifiedPeptides:bool) (treatModifiedPeptidesAs:bool) (outputDir:string) (instrumentOutput:string[]) =

//        let logger = Logging.createLogger "LabeledQuantification"

//        logger.Trace (sprintf "Input files: %A" instrumentOutput)
//        logger.Trace (sprintf "Output directory: %s" outputDir)

//        // initialize Reader and Transaction
//        let outFilePath =
//            let fileName = "LabeledQuant.txt"
//            Path.Combine [|outputDir;fileName|]
//        logger.Trace (sprintf "Result file path: %s" outFilePath)

//        let peptidesAndProteins = 
//            instrumentOutput
//            |> Array.map ( fun fp ->
//                Csv.CsvReader<ProteinAssignedQuantifiedIon>(SchemaMode=Csv.Fill).ReadFile(fp,'\t',false,1)
//                |> Array.ofSeq
//                |> Frame.ofRecords
//                )
//            |> Frame.mergeAll
        
//        let keyCols =         
//            [|
//                "FileName"      
//                "ProteinGroup"  
//                "StringSequence"
//                "PepSequenceID" 
//                "ModSequenceID" 
//                "Charge"        
//                "GlobalMod"     
//            |]
//        /// Step 1: Aggregate GlobalModifications
//        ///
//        let peptidesAndProteinsIndexed = Core.indexWithColumnValues keyCols peptidesAndProteins
//        /// Cols
//        let heavy       = peptidesAndProteinsIndexed |> Core.getColumn<float>"MeasuredApex_Heavy"
//        let light       = peptidesAndProteinsIndexed |> Core.getColumn<float>"MeasuredApex_Light"
//        let correlation = peptidesAndProteinsIndexed |> Core.getColumn<float>"Correlation_Light_Heavy"
//        let gMod        = peptidesAndProteinsIndexed |> Core.getColumn<bool>"GlobalMod"
//        /// Zipped
//        let N14ByN15 = Core.zip (fun x y -> x / y) light heavy
//        /// Filters
//        let correlationF = 
//            match correlationFilter with 
//            | Some c -> 
//                Core.createFilter (fun x -> x > c ) correlation
//            | None -> 
//                Core.createFilter (fun x -> x > 0. ) correlation
//        let filterOutLight = Core.createFilter id gMod
//        let filterOutHeavy = Core.createFilter (id >> not) gMod
//        /// Aggregated
//        let light_agg :Series<_,float> = aggregate Seq.mean aggregationLevel [] light 
//        let heavy_agg :Series<_,float> = aggregate Seq.mean aggregationLevel [] heavy 
//        let ratio_agg :Series<_,float> = aggregate Seq.mean aggregationLevel [] N14ByN15 
//        let correlation_agg :Series<_,float> = aggregate Seq.mean aggregationLevel [] correlation 
//        let patternLight_Mz_agg :Series<_,string> = aggregate (Seq.item 0)  aggregationLevel [filterOutHeavy] lightPatternMz
//        let patternLight_Int_agg :Series<_,string> = aggregate (Seq.item 0)  aggregationLevel [filterOutHeavy] lightPatternI
//        let patternHeavy_Mz_agg :Series<_,string> = aggregate (Seq.item 0)  aggregationLevel [filterOutLight] heavyPatternMz
//        let patternHeavy__Int_agg :Series<_,string> = aggregate (Seq.item 0)  aggregationLevel [filterOutLight] heavyPatternI
//        /// Assembeled
//        let assembeled = 
//            assemble 
//                [
//                "Light"  , light_agg :> ISeries<Key>
//                "Heavy"  , heavy_agg :> ISeries<Key>
//                "Ratio"  , ratio_agg :> ISeries<Key>
//                "CorrelationLightHeavy", correlation_agg :> ISeries<Key>
//                "lightPatternMz", patternLight_Mz_agg :> ISeries<Key>
//                "lightPatternI", patternLight_Int_agg :> ISeries<Key>
//                "heavyPatternMz", patternHeavy_Mz_agg :> ISeries<Key>
//                "heavyPatternI", patternHeavy__Int_agg :> ISeries<Key>
//                ]
//        /// Step 2: Aggregate Charges

//        /// Step 3: Aggregate Modifications

//        /// Step 3: Aggregate Peptides

//        let toSave = 
//            assembeled
//            |> pivot "FileName"
//            |> rowKeyToColumns 
//        toSave.Print()
//        toSave.SaveCsv((Path.Combine [|outDir; "Quantifications.txt"|]),separator='\t',includeRowKeys=false)