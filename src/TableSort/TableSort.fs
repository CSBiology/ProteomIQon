namespace ProteomIQon

open Deedle
open FSharpAux.IO
open FSharpAux
open FSharp.Stats
open BioFSharp


module TableSort =

    let createSchema (fieldType: string) (fields: string[]) =
        fields
        |> Array.map (fun s -> sprintf "%s=%s"s fieldType)
        |> String.concat ","

    let aggregationMethodSeries (agMethod: Domain.AggregationMethod): Series<'K0,'V1> -> float =
        match agMethod with
        | Domain.AggregationMethod.Sum    -> Stats.sum
        | Domain.AggregationMethod.Mean   -> Stats.mean
        | Domain.AggregationMethod.Median -> Stats.median

    let aggregationMethodArray (agMethod: Domain.AggregationMethod): float[] -> float =
        match agMethod with
        | Domain.AggregationMethod.Sum ->
            fun (x: float[]) ->
                if x.Length = 0 then nan
                else Array.sum x
        | Domain.AggregationMethod.Mean ->
            fun (x: float[]) ->
                if x.Length = 0 then nan
                else Array.average x
        | Domain.AggregationMethod.Median ->
            fun (x: float[]) ->
                if x.Length = 0 then nan
                else Array.median x

    let transform (method:Domain.Transform) =
        match method with
        |Domain.Transform.Log10       -> log10
        |Domain.Transform.Log2        -> log2
        |Domain.Transform.Ln          -> fun x -> System.Math.Log (x,System.Math.E)
        |Domain.Transform.NoTransform -> id

    let revertTransform (method:Domain.Transform) =
        match method with
        |Domain.Transform.Log10       -> fun x -> 10.**x
        |Domain.Transform.Log2        -> revLog2
        |Domain.Transform.Ln          -> fun x -> System.Math.E**x
        |Domain.Transform.NoTransform -> id

    let removeNan = Array.filter (System.Double.IsNaN >> not)

    let aggregateWithTukey ((tukeyC,method): float*Domain.Transform) (agMethod: Domain.AggregationMethod) (logger: NLog.Logger)(series: Series<'K1,'V2>): float =
        let values = 
            series.Values
            |> Seq.map float
            |> Array.ofSeq
            |> removeNan
        logger.Trace (sprintf "Unfiltered Values: %A" values)
        let transformedValues = 
            values
            |> Array.map (transform method)
        logger.Trace (sprintf "Transformed Values: %A" transformedValues)
        let borders =
            transformedValues
            |> FSharp.Stats.Testing.Outliers.tukey tukeyC
        logger.Trace (sprintf "Upper Border: %f" borders.Upper)
        logger.Trace (sprintf "Lower Border: %f" borders.Lower)
        let filteredValues =
            transformedValues
            |> Array.filter (fun v -> v <= (borders.Upper + 0.000001) && v >= (borders.Lower - 0.000001))
            |> Array.map (revertTransform method)
        logger.Trace (sprintf "Filtered Values: %A" filteredValues)
        let res =
            filteredValues 
            |> (aggregationMethodArray agMethod)
        logger.Trace (sprintf "Result: %f\n" res)
        res

    let seriesCV (series:Series<'R,float>) =
        series.Values
        |> Seq.toArray
        |> removeNan
        |> Seq.cv

    let getAggregatedPeptidesVals (peptidesPresent: Set<string>) (peptidesMapped:Series<string,string[]>) (data: Frame<string,string>) (columnName:string) (agMethod: Domain.AggregationMethod) (tukey: (string*float*Domain.Transform)[]): Series<string,float> =
        peptidesMapped
        |> Series.mapValues (fun peptides ->
            peptides
            |> Array.map (fun peptide ->
                if peptidesPresent.Contains(peptide) then
                    Some (try Some (data.Item(columnName).Item(peptide)) with _ -> None)
                else
                    None
            )
        )
        |> Series.mapValues (fun x ->
            x
            |> Array.choose id
            |> Array.choose id
            |> fun arr ->
                let tukeyField =
                    tukey
                    |> Array.tryFind (fun (fieldName,_,_) -> fieldName = columnName)
                match tukeyField with
                |Some (name,tukeyC,method)->
                    let arrTr = arr |> Array.map (transform method)
                    let borders = FSharp.Stats.Testing.Outliers.tukey tukeyC arrTr
                    arrTr
                    |> Array.filter (fun v -> v <  borders.Upper && v > borders.Lower)
                    |> Array.map (revertTransform method)
                    |> (aggregationMethodArray agMethod)
                |None ->
                    arr
                    |> (aggregationMethodArray agMethod)
        )

    let getAggregatedCVVals (peptidesPresent: Set<string>) (peptidesMapped:Series<string,string[]>) (data: Frame<string,string>) (columnName:string) (tukey: (string*float*Domain.Transform)[]) =
        let aggregatedValues =
            peptidesMapped
            |> Series.mapValues (fun peptides ->
                peptides
                |> Array.map (fun peptide ->
                    if peptidesPresent.Contains(peptide) then
                        Some (try Some (data.Item(columnName).Item(peptide)) with _ -> None)
                    else
                        None
                )
            )
        let tukeyField =
            tukey
            |> Array.tryFind (fun (fieldName, _,_) -> fieldName = columnName)
        let cvNoTukey =
            aggregatedValues
            |> Series.mapValues (fun x ->
                x
                |> Array.choose id
                |> Array.choose id
                |> Seq.cv
            )
        let cvTukey: Series<string,float> option=
            match tukeyField with
            |Some (name,tukeyC,method) ->
                let res = 
                    aggregatedValues
                    |> Series.mapValues (fun x ->
                        x
                        |> Array.choose id
                        |> Array.choose id
                        |> fun arr ->
                            let arrTr = arr |> Array.map (transform method)
                            let borders = FSharp.Stats.Testing.Outliers.tukey tukeyC arrTr
                            arrTr
                            |> Array.filter (fun v -> v <  borders.Upper && v > borders.Lower)
                            |> Array.map (revertTransform method)
                            |> Seq.cv
                    )
                Some res
            |None -> None
        {|CV = cvNoTukey; CVTukey = cvTukey|}

    let applyLevelWithException (levelSel:_ -> 'K) (ex: 'C[]) (op:_ -> 'T) (exOp:_ -> 'T) (frame:Frame<'R, 'C>) =
        let indexBuilder = Deedle.Indices.Linear.LinearIndexBuilder.Instance
        let vectorBuilder = Deedle.Vectors.ArrayVector.ArrayVectorBuilder.Instance
        frame.GetColumns<'T>()
        |> Series.map (fun c s ->
                if (Array.contains c ex) then
                    Series.applyLevel levelSel exOp s
                else
                    Series.applyLevel levelSel op s)
        |> FrameUtils.fromColumns indexBuilder vectorBuilder

    let ensureUniqueKeys (frame1: Frame<'a,string>) (frame2: Frame<'b,string>) (differentiator: string)=
        let colKeys1 = frame1.ColumnKeys |> Set.ofSeq
        let colKeys2 = frame2.ColumnKeys |> Set.ofSeq
        let sameKeys = Set.intersect colKeys1 colKeys2
        let newFrame2 = 
            frame2 
            |> Frame.mapColKeys (fun ck -> 
                if sameKeys.Contains ck then
                    ck+differentiator
                else ck
            )
        {|Frame1 = frame1; Frame2 = newFrame2; sameKeys = sameKeys|}

    let peptideEvidenceClassToFloat (ec: string) =
        match ec with
        | "Unknown" -> 0.
        | "C1a" -> 1.
        | "C1b" -> 2.
        | "C2a" -> 3.
        | "C2b" -> 4.
        | "C3a" -> 5.
        | "C3b" -> 6.
        | _ -> failwith "invalid peptide evidence class"

    let filterFrame (filterOn: Domain.FilterOnField[]) (frame: Frame<'a,string>) =
        let rec filterLoop (i: int) (frameAcc: Frame<'a,string>) =
            if i = filterOn.Length then
                frameAcc
            else
                let currentField = filterOn.[i]
                let frame' =
                    frameAcc
                    |> Frame.filterRowValues (fun s ->
                        if currentField.FieldName = "Class" then
                            match s.TryGetAs<string>(currentField.FieldName) with
                            | OptionalValue.Missing -> failwith "EvidenceClass can not be missing"
                            | OptionalValue.Present v ->
                                match (currentField.UpperBound, currentField.LowerBound) with
                                | None, None -> true
                                | Some uB, None -> (peptideEvidenceClassToFloat v) < uB
                                | None, Some lB -> (peptideEvidenceClassToFloat v) > lB
                                | Some uB, Some lB ->
                                    (peptideEvidenceClassToFloat v) < uB && (peptideEvidenceClassToFloat v) > lB
                        else
                            match s.TryGetAs<float>(currentField.FieldName) with
                            | OptionalValue.Missing -> true
                            | OptionalValue.Present v ->
                                match (currentField.UpperBound, currentField.LowerBound) with
                                | None, None -> true
                                | Some uB, None -> v < uB
                                | None, Some lB -> v > lB
                                | Some uB, Some lB -> v < uB && v > lB
                    )
                filterLoop (i+1) frame'
        filterLoop 0 frame

    let sortTables (quantFiles: string[]) (protFiles: string[]) outDirectory (param: Domain.TableSortParams) =

        let logger = Logging.createLogger "TableSort"
        // creates a schema that sets the column type of every column containing values to float
        let quantColsWithValues =
            Array.append param.QuantColumnsOfInterest (param.QuantFieldsToFilterOn |> Array.map (fun x -> x.FieldName))
            |> Array.distinct
        logger.Trace (sprintf "Prot columns of interest: %A" param.ProtColumnsOfInterest)
        logger.Trace (sprintf "Quant cloumns of interest: %A" param.QuantColumnsOfInterest)
        logger.Trace (sprintf "Quant cloumns with values that are kept for analysis: %A" quantColsWithValues)
        let quantColumnTypes =
            createSchema "float" quantColsWithValues
        let labeled = param.EssentialFields.Heavy.IsSome
        logger.Trace (sprintf "Labeled experiment: %b" labeled)
        let quantColumnsOfInterest =
            if labeled then 
                param.QuantColumnsOfInterest
                // appends "Ratio" column since it's not included in the original columns with values
                |> Array.append [|"Ratio"|]
            else
                param.QuantColumnsOfInterest

        let isQuantColumnOfInterest = 
            let s = quantColumnsOfInterest |> Set.ofArray
            fun x -> s.Contains x

        let tables =
            Array.map2 (fun (quantFile:string) protFile ->
                let loggerFile = Logging.createLogger (System.IO.Path.GetFileNameWithoutExtension quantFile)
                // reads quant and prot table
                loggerFile.Trace "Reading quant table"
                let quantTable: Frame<int,string> =
                    Frame.ReadCsv(path=quantFile ,separators=param.SeparatorIn, schema=quantColumnTypes)
                loggerFile.Trace "Reading prot table"
                let protTable : Frame<string*string,string> =
                    Frame.ReadCsv(path=protFile ,separators=param.SeparatorIn)
                    |> Frame.groupRowsBy param.EssentialFields.ProteinIDs
                    |> Frame.groupRowsBy param.EssentialFields.PepSequences
                    |> Frame.mapRowKeys (fun (ss,(id,_)) -> id,ss)
                loggerFile.Trace "Filtering quant table"
                loggerFile.Trace (sprintf "Calculating ratios: %b" labeled)
                let quantTableFiltered: Frame<int,string> =
                    // calculate 14N/15N ratios based on the corresponding fields from the quant table
                    // Ratio gets added before filtering, so that it can also be filtered on
                    if labeled then
                        let ratios =
                            let n14 = quantTable.GetColumn<float>param.EssentialFields.Light
                            let n15 = quantTable.GetColumn<float>param.EssentialFields.Heavy.Value
                            n14/n15
                        quantTable.AddColumn ("Ratio", ratios)

                    quantTable
                    |> filterFrame param.QuantFieldsToFilterOn
                loggerFile.Trace "Filtering prot table"
                let protTableFiltered: Frame<string*string,string> =
                    protTable
                    |> filterFrame param.ProtFieldsToFilterOn
                    |> Frame.expandRowsByKey (fun (prot,s) -> String.split ';' s |> Seq.map (fun pep -> prot,pep))
                loggerFile.Trace "Aggregating peptide versions"
                let quantTableAggregated: Frame<string,string> =
                    quantTableFiltered
                    // group all rows based on the peptide sequence. This groups different charges and global modifications of the same peptide
                    //|> Frame.groupRowsUsing (fun k ser -> (ser.GetAs<string>("StringSequence")))
                    |> Frame.groupRowsBy param.EssentialFields.PepSequence
                    // makes sure only columns with values are present (i.e. charge/global mod are removed)
                    |> Frame.filterCols (fun ck _ -> isQuantColumnOfInterest ck)
                    // aggregates the columns over the peptide sequence with a defined method (i.e. average,median,...)
                    |> applyLevelWithException (fun (sequence,index) -> sequence) ([|Some param.EssentialFields.Light;param.EssentialFields.Heavy|] |> Array.choose id)
                        (aggregationMethodSeries param.AggregatorFunction) (aggregationMethodSeries param.AggregatorFunctionIntensity)
                loggerFile.Trace "Ensuring unique Keys"
                let uniqueKeyTables = ensureUniqueKeys quantTableAggregated protTableFiltered "_Prot"
                loggerFile.Trace (sprintf "Keys in common that were changed: %A" uniqueKeyTables.sameKeys)
                let alignedTables =
                    Frame.align snd id uniqueKeyTables.Frame2 uniqueKeyTables.Frame1
                    |> Frame.mapRowKeys (fun (pep,(prot,_),(id)) -> prot,pep,id)
                    //|> Frame.applyLevel (fun (prot,pep,id) -> prot) (fun (s:Series<_,float>) -> s.Values |> Seq.mean)
                loggerFile.Trace "Aggregating peptides to proteins"
                let alignedAggTables =
                    let fieldsWoTukey = 
                        let allFields =
                            Array.append param.QuantColumnsOfInterest param.ProtColumnsOfInterest
                            |> fun x ->
                                if labeled then 
                                    Array.append x [|"Ratio"|]
                                else
                                    x
                        let tukeyFields =
                            param.Tukey
                            |> Array.map (fun (name,_,_) -> name)
                        allFields
                        |> Array.filter (fun x -> Array.contains x tukeyFields |> not)
                    let distinctPeptideCount =
                        if param.ProtColumnsOfInterest |> Array.contains "DistinctPeptideCount" && labeled then
                            [|"Ratio"|]
                            |> Array.map (fun name ->
                                alignedTables
                                |> Frame.sliceCols [name]
                                |> Frame.applyLevel (fun (prot,pep,id) -> prot) (fun (s: Series<string*string*string,float>) -> 
                                    let length = s.Values |> Seq.length
                                    float length
                                    )
                                |> Frame.mapColKeys (fun c -> "DistinctPeptideCount")
                            )
                        else
                            [||]
                    let tukeyColumns =
                        param.Tukey
                        |> Array.map (fun (name,tukeyC,method) ->
                            alignedTables
                            |> Frame.sliceCols [name]
                            |> Frame.applyLevel (fun (prot,pep,id) -> prot) (aggregateWithTukey (tukeyC, method) param.AggregatorPepToProt loggerFile)
                        )
                    let noTukeyColumns =
                        fieldsWoTukey
                        |> Array.map (fun name ->
                            alignedTables
                            |> Frame.sliceCols [name]
                            |> Frame.applyLevel (fun (prot,pep,id) -> prot) (aggregationMethodSeries param.AggregatorPepToProt)
                        )
                    let cvColumns =
                        param.CoefficientOfVariation
                        |> Array.map (fun name ->
                            alignedTables
                            |> Frame.sliceCols [name]
                            |> Frame.applyLevel (fun (prot,pep,id) -> prot) seriesCV
                            |> Frame.mapColKeys (fun c -> c + "_CV")
                        )
                    let combinedColumns =
                        [tukeyColumns;distinctPeptideCount;noTukeyColumns;cvColumns]
                        |> Array.concat
                        |> Frame.mergeAll
                    combinedColumns
                    |> Frame.sortColsByKey
                    |> Frame.mapRowKeys (fun prot -> prot, System.IO.Path.GetFileNameWithoutExtension protFile)
                alignedAggTables
            ) quantFiles protFiles
        tables
        |> Frame.mergeAll
        |> fun frame ->
            frame.SaveCsv (path=(outDirectory+(@"\TableSort.tab")), separator=(param.SeparatorOut))
            frame
            |> Frame.pivotTable
                (fun (proteinGroup, experiment) os -> proteinGroup)
                (fun (proteinGroup, experiment) os -> experiment)
                (fun f ->
                    f .GetRowAt<float> 0
                )
                |> Frame.expandAllCols 1
            |> fun framePiv -> framePiv.SaveCsv (path=(outDirectory+(@"\TableSort_horizontal.tab")), separator=param.SeparatorOut)