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

    let aggregateWithTukey ((tukeyC,method): float*Domain.Transform) (agMethod: Domain.AggregationMethod) (logger: NLog.Logger)(series: Series<'K1,float>): float =
        let values = 
            series.Values
            |> Seq.map float
            |> Array.ofSeq
            |> removeNan
        //logger.Trace (sprintf "Unfiltered Values: %A" values)
        let transformedValues = 
            values
            |> Array.map (transform method)
        //logger.Trace (sprintf "Transformed Values: %A" transformedValues)
        let borders =
            transformedValues
            |> FSharp.Stats.Testing.Outliers.tukey tukeyC
        //logger.Trace (sprintf "Upper Border: %f" borders.Upper)
        //logger.Trace (sprintf "Lower Border: %f" borders.Lower)
        let filteredValues =
            transformedValues
            |> Array.filter (fun v -> v <= borders.Upper && v >= borders.Lower)
            |> Array.map (revertTransform method)
        //logger.Trace (sprintf "Filtered Values: %A" filteredValues)
        let res =
            filteredValues 
            |> (aggregationMethodArray agMethod)
        //logger.Trace (sprintf "Result: %f\n" res)
        res

    let seriesCV (series:Series<'R,float>) =
        series.Values
        |> Seq.toArray
        |> removeNan
        |> Seq.cv

    let seriesStDev (series:Series<'R,float>) =
        series.Values
        |> Seq.toArray
        |> removeNan
        |> Seq.stDev

    let seriesSEM (series:Series<'R,float>) =
        series.Values
        |> Seq.toArray
        |> removeNan
        |> fun x -> 
            let stDev = Seq.stDev x
            stDev / (sqrt (float x.Length))

    let matchStatMeasurement (param:Domain.StatisticalMeasurement) =
        match param with
        |Domain.StatisticalMeasurement.CV    -> seriesCV, "_CV"
        |Domain.StatisticalMeasurement.StDev -> seriesStDev, "_StDev"
        |Domain.StatisticalMeasurement.SEM   -> seriesSEM, "_SEM"

    // works like applyLevel, but columns can be excluded and have a different function applied
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

    // checks for duplicate column keys in the tables and renames the second if duplicate is present
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
                    // determines which filed should be included in the final table, that aren't filtered with Tukey
                    let fieldsWoTukey = 
                        // all columns that are in the final table
                        let allFields =
                            let duplicateProtectedProtCols = 
                                param.ProtColumnsOfInterest
                                |> Array.map (fun col ->
                                    if uniqueKeyTables.sameKeys.Contains col then
                                        col + "_Prot"
                                    else
                                        col
                                    )
                            Array.append param.QuantColumnsOfInterest duplicateProtectedProtCols
                            |> fun x ->
                                if labeled then 
                                    Array.append x [|"Ratio"|]
                                else
                                    x
                        // fields that are in the final table but filtered with tukey
                        let tukeyFields =
                            param.Tukey
                            |> Array.map (fun (name,_,_) -> name)
                        allFields
                        |> Array.filter (fun x -> Array.contains x tukeyFields |> not)
                    // calculates distinct peptide count based on the amount of peptides used for ratio calculation before tukey filtering
                    loggerFile.Trace "calculating distinct peptide count"
                    let distinctPeptideCount =
                        if param.ProtColumnsOfInterest |> Array.contains "DistinctPeptideCount" && labeled then
                            let ratios = alignedTables.GetColumn<float>("Ratio")
                            ratios
                            |> Series.applyLevel (fun (prot,pep,id) -> prot) (Series.countValues)
                            |> Some
                        elif param.ProtColumnsOfInterest |> Array.contains "DistinctPeptideCount" then
                            let light = alignedTables.GetColumn<float>(param.EssentialFields.Light)
                            light
                            |> Series.applyLevel (fun (prot,pep,id) -> prot) (Series.countValues)
                            |> Some
                        else
                            None
                    loggerFile.Trace "processing columns with tukey property"
                    let tukeyColumns =
                        param.Tukey
                        |> Array.map (fun (name,tukeyC,method) ->
                            name,
                            alignedTables
                            |> Frame.getCol name
                            |> Series.applyLevel (fun (prot,pep,id) -> prot) (aggregateWithTukey (tukeyC, method) param.AggregatorPepToProt loggerFile)
                        )
                        |> Frame.ofColumns   
                    loggerFile.Trace "processing columns with no tukey property"
                    let noTukeyColumns =
                        fieldsWoTukey
                        |> Array.filter (fun x -> x <> "DistinctPeptideCount")
                        |> Array.map (fun name ->
                            name,
                            alignedTables
                            |> Frame.getCol name
                            |> Series.applyLevel (fun (prot,pep,id) -> prot) (aggregationMethodSeries param.AggregatorPepToProt)
                        )
                        |> Frame.ofColumns
                    loggerFile.Trace (sprintf "calculating stats for %A" param.StatisticalMeasurements)
                    let statsColumns =
                        param.StatisticalMeasurements
                        |> Array.map (fun (name, method) ->
                            let seriesMethod, columnNameExtension = matchStatMeasurement method
                            name + columnNameExtension,
                            alignedTables
                            |> Frame.getCol name
                            |> Series.applyLevel (fun (prot,pep,id) -> prot) seriesMethod
                        )
                        |> Frame.ofColumns
                    loggerFile.Trace "combining all columns"
                    let combinedColumns =
                        let s = (tukeyColumns |> Frame.addCol "DistinctPeptideCount" (distinctPeptideCount.Value))
                        [
                        s
                        noTukeyColumns
                        statsColumns
                        ]
                        |> Frame.mergeAll
                    combinedColumns
                    |> Frame.sortColsByKey
                    |> Frame.mapRowKeys (fun prot -> prot, System.IO.Path.GetFileNameWithoutExtension protFile)
                alignedAggTables
            ) quantFiles protFiles
        logger.Trace "Merging all results"
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