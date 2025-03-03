namespace ProteomIQon

open Deedle
open Microsoft
open Microsoft.ML
open Microsoft.ML.Data
open ProteomIQon.DTW'
open ProteomIQon.FDRControl'
open FSharpAux
open Plotly.NET
open FSharpAux.IO
open FSharpAux.IO.SchemaReader

module AlignmentBasedQuantStatistics =
    
    
    [<CLIMutable>]
    type PeptideForLearning = 
        {
            [<ColumnName("X_ApexIntensity")>]
            Source_ApexIntensity                 : float32
            [<ColumnName("X_Intensities")>]
            Source_Intensities                   : float32
            [<ColumnName("DtwDistanceBefore")>]
            DtwDistanceBefore               : float32
            [<ColumnName("Y_ApexIntensity")>]
            Target_ApexIntensity                 : float32
            [<ColumnName("Y_Intensities")>]
            Target_Intensity                     : float32
            [<ColumnName("ScanTimeDifference")>]
            ScanTimeDifference              : float32
            [<ColumnName("Label")>]
            Label                           : bool
        }

    [<CLIMutable>]
    type ClassPredition = 
        {
            [<ColumnName("Score")>]
            Score       : float32
            [<ColumnName("Probability")>]
            Probability : float32
        }

    let unzip4 sequence =
        let (lst1, lst2, lst3, lst4) = 
            Seq.foldBack (fun (a,b,c,d) (acc1, acc2, acc3, acc4) -> 
                a::acc1, b::acc2, c::acc3, d::acc4) sequence ([],[],[],[])
        (Array.ofList lst1, Array.ofList lst2, Array.ofList lst3, Array.ofList lst4)

    ///
    let downcastPipeline (x : IEstimator<_>) = 
        match x with 
        | :? IEstimator<ITransformer> as y -> y
        | _ -> failwith "downcastPipeline: expecting a IEstimator<ITransformer>"
    
    let toPeptideForLearning positiveSet (s: ObjectSeries<string>) : PeptideForLearning=
        {
            Source_ApexIntensity = s.GetAs<float32>("align_ApexIntensity_SourceFile")
            Source_Intensities   = s.GetAs<float32>("align_Quant_SourceFile")
            DtwDistanceBefore    =
                let source = Array.zip (s.GetAs<float[]>("align_RtTrace_SourceFile")) (s.GetAs<float[]>("align_IntensityTrace_SourceFile") |> zNorm)
                let target = 
                    Array.zip 
                        ( 
                            if s.GetAs<bool>("GlobalMod") then
                                s.GetAs<float[]>("RtTrace_Heavy")
                            else
                                s.GetAs<float[]>("RtTrace_Light")
                        )
                        ( 
                            if s.GetAs<bool>("GlobalMod") then
                                s.GetAs<float[]>("IntensityTrace_Corrected_Heavy") |> zNorm
                            else
                                s.GetAs<float[]>("IntensityTrace_Corrected_Light") |> zNorm
                        )
                distance None None None None None None (Array.map snd target) (Array.map snd source)
                |> float32
            Target_ApexIntensity =
                if s.GetAs<bool>("GlobalMod") then
                    s.GetAs<float32>("MeasuredApex_Heavy")
                else
                    s.GetAs<float32>("MeasuredApex_Light")
            Target_Intensity     =
                if s.GetAs<bool>("GlobalMod") then
                    s.GetAs<float32>("Quant_Heavy")
                else
                    s.GetAs<float32>("Quant_Light")
            ScanTimeDifference   = 
                if s.GetAs<bool>("GlobalMod") then
                    abs (s.GetAs<float>("align_ScanTime_SourceFile") - ((s.GetAs<float[]>("Params_Heavy"))).[1])
                    |> float32
                else
                    abs (s.GetAs<float>("align_ScanTime_SourceFile") - ((s.GetAs<float[]>("Params_Light"))).[1])
                    |> float32
            Label = positiveSet
        }

    let createDataToScore fullQuant align alignQuant=
        let quant =
            Csv.CsvReader<Dto.QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(fullQuant,'\t',false,1)
            |> Array.ofSeq
            |> Frame.ofRecords
            |> Frame.indexRowsUsing (fun s ->
                s.GetAs<string>("StringSequence"),
                s.GetAs<bool>("GlobalMod"),
                s.GetAs<int>("Charge"),
                s.GetAs<int>("PepSequenceID"),
                s.GetAs<int>("ModSequenceID")
            )
            |> Frame.mapColKeys(fun ck -> "quant_" + ck)

        let align = 
            Csv.CsvReader<Dto.AlignmentResult>(SchemaMode=Csv.Fill).ReadFile(align,'\t',false,1)
            |> Array.ofSeq
            |> Frame.ofRecords
            |> Frame.indexRowsUsing (fun s ->
                s.GetAs<string>("StringSequence"),
                s.GetAs<bool>("GlobalMod"),
                s.GetAs<int>("Charge"),
                s.GetAs<int>("PepSequenceID"),
                s.GetAs<int>("ModSequenceID")
            )
            |> Frame.mapColKeys(fun ck -> "align_" + ck)

        let alignedQuant = 
            Csv.CsvReader<Dto.QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(alignQuant,'\t',false,1)
            |> Array.ofSeq
            |> Frame.ofRecords
            |> Frame.filterRows (fun k s ->
                s.GetAs<Dto.QuantificationSource>("QuantificationSource") = Dto.QuantificationSource.Alignment
            )
            |> Frame.indexRowsUsing (fun s ->
                s.GetAs<string>("StringSequence"),
                s.GetAs<bool>("GlobalMod"),
                s.GetAs<int>("Charge"),
                s.GetAs<int>("PepSequenceID"),
                s.GetAs<int>("ModSequenceID")
            )

        let difference =
            align
            |> Frame.filterRows (fun k s ->
                quant.RowKeys
                |> Seq.contains k
                |> not
            )

        let differenceSet =
            Frame.join JoinKind.Inner difference alignedQuant

        let pepForLearningToTakeMap =
            differenceSet
            |> Frame.mapRows (fun rk s ->
                toPeptideForLearning true s
            )
            |> Series.observations
            |> Map.ofSeq
    
        alignedQuant, quant |> Frame.mapColKeys (fun ck -> ck |> String.replace"quant_" ""), pepForLearningToTakeMap

    let createTrainingsData fullQuant align alignQuant (alignScoreParams: Domain.AlignmentBasedQuantStatisticsParams) =
        let quant =
            Csv.CsvReader<Dto.QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(fullQuant,'\t',false,1)
            |> Array.ofSeq
            |> Frame.ofRecords
            |> Frame.indexRowsUsing (fun s ->
                s.GetAs<string>("StringSequence"),
                s.GetAs<bool>("GlobalMod"),
                s.GetAs<int>("Charge"),
                s.GetAs<int>("PepSequenceID"),
                s.GetAs<int>("ModSequenceID")
            )
            |> Frame.mapColKeys(fun ck -> "quant_" + ck)

        let align = 
            Csv.CsvReader<Dto.AlignmentResult>(SchemaMode=Csv.Fill).ReadFile(align,'\t',false,1)
            |> Array.ofSeq
            |> Frame.ofRecords
            |> Frame.indexRowsUsing (fun s ->
                s.GetAs<string>("StringSequence"),
                s.GetAs<bool>("GlobalMod"),
                s.GetAs<int>("Charge"),
                s.GetAs<int>("PepSequenceID"),
                s.GetAs<int>("ModSequenceID")
            )
            |> Frame.mapColKeys(fun ck -> "align_" + ck)

        let alignedQuant = 
            Csv.CsvReader<Dto.QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(alignQuant,'\t',false,1)
            |> Array.ofSeq
            |> Frame.ofRecords
            |> Frame.filterRows (fun k s ->
                s.GetAs<Dto.QuantificationSource>("QuantificationSource") = Dto.QuantificationSource.Alignment
            )
            |> Frame.indexRowsUsing (fun s ->
                s.GetAs<string>("StringSequence"),
                s.GetAs<bool>("GlobalMod"),
                s.GetAs<int>("Charge"),
                s.GetAs<int>("PepSequenceID"),
                s.GetAs<int>("ModSequenceID")
            )

        let overlap =
            align
            |> Frame.filterRows (fun k s ->
                quant.RowKeys
                |> Seq.contains k
            )
            
        let quantMap =
            let quantMzHeavy: Map<string*bool*int*int*int,float> =
                quant
                |> Frame.getCol "quant_QuantMz_Heavy"
                |> Series.observations
                |> Map.ofSeq
            let quantHeavy: Map<string*bool*int*int*int,float> =
                quant
                |> Frame.getCol "quant_Quant_Heavy"
                |> Series.observations
                |> Map.ofSeq
            let quantMzLight: Map<string*bool*int*int*int,float> =
                quant
                |> Frame.getCol "quant_QuantMz_Light"
                |> Series.observations
                |> Map.ofSeq
            let quantLight: Map<string*bool*int*int*int,float> =
                quant
                |> Frame.getCol "quant_Quant_Light"
                |> Series.observations
                |> Map.ofSeq
            {|QuantMzHeavy = quantMzHeavy; QuantHeavy = quantHeavy; QuantMzLight = quantMzLight; QuantLight = quantLight|}

        let overlapSet =
            Frame.join JoinKind.Inner overlap alignedQuant

        let trainingSet = 
            overlapSet
            |> Frame.mapRows (fun rk s ->
                if (rk |> fun (_,b,_,_,_) -> b) then
                    let originalQuant = 
                        quantMap.QuantHeavy
                        |> Map.tryFind rk
                    match originalQuant with
                    | Some o ->
                        if 
                            abs (quantMap.QuantMzHeavy.[rk] - s.GetAs<float>("QuantMz_Heavy")) < (quantMap.QuantMzHeavy.[rk] * (alignScoreParams.PositiveQuantMzCutoff)) &&
                            abs (quantMap.QuantHeavy.[rk] - s.GetAs<float>("Quant_Heavy")) < (quantMap.QuantHeavy.[rk] * (alignScoreParams.PositiveQuantCutoff)) then
                            Some (toPeptideForLearning true s)
                        elif
                            abs (quantMap.QuantMzHeavy.[rk] - s.GetAs<float>("QuantMz_Heavy")) > (quantMap.QuantMzHeavy.[rk] * (1. - alignScoreParams.NegativeQuantMzCutoff)) ||
                            abs (quantMap.QuantHeavy.[rk] - s.GetAs<float>("Quant_Heavy")) > (quantMap.QuantHeavy.[rk] * (1. - alignScoreParams.NegativeQuantCutoff)) then
                            Some (toPeptideForLearning false s)
                        else
                            Some (toPeptideForLearning false s)
                    | None -> None
                else
                    let originalQuant = 
                        quantMap.QuantLight
                        |> Map.tryFind rk
                    match originalQuant with
                    | Some o ->
                        if
                            abs (quantMap.QuantMzLight.[rk] - s.GetAs<float>("QuantMz_Light")) < (quantMap.QuantMzLight.[rk] * (alignScoreParams.PositiveQuantMzCutoff)) &&
                            abs (quantMap.QuantLight.[rk] - s.GetAs<float>("Quant_Light")) < (quantMap.QuantLight.[rk] * (alignScoreParams.PositiveQuantCutoff)) then
                            Some (toPeptideForLearning true s)
                        elif
                            abs (quantMap.QuantMzLight.[rk] - s.GetAs<float>("QuantMz_Light")) > (quantMap.QuantMzLight.[rk] * (1. - alignScoreParams.NegativeQuantMzCutoff)) ||
                            abs (quantMap.QuantLight.[rk] - s.GetAs<float>("Quant_Light")) > (quantMap.QuantLight.[rk] * (1. - alignScoreParams.NegativeQuantCutoff)) then
                            Some (toPeptideForLearning false s)
                        else
                            Some (toPeptideForLearning false s)
                    | None -> None
            )
            |> Series.values
            |> Seq.choose id
    
        trainingSet
        
    let assignScoreAndQValue ((quantFile,alignFile,alignQuantFile): (string*string*string)) (matchedFilesLearning: (string*string*string)[]) diagnosticCharts (alignScoreParams: Domain.AlignmentBasedQuantStatisticsParams) outputDirectory =

        let logger = Logging.createLogger (System.IO.Path.GetFileNameWithoutExtension quantFile)

        let chartDirectory = System.IO.Path.Combine(outputDirectory, System.IO.Path.GetFileNameWithoutExtension quantFile)
        System.IO.Directory.CreateDirectory(chartDirectory) |> ignore

        let trainingsData =
            matchedFilesLearning
            |> Array.map (fun (quantFilePath,alignfilePath,alignQuantFilePath) -> createTrainingsData quantFilePath alignfilePath alignQuantFilePath alignScoreParams)

        let alignedQuants, quants, pepForLearningToTakeMap =
            createDataToScore quantFile alignFile alignQuantFile

        let trainingsData' =
            trainingsData
            |> Seq.concat
            |> Array.ofSeq

        //Create the MLContext to share across components for deterministic results
        let mlContext = MLContext(seed = 1) // Seed set to any number

        //let fullData = mlContext.Data.LoadFromEnumerable comb
        // STEP 1: Common data loading configuration   
        let fullData = mlContext.Data.LoadFromEnumerable trainingsData'

        let ttdata = mlContext.Data.TrainTestSplit(data=fullData,testFraction = 0.8) 
      
        //STEP 2: Process data, create and train the model 
        let pipeline = 
            let trainer = mlContext.BinaryClassification.Trainers.FastTree(featureColumnName="Features",labelColumnName="Label")
            // Process data transformations in pipeline
            (
                mlContext.Transforms.Concatenate(
                    "Features",
                    "X_ApexIntensity",
                    "X_Intensities",
                    "DtwDistanceBefore",
                    "Y_ApexIntensity",
                    "Y_Intensities",
                    "ScanTimeDifference"
                )
                |> downcastPipeline
            ).Append(trainer)

        let model = pipeline.Fit ttdata.TrainSet

        //// STEP3: Run the prediciton on the test data
        let predictions =
            model.Transform ttdata.TestSet

        let metrics = 
            mlContext.BinaryClassification.Evaluate(predictions)

        logger.Trace $"{metrics.ConfusionMatrix.GetFormattedConfusionTable()}"    

        let predF = mlContext.Model.CreatePredictionEngine<PeptideForLearning,ClassPredition>(model)

        let predict qp = 
            qp
            |> predF.Predict

        let qValueStorey =
            calculateQValueStorey trainingsData' (fun x -> x.Label |> not) (fun x -> float (predict x).Score) (fun x -> float (predict x).Score)
        
        if diagnosticCharts then
            let setPositive,setNegative =
                trainingsData'
                |> Array.partition (fun x -> x.Label)
            let positive =
                setPositive
                |> Array.map predict
                |> Array.map (fun x -> x.Probability)
                |> Chart.Histogram
                |> Chart.withTraceName "Positive"
    
            let negative =
                setNegative
                |> Array.map predict
                |> Array.map (fun x -> x.Probability)
                |> Chart.Histogram
                |> Chart.withTraceName "Negative"

            [
                positive
                negative
            ]
            |> Chart.Combine
            |> Chart.withX_AxisStyle("Probability")
            |> Chart.withY_AxisStyle("Count")
            |> Chart.SaveHtmlAs (System.IO.Path.Combine(chartDirectory,"ProbabilityHistogram"))

            let positiveQVal =
                setPositive
                |> Array.map predict
                |> Array.map (fun x -> qValueStorey(float x.Score))
                |> Chart.Histogram
                |> Chart.withTraceName "Positive"
                
            let negativeQVal =
                setNegative
                |> Array.map predict
                |> Array.map (fun x -> qValueStorey(float x.Score))
                |> Chart.Histogram
                |> Chart.withTraceName "Negative"
                
            [
                positiveQVal
                negativeQVal
            ]
            |> Chart.Combine
            |> Chart.withX_AxisStyle("Q-Value")
            |> Chart.withY_AxisStyle("Count")
            |> Chart.SaveHtmlAs (System.IO.Path.Combine(chartDirectory,"QValueDistribution"))

        alignedQuants
        |> Frame.filterRows (fun rk s ->
            pepForLearningToTakeMap
            |> Map.containsKey rk
        )
        |> fun frame ->
            let scoreSeries =
                frame
                |> Frame.mapRows (fun rk s ->
                    let prediction =
                        pepForLearningToTakeMap
                        |> fun map -> map.[rk]
                        |> predict
                    float prediction.Score
                )
            let qValSeries =
                scoreSeries
                |> Series.map (fun rk v ->
                    let qValue =
                        v
                        |> float
                        |> qValueStorey
                    qValue
                )
            frame
            |> Frame.dropCol "AlignmentScore"
            |> Frame.dropCol "AlignmentQValue"
            |> Frame.addCol "AlignmentScore" scoreSeries
            |> Frame.addCol "AlignmentQValue" qValSeries
        |> Frame.merge quants
        |> Frame.mapValues (fun (v: obj) -> 
            if v :? System.Double [] then
                v :?> System.Double []
                |> Array.map string
                |> String.concat ";"
            else
                string v
        )
        |> fun frame ->
            logger.Trace $"Original Count: {quants.RowCount}"
            logger.Trace $"New Count: {frame.RowCount}"
            frame.SaveCsv(System.IO.Path.Combine(outputDirectory, quantFile |> System.IO.Path.GetFileName), includeRowKeys = false, separator = '\t')
    
            
            
            