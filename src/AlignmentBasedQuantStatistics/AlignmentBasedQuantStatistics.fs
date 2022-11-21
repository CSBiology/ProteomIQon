namespace ProteomIQon

open Deedle
open Microsoft
open Microsoft.ML
open Microsoft.ML.Data
open ProteomIQon.DTW'
open ProteomIQon.FDRControl'
open FSharpAux
open Plotly.NET

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
    ///
    let downcastPipeline (x : IEstimator<_>) = 
        match x with 
        | :? IEstimator<ITransformer> as y -> y
        | _ -> failwith "downcastPipeline: expecting a IEstimator<ITransformer>"

    let stringToArray s =
        s
        |> String.split ';'
        |> Array.map float
    
    let toPeptideForLearning positiveSet (s: ObjectSeries<string>) : PeptideForLearning=
        {
            Source_ApexIntensity = s.GetAs<float32>("align_ApexIntensity_SourceFile")
            Source_Intensities   = s.GetAs<float32>("align_Quant_SourceFile")
            DtwDistanceBefore    =
                let source = Array.zip (s.GetAs<string>("align_RtTrace_SourceFile") |> stringToArray) (s.GetAs<string>("align_IntensityTrace_SourceFile") |> stringToArray |> zNorm)
                let target = 
                    Array.zip 
                        ( 
                            if s.GetAs<bool>("GlobalMod") then
                                s.GetAs<string>("RtTrace_Heavy") |> stringToArray
                            else
                                s.GetAs<string>("RtTrace_Light") |> stringToArray
                        )
                        ( 
                            if s.GetAs<bool>("GlobalMod") then
                                s.GetAs<string>("IntensityTrace_Corrected_Heavy") |> stringToArray |> zNorm
                            else
                                s.GetAs<string>("IntensityTrace_Corrected_Light") |> stringToArray |> zNorm
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
                    abs (s.GetAs<float>("align_ScanTime_SourceFile") - (stringToArray (s.GetAs<string>("Params_Heavy"))).[1])
                    |> float32
                else
                    abs (s.GetAs<float>("align_ScanTime_SourceFile") - (stringToArray (s.GetAs<string>("Params_Light"))).[1])
                    |> float32
            Label = positiveSet
        }

    let createTrainingsData fullQuant align alignQuant=
        let quant =
            Frame.ReadCsv(fullQuant, true, separators = "\t")
            |> Frame.indexRowsUsing (fun s ->
                s.GetAs<string>("StringSequence"),
                s.GetAs<bool>("GlobalMod"),
                s.GetAs<int>("Charge"),
                s.GetAs<int>("PepSequenceID"),
                s.GetAs<int>("ModSequenceID")
            )
            |> Frame.mapColKeys(fun ck -> "quant_" + ck)

        let align = 
            Frame.ReadCsv(align, true, separators = "\t")
            |> Frame.indexRowsUsing (fun s ->
                s.GetAs<string>("StringSequence"),
                s.GetAs<bool>("GlobalMod"),
                s.GetAs<int>("Charge"),
                s.GetAs<int>("PepSequenceID"),
                s.GetAs<int>("ModSequenceID")
            )
            |> Frame.mapColKeys(fun ck -> "align_" + ck)

        let alignedQuant = 
            Frame.ReadCsv(alignQuant, true, separators = "\t")
            |> Frame.filterRows (fun k s ->
                s.GetAs<string>("QuantificationSource") = "Alignment"
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

        let difference =
            align
            |> Frame.filterRows (fun k s ->
                quant.RowKeys
                |> Seq.contains k
                |> not
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

        let differenceSet =
            Frame.join JoinKind.Inner difference alignedQuant

        let finalSet =
            quant
            |> Frame.mapColKeys (fun ck -> ck.Replace("quant_",""))

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
                            abs (quantMap.QuantMzHeavy.[rk] - s.GetAs<float>("QuantMz_Heavy")) < (quantMap.QuantMzHeavy.[rk] * 0.01) &&
                            abs (quantMap.QuantHeavy.[rk] - s.GetAs<float>("Quant_Heavy")) < (quantMap.QuantHeavy.[rk] * 0.1) then
                            Some (toPeptideForLearning true s)
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
                            abs (quantMap.QuantMzLight.[rk] - s.GetAs<float>("QuantMz_Light")) < (quantMap.QuantMzLight.[rk] * 0.01) &&
                            abs (quantMap.QuantLight.[rk] - s.GetAs<float>("Quant_Light")) < (quantMap.QuantLight.[rk] * 0.1) then
                            Some (toPeptideForLearning true s)
                        else
                            Some (toPeptideForLearning false s)
                    | None -> None
            )
            |> Series.values
            |> Seq.choose id
            
        let pepForLearningToTakeMap =
            differenceSet
            |> Frame.mapRows (fun rk s ->
                toPeptideForLearning true s
            )
            |> Series.observations
            |> Map.ofSeq
    
        trainingSet, finalSet, pepForLearningToTakeMap
        
    let assignScoreAndQValue (matchedFiles: (string*string*string)[]) (logger: NLog.Logger) parallelismLevel diagnosticCharts outputDirectory =
        let trainingsData, allAlignments, pepForLearningToTakeMap =
            matchedFiles
            |> FSharpAux.PSeq.map (fun (quantFilePath,alignfilePath,alignQuantFilePath) -> createTrainingsData quantFilePath alignfilePath alignQuantFilePath)
            |> FSharpAux.PSeq.withDegreeOfParallelism parallelismLevel
            |> Seq.unzip3

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
            |> Chart.SaveHtmlAs (System.IO.Path.Combine(outputDirectory,"ProbabilityHistogram"))

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
            |> Chart.SaveHtmlAs (System.IO.Path.Combine(outputDirectory,"QValueDistribution"))

        allAlignments
        |> Seq.mapi (fun i file ->
            file
            |> Frame.filterRows (fun rk s ->
                pepForLearningToTakeMap
                |> Seq.item i
                |> Map.containsKey rk ||
                (s.GetAs<string>("QuantificationSource") = "PSM")
            )
            |> Frame.mapRows (fun rk s ->
                if s.GetAs<string>("QuantificationSource") = "Alignment" then
                    let prediction =
                        pepForLearningToTakeMap
                        |> Seq.item i
                        |> fun map -> map.[rk]
                        |> predict
                    let qValue =
                        prediction.Score
                        |> float
                        |> qValueStorey
                    s 
                    |> Series.replace "AlignmentScore" prediction.Score
                    |> Series.replace "AlignmentQValue" qValue
                else
                    s
            )
            |> Frame.ofRows
            |> fun frame ->
                frame.SaveCsv(System.IO.Path.Combine(outputDirectory, matchedFiles.[i] |> (fun (quantFilePath,alignfilePath,alignQuantFilePath) -> quantFilePath) |> System.IO.Path.GetFileName), includeRowKeys = false)
        )
    
            
            
            