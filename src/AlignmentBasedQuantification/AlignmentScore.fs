namespace ProteomIQon

open System.IO
open Argu
open System
open ProteomIQon.Core
open Core.MzIO
open Dto
open FSharp.Stats
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute
open Plotly.NET
open BioFSharp
open MzIO.Processing
open BioFSharp.Mz.SearchDB
open System.Data
open Microsoft.ML
open Microsoft.ML.Data

module AlignmentScore =

    [<CLIMutable>]
    type PeptideForLearning = 
        {
            [<ColumnName("X_ApexIntensity")>]
            Source_ApexIntensity                 : float32
            [<ColumnName("X_Intensities")>]
            Source_Intensities                   : float32
            [<ColumnName("X_Stabw")>]
            Source_Stabw                         : float32
            [<ColumnName("DtwDistanceBefore")>]
            DtwDistanceBefore               : float32
            [<ColumnName("Y_ApexIntensity")>]
            Target_ApexIntensity                 : float32
            [<ColumnName("Y_ReQuant")>]
            Target_ReQuant                       : float32
            [<ColumnName("Y_Intensity")>]
            Target_Intensity                     : float32
            [<ColumnName("DistancePredFinal")>]
            DistancePredFinal               : float32
            [<ColumnName("yDistance")>]
            yDistance                       : float32
            [<ColumnName("xDistance")>]
            xDistance                       : float32
            [<ColumnName("ydiff")>]
            ydiff                           : float32
            [<ColumnName("xdiff")>]
            xdiff                           : float32 
            [<ColumnName("Label")>]
            Label                           : bool
        }

    //[<CLIMutable>]
    //type ClassPredition = 
    //    {
    //        [<ColumnName("Score")>]
    //        Score       : float32
    //        [<ColumnName("Probability")>]
    //        Probability : float32
    //    }
    /////
    //let downcastPipeline (x : IEstimator<_>) = 
    //    match x with 
    //    | :? IEstimator<ITransformer> as y -> y
    //    | _ -> failwith "downcastPipeline: expecting a IEstimator<ITransformer>"


    ////Create the MLContext to share across components for deterministic results
    //let mlContext = MLContext(seed = Nullable 1) // Seed set to any number so you

    ////let fullData = mlContext.Data.LoadFromEnumerable comb
    //// STEP 1: Common data loading configuration   
    //let fullData = mlContext.Data.LoadFromEnumerable quantifiedV'

    //let ttdata = mlContext.Data.TrainTestSplit(data=fullData,testFraction = 0.8) 
      
    ////STEP 2: Process data, create and train the model 
    //let pipeline = 
    //    let trainer = mlContext.BinaryClassification.Trainers.FastTree(featureColumnName="Features",labelColumnName="Label")
    //    // Process data transformations in pipeline
    //    (mlContext.Transforms.Concatenate(
    //            "Features",
    //            "X_ApexIntensity",
    //            "X_Intensities",
    //            //"X_Test",
    //            "X_Stabw",
    //            "DtwDistanceBefore",
    //            "Y_ApexIntensity",
    //            "Y_ReQuant",
    //            "DistancePredFinal",
    //            "xDistance",
    //            "yDistance",
    //            "ydiff",
    //            "xdiff"
    //            )
    //        |> downcastPipeline
    //        ).Append(trainer)

    //let model = pipeline.Fit ttdata.TrainSet

    ////// STEP3: Run the prediciton on the test data
    //let predictions =
    //    model.Transform ttdata.TestSet

    //let metrics = 
    //    mlContext.BinaryClassification.Evaluate(predictions)
    //metrics
    ////PREDICTED || positive | negative | Recall
    ////TRUTH     ||======================
    //// positive ||    5,636 |      648 | 0.8969
    //// negative ||      776 |    5,613 | 0.8785
    ////          ||======================
    ////Precision ||   0.8790 |   0.8965 |
    ////"
    //metrics.ConfusionMatrix.GetFormattedConfusionTable()    


    //let predF = mlContext.Model.CreatePredictionEngine<PeptideForLearning,ClassPredition>(model)

    //let predict qp = 
    //    qp
    //    |> predF.Predict