// Include CsbScaffold
#r "netstandard"
#load "../../.env/CsbScaffold.fsx"
#load "../../.env/packages/FSharpML/lib/netstandard2.0/FSharpML.fsx"
// If you want to use the wrappers for unmanaged LAPACK functions from of FSharp.Stats 
// include the path to the .lib folder manually to your PATH environment variable and make sure you set FSI to 64 bit
// use the following lines of code to ensure that LAPACK functionalities are enabled if you want to use them
// fails with "MKL service either not available, or not started" if lib folder is not included in PATH.
//open FSharp.Stats
//FSharp.Stats.Algebra.LinearAlgebra.Service()
 
open BioFSharp
open FSharpML

open Deedle
open FSharpAux
open FSharp.Stats
open System.IO
open Microsoft.ML
open Microsoft.ML.Data
open BioFSharp
open Microsoft.ML
open Microsoft.ML.Data   
open Microsoft.Data
open Microsoft.Data.DataView
open FSharp.Plotly

///
type Result = 
    {
    StringSequence               : string;
    GlobalMod                    : int;
    Charge                       : int;
    PrecursorMZ                  : float;
    MeasuredMass                 : float; 
    TheoMass                     : float;
    AbsDeltaMass                 : float;
    MeanPercolatorScore          : float;
    QValue                       : float;
    PEPValue                     : float;
    ProteinNames                 : string;
    N14QuantMz                   : float
    N14Quant                     : float
    N14Seo                       : float
    N14Params                    : string
    N15QuantMz                   : float
    N15Quant                     : float
    N15Seo                       : float
    N15Params                    : string
    N15Minus1QuantMz             : float
    N15Minus1Quant               : float
    N15Minus1Seo                 : float
    N15Minus1Params              : string
    }  
    
type AlignmentItem<'a> = {
    Sequence             : string
    GlobalMod            : bool
    Charge               : int
    Mz                   : float          
    MeanScanTime         : 'a
    Stabw                : 'a
    Intensity            : 'a
    CorrespondingScanTime: float
    N14Quant             : 'a
    N15Quant             : 'a
    Ratio                : float
    
    }

type Peptide = {
    Sequence             : string
    GlobalMod            : bool
    Charge               : int
    }

let createPeptideof (alignmentItem:AlignmentItem<'a>) = 
    {
    Sequence   = alignmentItem.Sequence
    GlobalMod  = alignmentItem.GlobalMod
    Charge     = alignmentItem.Charge
    }

let rawDataPath1 = Path.Combine [|__SOURCE_DIRECTORY__; "../" ;"data/181107_CK2_15m_GD2_01_9222.quant";|]
let rawDataPath2 = Path.Combine [|__SOURCE_DIRECTORY__; "../" ;"data/181107_CK2_15m_GD2_01_9221.quant";|]
let rawDataPath3 = Path.Combine [|__SOURCE_DIRECTORY__; "../" ;"data/181107_CK1_15m_GD1_01_9231.quant";|]
let rawDataPath4 = Path.Combine [|__SOURCE_DIRECTORY__; "../" ;"data/181107_CK1_15m_GD1_01_9230.quant";|]

let test = Frame.ReadCsv(rawDataPath1,separators="\t") 

let fileNameOf path = FSharpAux.IO.PathFileName.fileNameOfPath path
let files = 
    [
    rawDataPath1
    rawDataPath2
    rawDataPath3
    rawDataPath4        
    ]

let parseString s = 
    if s = "nan" then nan else float s

let getAlignmentItems (dataPath:string) = 
    Frame.ReadCsv(dataPath,separators="\t") 
    |> Frame.mapRowValues (fun x -> 
                            let sequence     = x.GetAs<string>("StringSequence")
                            let globalMod    = x.GetAs<bool>("GlobalMod")
                            let charge       = x.GetAs<int>("Charge")
                            let precMZ       = x.GetAs<float>("PrecursorMZ")
                            let N14Quant     = x.GetAs<string>("N14Quant")
                            let N15Quant     = x.GetAs<string>("N15Quant")
                           
                            let params       = 
                                if globalMod then 
                                    try
                                    x.GetAs<string>("N15Params") |> String.split ';'(* |> Array.map float*)
                                    |> Some
                                    with                                 
                                    | ex -> 
                                        None
                                else 
                                    try
                                    x.GetAs<string>("N14Params") |> String.split ';'(* |> Array.map float*)
                                    |> Some
                                    with
                                    | ex -> 
                                        None
                            match params with
                            | None -> None
                            | Some params ->
                                let meanScanTime = params.[1]
                                let stabw        = params.[2]
                                let intensity    = params.[0]

                                {
                                Sequence     = sequence
                                GlobalMod    = globalMod
                                Charge       = charge
                                Mz           = precMZ
                                MeanScanTime = meanScanTime
                                Stabw        = stabw
                                Intensity    = intensity
                                CorrespondingScanTime = nan
                                N14Quant = N14Quant 
                                N15Quant = N15Quant
                                Ratio = nan
                                }
                                |> Some 
                          ) 
    |> Series.values
    |> Seq.choose id
    |> Seq.choose (fun x -> 

                        try
                        {
                        Sequence     = x.Sequence     
                        GlobalMod    = x.GlobalMod    
                        Charge       = x.Charge       
                        Mz           = x.Mz           
                        MeanScanTime = float x.MeanScanTime 
                        Stabw        = float x.Stabw        
                        Intensity    = float x.Intensity    
                        CorrespondingScanTime = nan
                        N14Quant = parseString x.N14Quant 
                        N15Quant = parseString x.N15Quant
                        Ratio = parseString x.N14Quant / parseString x.N15Quant
                        } |> Some
                        with 
                        | _-> 
                            printfn "%s x.N14Quant, %s x.N15Quant" x.N14Quant x.N15Quant
                            None
                 )
    |> Array.ofSeq

let allAlignmentItems = 
    files
    |> Seq.map (fun x -> fileNameOf x, getAlignmentItems x)



let allIdentifiedPeptides = 
    allAlignmentItems
    |> Seq.map snd
    |> Seq.concat
    |> Seq.map createPeptideof
    |> Seq.distinct

let allSharedPeptides =
    allAlignmentItems
    |> Seq.map snd
    |> Seq.map (Seq.map createPeptideof)
    |> Seq.map (Set.ofSeq)
    |> Set.intersectMany

let filteredForSharedAligmentItems = 
    allAlignmentItems
    |> Seq.map snd
    |> Seq.map (Seq.filter (fun item -> allSharedPeptides.Contains (createPeptideof item)) >> Array.ofSeq)
    |> Array.ofSeq


let getAlignmentRanking = 
    [
    for i = 0 to filteredForSharedAligmentItems.Length-1 do 
        
        for j = 0 to filteredForSharedAligmentItems.Length-1 do 
            printfn "%i %i" i j 
            if j = i then 
                () 
            else
                let x = filteredForSharedAligmentItems.[i]
                printfn "peptides = %i" x.Length
                let y = filteredForSharedAligmentItems.[j]
                let (x',y') = 
                    x
                    |> Array.map (fun x -> x , Array.find (fun y -> createPeptideof x = createPeptideof y ) y)
                    |> Array.filter (fun x -> nan.Equals((fst x).MeanScanTime) |> not && nan.Equals((snd x).MeanScanTime) |> not)
                    |> Array.filter (fun x -> (fst x).Ratio < 5. && (fst x).Ratio > 0.01  && (snd x).Ratio < 5. && (snd x).Ratio > 0.01 )
                    |> Array.unzip
                printfn "peptides = %i" x'.Length
                printfn "calcPearson"                
                yield (i,j, Array.map2 (fun x y -> abs (x-y) ) (x' |> Array.map (fun x -> x.MeanScanTime)) (y'|> Array.map (fun x -> x.MeanScanTime)) |> Array.median) 
    ]
    
getAlignmentRanking    
|> List.groupBy (fun (t,s,dis) -> t)    

FSharp.Stats.Correlation.Seq.pearson [1.;2.;3.] [4.;5.;6.]

5+5
// let itemsY = getAlignmentItems rawData2

// let mapY = 
//     itemsY 
//     |> Array.map (fun x -> (x.Sequence,x.GlobalMod,x.Charge),x)
//     |> Map.ofArray

// let itemsForAlignment = 
//     itemsX
//     |> Array.choose (fun x -> 
//                         match mapY |> Map.tryFind ((x.Sequence,x.GlobalMod,x.Charge)) with
//                         | Some y -> {x with CorrespondingScanTime = y.MeanScanTime}  |> Some
//                         | None -> None
//                     )

// itemsForAlignment
// |> Array.map (fun x -> x.Stabw)
// |> Chart.BoxPlot
// |> Chart.Show

// itemsForAlignment
// |> Array.map (fun x -> x.MeanScanTime, x.CorrespondingScanTime)
// |> Chart.Point
// |> Chart.Show

// itemsForAlignment.Length
// /// Start MachineLearning
// let ctx = new MLContext()

// let trainTestSplit = 
//     ctx.Data.ReadFromEnumerable(itemsForAlignment)
//     |> FSharpML.Data.Regression.initTrainTestSplit(ctx,Testfraction=0.2)
// trainTestSplit.TestData.Preview().ColumnView
// 5+5    
// let model = 
//     EstimatorModel.create ctx
//     // Process data transformations in pipeline
//     |> EstimatorModel.appendBy (fun mlc -> mlc.Transforms.Conversion.ConvertType(DefaultColumnNames.Label,"CorrespondingScanTime",DataKind.R4))
//     |> EstimatorModel.appendBy (fun mlc -> 
//                                     mlc.Transforms.Concatenate(
//                                                                 "FeaturesRaw" , 
//                                                                 "Mz", 
//                                                                 "MeanScanTime", 
//                                                                 "Stabw", 
//                                                                 "Intensity"
//                                                                ) )
//     |> EstimatorModel.appendBy (fun mlc -> mlc.Transforms.Conversion.ConvertType(DefaultColumnNames.Features,"FeaturesRaw",DataKind.R4))
//     //|>  EstimatorModel.appendBy (fun mlc -> mlc.Transforms.Normalize(DefaultColumnNames.Features,DefaultColumnNames.Features,Microsoft.ML.Transforms.Normalizers.NormalizingEstimator.NormalizerMode.MeanVariance))
//     // Create the model
//     |> EstimatorModel.appendBy (fun mlc -> 
//                                     mlc.Regression.Trainers.FastTree(
//                                                                     labelColumn = DefaultColumnNames.Label,
//                                                                     featureColumn = DefaultColumnNames.Features,
//                                                                     numTrees=200,numLeaves=2000,minDatapointsInLeaves=1
                                                                    
//                                                                     ) )
//     |> EstimatorModel.appendCacheCheckpoint
//     // Train the model
//     |> EstimatorModel.fit trainTestSplit.TrainingData


// let metrics = 
//     model
//     |> Evaluation.Regression.InitEvaluate() trainTestSplit.TrainingData

// let preds = 
//     model
//     |> TransformerModel.transform trainTestSplit.TestData

// let yy = 
//     preds.GetColumn<float>(ctx,"CorrespondingScanTime")

// let yyHat = 
//     preds.GetColumn<float32>(ctx,"Score")

// [
// itemsForAlignment
// |> Array.map (fun x -> x.MeanScanTime-x.CorrespondingScanTime)
// |> Chart.BoxPlot;
// Seq.map2 (fun x y -> x - float y) yy yyHat
// |> Chart.BoxPlot
// ]
// |> Chart.Combine
// |> Chart.Show

// 5+5

// Chart.Point(yy,yyHat)
// |> Chart.Show

// //let permMetrix = model.Context.Regression.PermutationFeatureImportance(model.TransformerChain.LastTransformer |> FSharpML.Transformer.downcastTransformer,preds,permutationCount=5)
// //permMetrix.Length
// //preds

// //let x = 
// //    preds.Preview().ColumnView 
// //    |> Seq.map (fun x -> x.Column.Name)
// //    |> Array.ofSeq
// //    |> fun x -> x.[3..6]


// //x
// //|> Array.mapi (fun i x -> x, permMetrix.[i].RSquared.Mean)
// //|> Chart.Column
// //|> Chart.Show

// preds.Schema
// |> Seq.map (fun x -> x.Name)
// |> Array.ofSeq
// metrics.Rms

// metrics.L1 

// metrics.L2

// //Key #1	Key #2	Key #3	ScanTime_20180927 MS SH116sbpase001	AverageSequestScore_20180927 MS SH116sbpase001	
// //ExperimentalMass_20180927 MS SH116sbpase001	N14MZ_20180927 MS SH116sbpase001	N14Quant_20180927 MS SH116sbpase001	
// //N15MZ_20180927 MS SH116sbpase001	N15Quant_20180927 MS SH116sbpase001	N15Minus1MZ_20180927 MS SH116sbpase001	
// //N15Minus1Quant_20180927 MS SH116sbpase001	ScanTime_20180927 MS SH116sbpase003	AverageSequestScore_20180927 MS SH116sbpase003	
// //ExperimentalMass_20180927 MS SH116sbpase003	N14MZ_20180927 MS SH116sbpase003	N14Quant_20180927 MS SH116sbpase003	N15MZ_20180927 
// //MS SH116sbpase003	N15Quant_20180927 MS SH116sbpase003	N15Minus1MZ_20180927 MS SH116sbpase003	N15Minus1Quant_20180927 MS SH116sbpase003

// rawData.GetColumn ("N14Quant_20180927 MS SH116sbpase001") 

// (rawData.GetColumn ("N14Quant_20180927 MS SH116sbpase001"))
// let ctx = new MLContext()

