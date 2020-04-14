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
let parseString s = 
    if s = "nan" then nan else float s

///
let fileNameOf path = FSharpAux.IO.PathFileName.fileNameOfPath path

///
let ctx = new MLContext()

///
type AlignmentResult = 
    {
    StringSequence               : string;
    GlobalMod                    : int;
    Charge                       : int;
    N14QuantMz                   : float
    N15QuantMz                   : float
    N15Minus1QuantMz             : float
    PredictedScanTime            : float
    ProteinNames                 : string;
    }  
    
type AlignmentItem<'a> = {
    Sequence             : string
    GlobalMod            : bool
    Charge               : int
    Mz                   : float          
    MeanScanTime         : 'a
    Stabw                : 'a
    Intensity            : 'a
    N14Quant             : 'a
    N15Quant             : 'a
    N14QuantMz           : 'a
    N15QuantMz           : 'a
    N15Minus1QuantMz     : 'a
    Ratio                : float
    }

type AlignmentItemForLearning<'a> = {
    Sequence             : string
    GlobalMod            : bool
    Charge               : int
    Mz                   : float          
    MeanScanTime         : 'a
    CorrespondingScanTime: 'a
    Stabw                : 'a
    Intensity            : 'a
    N14Quant             : 'a
    N15Quant             : 'a
    N14QuantMz           : 'a
    N15QuantMz           : 'a
    N15Minus1QuantMz     : 'a
    Ratio                : float
    }

let alignmentItemForLearningOf correspondingScanTime (aItem:AlignmentItem<float>) =
    {
    Sequence             = aItem.Sequence             
    GlobalMod            = aItem.GlobalMod            
    Charge               = aItem.Charge               
    Mz                   = aItem.Mz                            
    MeanScanTime         = aItem.MeanScanTime         
    CorrespondingScanTime= correspondingScanTime
    Stabw                = aItem.Stabw                
    Intensity            = aItem.Intensity            
    N14Quant             = aItem.N14Quant             
    N15Quant             = aItem.N15Quant             
    N14QuantMz           = aItem.N14QuantMz           
    N15QuantMz           = aItem.N15QuantMz           
    N15Minus1QuantMz     = aItem.N15Minus1QuantMz     
    Ratio                = aItem.Ratio                
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

type AlignmentFile<'a> = {
    FileName                   : string
    Peptides                   : array<AlignmentItem<'a>> 
    PeptidesSharedWithAllFiles : array<AlignmentItem<'a>> 
    MissingPeptides            : array<Peptide>
    }
    
///
let dataPath = Path.Combine [|__SOURCE_DIRECTORY__; "../data/" |]

///
let files = 
    FSharpAux.IO.FileIO.filesInDir (FSharpAux.IO.FileIO.directoryInfo dataPath)
    |> Array.filter (fun file -> file.Extension = ".quant")
    |> Array.map (fun file -> file.FullName)


///
let createModel (ctx:MLContext) dataForFitAndEval =
    let trainTestSplit = 
        ctx.Data.ReadFromEnumerable(dataForFitAndEval)
        |> FSharpML.Data.Regression.initTrainTestSplit(ctx,Testfraction=0.1)
   
    let model = 
        EstimatorModel.create ctx
        // Process data transformations in pipeline
        |> EstimatorModel.appendBy (fun mlc -> mlc.Transforms.Conversion.ConvertType(DefaultColumnNames.Label,"CorrespondingScanTime",DataKind.R4))
        |> EstimatorModel.appendBy (fun mlc -> 
                                        mlc.Transforms.Concatenate(
                                                                    "FeaturesRaw" , 
                                                                    "Mz", 
                                                                    "MeanScanTime", 
                                                                    "Stabw", 
                                                                    "Intensity"
                                                                ) )
        |> EstimatorModel.appendBy (fun mlc -> mlc.Transforms.Conversion.ConvertType(DefaultColumnNames.Features,"FeaturesRaw",DataKind.R4))
        //|>  EstimatorModel.appendBy (fun mlc -> mlc.Transforms.Normalize(DefaultColumnNames.Features,DefaultColumnNames.Features,Microsoft.ML.Transforms.Normalizers.NormalizingEstimator.NormalizerMode.MeanVariance))
        // Create the model
        |> EstimatorModel.appendBy (fun mlc -> 
                                        mlc.Regression.Trainers.FastTree(
                                                                        labelColumn = DefaultColumnNames.Label,
                                                                        featureColumn = DefaultColumnNames.Features(*,*)
                                                                        //numTrees=200,numLeaves=2000,minDatapointsInLeaves=1
                                                                    
                                                                        ) )
        |> EstimatorModel.appendCacheCheckpoint
        // Train the model
        |> EstimatorModel.fit trainTestSplit.TrainingData

    let preds = 
        model
        |> TransformerModel.transform trainTestSplit.TestData

    let y = 
        preds.GetColumn<float>(ctx,"CorrespondingScanTime")

    let yHat = 
        preds.GetColumn<float32>(ctx,"Score")

    let evalPlot = 
        Seq.map2 (fun y yHat -> y - float yHat) y yHat
        |> Chart.BoxPlot
    evalPlot, 
    (fun (data:AlignmentItemForLearning<float> []) -> 
        let dataPrep = ctx.Data.ReadFromEnumerable(data)
        let transFormedData = model |> TransformerModel.transform dataPrep
        let yHat = transFormedData.GetColumn<float32>(ctx,"Score")
         
        Seq.map2 (fun (alignmentItem:AlignmentItemForLearning<float> ) predScanTime -> 
                        {
                        StringSequence     = alignmentItem.Sequence
                        GlobalMod          = if alignmentItem.GlobalMod then 1 else 0
                        Charge             = alignmentItem.Charge
                        N14QuantMz         = alignmentItem.N14QuantMz
                        N15QuantMz         = alignmentItem.N15QuantMz
                        N15Minus1QuantMz   = alignmentItem.N15Minus1QuantMz
                        PredictedScanTime  = float predScanTime
                        ProteinNames       = alignmentItem.Sequence
                        }  
                  ) data yHat
    )

///
let getAlignmentItems (dataPath:string) = 
    Frame.ReadCsv(dataPath,separators="\t") 
    |> Frame.mapRowValues (fun x -> 
                            let sequence     = x.GetAs<string>("StringSequence")
                            let globalMod    = x.GetAs<bool>("GlobalMod")
                            let charge       = x.GetAs<int>("Charge")
                            let precMZ       = x.GetAs<float>("PrecursorMZ")
                            let N14Quant     = x.GetAs<string>("N14Quant")
                            let N15Quant     = x.GetAs<string>("N15Quant")
                            let N14QuantMz     = x.GetAs<string>("N14QuantMz")
                            let N15QuantMz     = x.GetAs<string>("N15QuantMz")
                            let N15Minus1QuantMz = x.GetAs<string>("N15Minus1QuantMz")
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
                                N14Quant = N14Quant 
                                N15Quant = N15Quant
                                N14QuantMz = N14QuantMz 
                                N15QuantMz = N15QuantMz
                                N15Minus1QuantMz = N15Minus1QuantMz
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
                        N14Quant = parseString x.N14Quant 
                        N15Quant = parseString x.N15Quant 
                        N14QuantMz = parseString x.N14QuantMz
                        N15QuantMz = parseString x.N15QuantMz
                        N15Minus1QuantMz = parseString x.N15Minus1QuantMz
                        Ratio = parseString x.N14Quant / parseString x.N15Quant
                        } |> Some
                        with 
                        | _-> 
                            printfn "%s x.N14Quant, %s x.N15Quant" x.N14Quant x.N15Quant
                            None
                 )
    |> Array.ofSeq

///
let allAlignmentItems = 
    files
    |> Seq.map (fun x -> 
                    {
                    FileName                    = fileNameOf x
                    Peptides                    = getAlignmentItems x  
                    PeptidesSharedWithAllFiles  = Array.empty
                    MissingPeptides             = Array.empty
                    }
               )

///
let allIdentifiedPeptides = 
    allAlignmentItems
    |> Seq.map (fun x -> x.Peptides)
    |> Seq.concat
    |> Seq.map createPeptideof
    |> Seq.distinct
    |> Array.ofSeq

///    
let allSharedPeptides =
    allAlignmentItems
    |> Seq.map (fun x -> x.Peptides)
    |> Seq.map (Seq.map createPeptideof)
    |> Seq.map (Set.ofSeq)
    |> Set.intersectMany

///
let filteredForSharedAligmentItems = 
    allAlignmentItems
    |> Seq.map (fun quantFile -> 
                    let sharedPeptides  = 
                        quantFile.Peptides 
                        |> Array.filter (fun item -> allSharedPeptides.Contains (createPeptideof item)) 
                    let allMissingPeptides =
                        allIdentifiedPeptides
                        |> Array.choose (fun identifiedPep -> 
                                        match quantFile.Peptides |> Seq.tryFind (fun peptideInFile -> createPeptideof peptideInFile = identifiedPep) with 
                                        | Some x -> None 
                                        | None   -> Some identifiedPep
                                   )
                    {quantFile with PeptidesSharedWithAllFiles = sharedPeptides ;MissingPeptides = allMissingPeptides}
               )
    |> Array.ofSeq

///
let alignmentRanking = 
    [|
    for i = 0 to filteredForSharedAligmentItems.Length-1 do         
        for j = 0 to filteredForSharedAligmentItems.Length-1 do 
            printfn "%i %i" i j 
            if j = i then 
                () 
            else
                let x = filteredForSharedAligmentItems.[i].PeptidesSharedWithAllFiles
                printfn "peptides = %i" x.Length
                let y = filteredForSharedAligmentItems.[j].PeptidesSharedWithAllFiles
                let (x',y') = 
                    x
                    |> Array.map (fun x -> x , Array.find (fun y -> createPeptideof x = createPeptideof y ) y)
                    |> Array.filter (fun x -> nan.Equals((fst x).MeanScanTime) |> not && nan.Equals((snd x).MeanScanTime) |> not)
                    |> Array.filter (fun x -> (fst x).Ratio < 5. && (fst x).Ratio > 0.01  && (snd x).Ratio < 5. && (snd x).Ratio > 0.01 )
                    |> Array.unzip
                printfn "peptides = %i" x'.Length
                yield (filteredForSharedAligmentItems.[i],filteredForSharedAligmentItems.[j], Array.map2 (fun x y -> abs (x-y) ) (x' |> Array.map (fun x -> x.MeanScanTime)) (y'|> Array.map (fun x -> x.MeanScanTime)) |> Seq.mean) 
    |]
    |> Array.groupBy (fun (t,s,dis) -> t)  
    |> Array.map (fun (t,ranking) -> t,ranking |> Array.map (fun (t,s,dis) -> s,dis) |> Array.sortBy (fun (s,dis) -> dis) )
    

///
let performAlignments (targetFile:AlignmentFile<'a>) (sourceFiles:AlignmentFile<'a> []) =
    let missingPeptides = targetFile.MissingPeptides
    let rec loop targetmissingPeptides counter pepsToRequantify qualityPlots = 
        if counter = sourceFiles.Length then
            pepsToRequantify,qualityPlots |> Chart.Combine
        else 
            let sourceFile     = sourceFiles.[counter]
            let sharedPeptides = 
                targetFile.Peptides
                |> Array.choose (fun targetPep -> 
                                    match sourceFile.Peptides |> Array.tryFind (fun sourcePep -> createPeptideof sourcePep = createPeptideof targetPep) with
                                    | Some sourcePep -> 
                                        alignmentItemForLearningOf targetPep.MeanScanTime sourcePep 
                                        |> Option.Some 
                                    | None -> None )
            let sourceFileSpecificPeptides = 
                sourceFile.Peptides
                |> Array.filter (fun sourcePep -> 
                                    targetmissingPeptides |> Array.exists (fun targetPep -> targetPep = createPeptideof sourcePep )
                                )
            let sourceFileSpecPeptidesForPred =
                sourceFileSpecificPeptides 
                |> Array.map (alignmentItemForLearningOf nan)
            let plot,model = createModel ctx sharedPeptides
            let alignmentResults = model (sourceFileSpecPeptidesForPred )
            let targetMissingPeptides' =
                targetmissingPeptides
                |> Array.filter (fun x -> sourceFileSpecificPeptides |> Array.exists (fun sourceFilePeptide -> createPeptideof sourceFilePeptide = x) |> not)
            printfn "still missing peptides: %i" (targetMissingPeptides'.Length)
            loop targetMissingPeptides' (counter+1) (Seq.append alignmentResults pepsToRequantify) (Seq.append [plot] qualityPlots) 
    printfn "still missing peptides: %i" (missingPeptides.Length)
    loop missingPeptides 0 Seq.empty Seq.empty

#time
let res = 
    alignmentRanking 
    |> Array.map (fun (targetFile,sourceFiles) -> 
                    let alignments,chart = performAlignments targetFile (sourceFiles |> Array.map fst)
                    alignments
                    |> FSharpAux.IO.SeqIO.Seq.toCSV "\t" true
                    |> FSharpAux.IO.SeqIO.Seq.write 
                    
                 )


res.[0] |> fst     |> Seq.length                 
//             )

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

