#r "nuget: MzIO, 0.1.7"
#r "nuget: MzIO.Processing, 0.1.7"
#r "nuget: MzLite, 1.1.0"
#r "nuget: MzIO.SQL, 0.1.9"
#r "nuget: Plotly.NET, 5.0.0"
#r "nuget: Plotly.NET.Interactive, 5.0.0"
#r "nuget: MzIO.MzML, 0.1.10"
#r "nuget: Deedle, 2.5.0"
#r "nuget: BioFSharp.ImgP, 2.0.0-beta5"
#r "nuget: BioFSharp.Mz, 0.1.5-beta"
#r "nuget: FsMath, 0.0.2"
#r "nuget: FSharpAux, 2.1.0"
#r "nuget: ProteomIQon, 0.0.10"
#r "nuget: BioFSharp.Mz, 0.1.5-beta"

open Plotly.NET
open System
open FSharpAux 
open BioFSharp.ImgP
open BioFSharp.Mz.Quantification
open BioFSharp.Mz
open System
open FsMath
open FSharpAux
open FSharp.Stats
open MzIO
open MzIO.Processing
open MzLite
open MzIO.MzSQL
open System.Data
// open Plotly.NET
// open Plotly.NET.TraceObjects
// open Plotly.NET.LayoutObjects
// open Plotly.NET.Interactive
open MzIO.Commons.Arrays
open MzIO.Processing.MzIOLinq
open MzIO.Binary
open Deedle
open System.Data.SQLite
open System.Linq
open MzIO.IO.MzML
open FSharp.Stats.Signal
open System.IO

let filePathQuantGabor = ("/home/paulinehans/Dokumente/testRunGabor/GaborOutput/binned_spectra_1.000.quant")
let quantFileGabor = Frame.ReadCsv (filePathQuantGabor, hasHeaders =true, separators = "\t")

let filePathOther = ("/home/paulinehans/Dokumente/PaulineTIMSDaten/runs/TIMsTest/out/singleFilequant/binned_spectra_1.000.quant")
let quantFileOther = Frame.ReadCsv (filePathOther, hasHeaders =true, separators = "\t")

let quantFileFP = Frame.ReadCsv ("/home/paulinehans/Dokumente/PaulineTIMSDaten/runs/TIMsTest/FragPipe/14N2/combined_ion.tsv", hasHeaders =true, separators = "\t")


let quantifiedRangeG: Frame<string*string,string> = 
    quantFileGabor
    |> Frame.indexRowsUsing (fun s->
        s.GetAs<string> "StringSequence",
        s.GetAs<string> "Charge"
    )
    |> Frame.sliceCols ["Quant_Light"]
    |> Frame.mapColKeys (fun s -> s + "_Gabor")

quantifiedRangeG    
let quantifiedRangeO: Frame<string*string,string> = 
    quantFileOther
    |> Frame.indexRowsUsing (fun s->
        s.GetAs<string> "StringSequence",
        s.GetAs<string> "Charge"
    )
    |> Frame.sliceCols ["Quant_Light"]
    |> Frame.mapColKeys (fun s -> s + "_Other")
quantifiedRangeO.RowCount


let quantifiedRangeFP : Frame<string * string, string> =
    quantFileFP
    |> Frame.aggregateRowsBy
        ["Peptide Sequence"; "Charge"]
        ["1 Intensity"]
        Stats.mean
    |> Frame.indexRowsUsing (fun row ->
        row.GetAs<string> "Peptide Sequence",
        row.GetAs<string> "Charge"
    )
    |> Frame.sliceCols ["1 Intensity"]
    |> Frame.mapColKeys (fun _ -> "Intensity_FP")

quantifiedRangeFP.RowCount
    

// let quantifiedRangeFP: Frame<string*string,string> = 
//     quantFileFP
//     |> Frame.indexRowsUsing (fun s->
//         s.GetAs<string> "Peptide Sequence",
//         s.GetAs<string> "Charge",
//         s.GetAs<string> "1 Apex Retention Time",
//         s.GetAs<string> "Modified Sequence"
//     )
//     |> Frame.sliceCols ["1 Intensity"]
//     |> Frame.mapValues float
//     |> Frame.applyLevel (fun (s,c,_,_) -> s,c) Stats.mean
// quantifiedRangeFP.RowCount
let both = Frame.join JoinKind.Inner quantifiedRangeG quantifiedRangeO
let both2 = Frame.join JoinKind.Inner quantifiedRangeG quantifiedRangeFP

let gabor: float [] = both2 |> Frame.getCol "Quant_Light_Gabor" |> Series.valuesAll |> Seq.toArray |> Array.map (fun x -> match x with | Some v -> log2 v | None -> 0.0)
let other: float [] = both |> Frame.getCol "Quant_Light_Other" |> Series.valuesAll |> Seq.toArray |> Array.map (fun x -> match x with | Some v -> log2 v | None -> 0.0)
let fp: float [] = both2 |> Frame.getCol "Intensity_FP" |> Series.valuesAll |> Seq.toArray |> Array.map (fun x -> match x with | Some v -> log2 v | None -> 0.0)

gabor.Length
other.Length
fp.Length

Array.zip gabor fp
|> Chart.Point
|> Chart.withXAxisStyle ("Gabor", ShowGrid = true)
|> Chart.withYAxisStyle ("Other", ShowGrid = true)
|> Chart.saveHtml ("/home/paulinehans/Dokumente/testRunGabor/quantificationComparisonWithFP.html")

// let pointOther = 
//     quantifiedRangeO 
//     |> Array.mapi (fun i x -> (i, x))
//     |> Chart.Point
//     //|> Chart.saveHtml ("/home/paulinehans/Dokumente/testRunGabor/otherLOG.html")

// let pointGabor = 
//     quantifiedRangeG 
//     |> Array.mapi (fun i x -> (i, x))
//     |> Chart.Point
//     //|> Chart.saveHtml ("/home/paulinehans/Dokumente/testRunGabor/gaborLOG.html")

// let visuSimilarity =  
//     [
//         pointOther |> Chart.withTraceInfo "Other"
//         pointGabor |> Chart.withTraceInfo "Gabor"
//     ]
//     |> Chart.combine
//     |> Chart.withTitle "Quantification results"
//     |> Chart.withXAxisStyle ("PSM index", ShowGrid = true)
//     |> Chart.withYAxisStyle ("Quantification value", ShowGrid = true)
//     |> Chart.saveHtml ("/home/paulinehans/Dokumente/testRunGabor/quantificationComparison.html")



let inReader = new MzSQL(@"")
inReader.Connection.Open()
let inRunID  = "sample=0"
let mzlite = "/home/paulinehans/Dokumente/testRunGabor/binned_spectra_1.000.mzlite"

