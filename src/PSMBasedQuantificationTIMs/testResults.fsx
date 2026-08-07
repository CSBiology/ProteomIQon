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
#r "nuget: FSharp.Stats, 0.4.0"
#r "nuget: BioFSharp, 2.0.0-beta4"
#r "nuget: BioFSharp.Mz, 0.1.5-beta"
#r "nuget: FSharpAux, 1.0.0"
#r "nuget: FSharpAux.IO, 1.0.0"
#r "nuget: ProteomIQon"

open System.IO
open System.Data.SQLite
open ProteomIQon
open ProteomIQon.Core
open ProteomIQon.Core.MzIO
open ProteomIQon.Dto
open FSharp.Stats
open BioFSharp.Mz.Quantification
open BioFSharp.Mz
open FSharpAux.IO.SchemaReader
open Plotly.NET
open BioFSharp
open MzIO.Processing
open BioFSharp.Mz.SearchDB

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

do fsi.AddPrinter(fun (printer:Deedle.Internal.IFsiFormattable) -> "\n" + (printer.Format()))

let inReader = new MzSQL(@"/home/paulinehans/Dokumente/testRunGabor/binned_spectra_1.000.mzlite")
inReader.Connection.Open()
let inRunID  = "sample=0"
let inTr = inReader.BeginTransaction()

//let ms2 = inReaderMS.ReadMassSpectrum("merged=377907 frame=45978 scanStart=1 scanEnd=918")

let ms = inReader.ReadMassSpectrum("merged=475568 frame=53975 scanStart=1 scanEnd=918")
ms |> MzIO.Processing.MassSpectrum.getScanTime



let unzipIMzliteArray (a:IMzIOArray<Peak1D>) = 
   let mzData = Array.zeroCreate a.Length
   let intensityData = Array.zeroCreate a.Length
   for i = 0 to a.Length-1 do 
       let peak = a.[i]
       mzData.[i] <- peak.Mz
       intensityData.[i] <- peak.Intensity
   mzData,intensityData



let initRTProfile (readspecPeaks:string -> Peak1DArray)  (rtIndex: IMzIOArray<RtIndexEntry>) (rtRange: RangeQuery) (mzRange: RangeQuery) (ionMobilityRange: RangeQuery)=
    let entries = RtIndexEntry.Search(rtIndex, rtRange)
    entries
    |>Seq.map (fun entry ->
        let peaks = (readspecPeaks entry.SpectrumID).Peaks
        let p = 
            peaks
            |> Seq.filter (fun x ->
                x.IonMobility.Value < ionMobilityRange.HighValue && x.IonMobility.Value > ionMobilityRange.LowValue &&
                x.Mz < mzRange.HighValue && x.Mz > mzRange.LowValue
            )
            |> fun s ->
                if (s |> Seq.tryHead).IsNone then
                    RtIndexEntry.AsPeak2D ((new Peak1D(0., mzRange.LockValue, ionMobilityRange.LockValue)),entry.Rt)
                else
                    s
                    |> Seq.groupBy (fun x -> x.Mz)
                    |> Seq.map (fun (mz,peaks) -> 
                        let i = peaks |> Seq.sumBy (fun x -> x.Intensity)
                        new Peak1D(i, mz, ionMobilityRange.LockValue)
                    )
                    |> fun x -> 
                        RtIndexEntry.ClosestMz (x, mzRange.LockValue)
                    |> fun x -> RtIndexEntry.AsPeak2D (x, entry.Rt)
        p
    )
    |> Array.ofSeq


let initRTProfile2
    (readspecPeaks: string -> Peak1DArray)
    (rtIndex: IMzIOArray<RtIndexEntry>)
    (rtRange: RangeQuery)
    (mzRange: RangeQuery)
    (ionMobilityRange: RangeQuery) =

    let entries = RtIndexEntry.Search(rtIndex, rtRange)
    entries
    |> Seq.collect (fun entry ->
        let peaks =
            (readspecPeaks entry.SpectrumID).Peaks
        // Keep only peaks inside the requested m/z and ion mobility window.
        let filteredPeaks =
            peaks
            |> Seq.filter (fun x ->
                x.IonMobility.IsSome
                &&
                x.IonMobility.Value < ionMobilityRange.HighValue
                &&
                x.IonMobility.Value > ionMobilityRange.LowValue
                &&
                x.Mz < mzRange.HighValue
                &&
                x.Mz > mzRange.LowValue
            )
            |> Seq.toArray
        if filteredPeaks.Length = 0 then
            // No signal was found for this RT.
            // Keep a zero-intensity placeholder at the expected position.
            Seq.singleton (
                RtIndexEntry.AsPeak2D(
                    new Peak1D(
                        0.,
                        mzRange.LockValue,
                        ionMobilityRange.LockValue
                    ),
                    entry.Rt
                )
            )
        else
            // IMPORTANT:
            // Treat every measured ion mobility separately.
            //
            // Before this change, all mobility values inside the search
            // window were collapsed into one value.
            filteredPeaks
            |> Seq.groupBy (fun peak ->
                peak.IonMobility.Value)
            |> Seq.map (fun (mobility, mobilityPeaks) ->
                // Within ONE mobility position there can still be multiple
                // peaks with the same m/z. Sum only those intensities.
                let mzPeaks =
                    mobilityPeaks
                    |> Seq.groupBy (fun peak -> peak.Mz)
                    |> Seq.map (fun (mz, sameMzPeaks) ->
                        let intensity =
                            sameMzPeaks
                            |> Seq.sumBy (fun peak -> peak.Intensity)
                        // Keep the ACTUAL measured ion mobility.
                        new Peak1D(intensity,mz,mobility
                        )
                    )
                // Find the m/z closest to the expected precursor m/z,
                // but only within this particular mobility position.
                let closestMzPeak =
                    RtIndexEntry.ClosestMz(
                        mzPeaks,
                        mzRange.LockValue
                    )
                // Add the RT coordinate.
                RtIndexEntry.AsPeak2D(
                    closestMzPeak,
                    entry.Rt
                )
            )
    )
    |> Array.ofSeq

let filePathPSMS = ("/home/paulinehans/Dokumente/testRunGabor/binned_spectra_1.000.mzlite.qpsm")
let psms = Frame.ReadCsv (filePathPSMS, hasHeaders =true, separators = "\t")
let im = 
    psms
    |> Frame.mapRows (fun j s -> 
        {|
            IonMobility = s.GetAs<float>("IonMobility")
            RTTraceLight = s.GetAs<float>("ScanTime")
            MassLight = s.GetAs<float>("PrecursorMZ")
            Hyper = s.GetAs<float>("Hyperscore")
        |}) 
        |> Series.values
        |> Seq.toArray
        |> Array.sortByDescending (fun x-> x.Hyper)


let retTimeIdxed = Query.getMS1RTIdx inReader inRunID
#time
let allQuerys =
    im 
    |> Array.take 5
    |> Array.map (fun r -> 
        let rt = RangeQuery (r.RTTraceLight, 2.0)
        let mz = RangeQuery (r.MassLight, 1.0)
        let im = RangeQuery (r.IonMobility, 0.05)
        initRTProfile2 inReader.ReadSpectrumPeaks retTimeIdxed  rt mz im
    )      
    
let transform  = 
    let waveletData = 
        allQuerys 
        |> Array.map (fun x -> 
            x   
            |> Array.filter (fun x -> 
                x.IonMobility.IsSome
            ) 
            |> Array.map (fun x -> x.Rt, x.IonMobility.Value, x.Intensity)  
        )
    waveletData
    |> Array.concat
    //|> Array.filter (fun (rt, im, i) -> i > 0.0)


transform

let visuPeak= 
    transform
    |> Chart.Point3D
    |> Chart.withMarkerStyle (Size = 1)
    |> Chart.saveHtml ("/home/paulinehans/Dokumente/testRunGabor/peakshape1000.html")
    
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

