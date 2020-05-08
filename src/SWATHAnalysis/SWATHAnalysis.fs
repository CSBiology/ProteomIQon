namespace ProteomIQon

open Argu
open System
open System.IO
open ProteomIQon
open BioFSharp
open MzIO.IO
open MzIO.Processing
open MzIO.Processing.MzIOLinq
open MzIO.MetaData.PSIMSExtension
open MzIO.MetaData.ParamEditExtension
open MzIO.MetaData.UO.UO
open MzIO.Binary
open MzIO.Model
open MzIO.Model.CvParam
open MzIO.Commons.Arrays
open System.Linq
open FSharpAux
open FSharpAux.IO

module SWATHAnalysis = 
  
    type LibraryEntry =
        {
            Charge: float
            Iontype: string
            MassOverCharge: float
            Number: int
            Intensity: float
            ModSequenceID: int
            PSMId: string
            PrecursorMZ: float
            ScanTime: float
            Count: float
            Version: int
            Sequence: string
        }

    let createPeptideQuery (peptide:string) (library: LibraryEntry[]) matchingTolerance=
        let entries =
            library
            |> Array.filter (fun entry -> entry.Sequence = peptide)
        let offset, averageRT =
            let upperRT =
                entries
                |> Array.maxBy (fun entry -> entry.ScanTime)
                |> (fun entry -> entry.ScanTime)
            let lowerRT =
                entries
                |> Array.minBy (fun entry -> entry.ScanTime)
                |> (fun entry -> entry.ScanTime)
            let averageRT =
                entries
                |> Array.averageBy (fun entry -> entry.ScanTime)
            let offset =
                if abs (averageRT-upperRT) > abs (averageRT-lowerRT) then
                    abs (averageRT-upperRT) + 10.
                else
                    abs (averageRT-lowerRT) + 10.
            offset, averageRT
        let rangeQueryRt = MzIO.Processing.Query.createRangeQuery averageRT offset
        let rangeQueriesMz =
            entries
            |> Array.map (fun entry ->
                MzIO.Processing.Query.createRangeQuery entry.MassOverCharge (Mass.deltaMassByPpm matchingTolerance entry.MassOverCharge)
            )
        let swathQuery = MzIO.Processing.Query.createSwathQuery entries.[0].PrecursorMZ rangeQueryRt rangeQueriesMz
        averageRT, swathQuery

    let getClosestMz (peaks: Peak1DArray, mzRange: RangeQuery) =
        let p1d = 
            peaks.Peaks
                .DefaultIfEmpty(new Peak1D(0., mzRange.LockValue))
                .ItemAtMin(fun x -> Math.Abs(x.Mz - mzRange.LockValue))
        if p1d.Mz > mzRange.HighValue || p1d.Mz < mzRange.LowValue then
            new Peak1D(0., mzRange.LockValue)
        else
            p1d

    let getRTProfiles (swathIndexer: SwathIndexer.SwathIndexer) (reader: IMzIODataReader) (swathQuery: SwathQuery) =
        swathIndexer.GetRTProfiles(reader, swathQuery, false, getClosestMz)

    let quantify (library: string) (swathAnalysisParams: Domain.SWATHAnalysisParams) (instrumentOutput: string) (outFolder: string)=
        
        let consensusLibrary =
            Seq.fromFileWithCsvSchema<LibraryEntry>(library, '\t', true,schemaMode = FSharpAux.IO.SchemaReader.Csv.SchemaModes.Fill)
            |> Seq.toArray

        let inReader = Core.MzIO.Reader.getReader instrumentOutput
        let tr = inReader.BeginTransaction()
        let runID = Core.MzIO.Reader.getDefaultRunID inReader
        let swathIdx = SwathIndexer.SwathIndexer.Create(inReader, runID)

        let peptideList =
            match swathAnalysisParams.PeptideList with
            | Some x -> x
            | None ->
                consensusLibrary
                |> Array.map (fun x -> x.Sequence)
                |> Array.distinct

        let queries =
            peptideList
            |> Array.map (fun peptide -> peptide, createPeptideQuery peptide consensusLibrary swathAnalysisParams.MatchingTolerancePPM)

        let quant =
            queries
            |> Array.map (fun (peptide,(rt,query)) ->
                let rtProfiles =
                    getRTProfiles swathIdx inReader query
                    |> Array2D.toJaggedArray
                    |> Array.map (fun peak2DArr -> 
                        peak2DArr
                        |> Array.map (fun peak2D ->
                            peak2D.Rt, peak2D.Intensity
                        )
                    )
                let quants = 
                    rtProfiles
                    |> Array.map (fun retInt ->
                        retInt
                        |> Array.unzip
                        |> fun (ret, intensity) ->
                            FSharp.Stats.Signal.PeakDetection.SecondDerivative.getPeaks 0.1 2 13 ret intensity
                        |> fun peaks ->
                            BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks rt
                        |> BioFSharp.Mz.Quantification.HULQ.quantifyPeak
                        |> fun x -> x.Area
                    )
                {|Peptide = peptide; Area = (quants |> Array.sum)|}
            )
        FSharpAux.IO.SeqIO.Seq.CSV "\t" true false quant
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFolder + "/" + (System.IO.Path.GetFileNameWithoutExtension instrumentOutput) + ".swath")