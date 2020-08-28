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
open FSharp.Stats

module SWATHAnalysis = 

    type SwathIndexer.SwathIndexer with
        
        member this.GetRTProfiles2(dataReader:IMzIODataReader, query: SwathQuery, getLockMz: bool, ?mzRangeSelector: Peak1DArray * RangeQuery -> Peak1D) =

            let getClosestMz (peaks: Peak1DArray, mzRange: RangeQuery) =
                peaks.Peaks
                    .DefaultIfEmpty(new Peak1D(0., mzRange.LockValue))
                    .ItemAtMin(fun x -> Math.Abs(x.Mz - mzRange.LockValue))

            let mzRangeSelector = defaultArg mzRangeSelector getClosestMz

            let swathSpectra = 
                this.SwathList.SearchAllTargetMz(query.TargetMz)
                    .Take(1)
                    .SelectMany(fun x -> x.SearchAllRt(query))
                    .ToArray()
            if swathSpectra.Length > 0 then

                let profile = Array2D.create query.CountMS2Masses swathSpectra.Length (new Peak2D())

                for specIdx = 0 to swathSpectra.Length - 1 do

                    let swathSpec = swathSpectra.[specIdx]
                    let pa = dataReader.ReadSpectrumPeaks(swathSpec.SpectrumID)

                    for ms2MassIndex = 0 to query.CountMS2Masses - 1 do
                        
                        let mzRange = query.Ms2Masses.[ms2MassIndex]
                        let p = mzRangeSelector(pa, mzRange)

                        if getLockMz then
                            profile.[ms2MassIndex, specIdx] <- new Peak2D(p.Intensity, mzRange.LockValue, swathSpec.Rt)
                        else
                            profile.[ms2MassIndex, specIdx] <- new Peak2D(p.Intensity, p.Mz, swathSpec.Rt)

                Some profile
            else
                None
  
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

    let createPeptideQuery (peptide:string) (library: LibraryEntry[]) matchingTolerance offsetRange =
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
                    abs (averageRT-upperRT) + offsetRange
                else
                    abs (averageRT-lowerRT) + offsetRange
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
        swathIndexer.GetRTProfiles2(reader, swathQuery, false, getClosestMz)

    let optimizeWindowWidth polOrder (windowWidthToTest:int[]) noiseAutoCorr (signalOfInterest:float[]) =
        let signalOfInterest' = signalOfInterest |> vector
        //let noiseAutoCorr = Correlation.Vector.autoCorrelation 1 (blankSignal |> vector)
        let filterF w yData = FSharp.Stats.Signal.Filtering.savitzky_golay w polOrder 0 0 yData
        let windowWidthToTest' = windowWidthToTest |> Array.filter (fun x -> x%2 <> 0)
        let optimizedWindowWidth =
            windowWidthToTest'
            |> Array.map (fun w ->
                let smoothedY = filterF w signalOfInterest
                let noise = smoothedY - (signalOfInterest')
                w, Correlation.Vector.autoCorrelation 1 noise
                )
            |> Array.minBy (fun (w,ac) -> (ac - noiseAutoCorr) |> abs )
            |> fst
        optimizedWindowWidth

    ///
    let initGetWindowWidth (windowEst:Domain.WindowSize) polynomOrder (windowWidthToTest:int[]) =
        match windowEst with
        | Domain.WindowSize.Fixed w  -> fun yData -> w
        | Domain.WindowSize.EstimateUsingAutoCorrelation noiseAutoCorr -> fun yData -> optimizeWindowWidth polynomOrder windowWidthToTest noiseAutoCorr yData


    ///
    let initIdentifyPeaks (peakDetectionParams:Domain.XicProcessing) =
        match peakDetectionParams with 
        | Domain.XicProcessing.SecondDerivative parameters ->
            let getWindowWidth = initGetWindowWidth parameters.WindowSize parameters.PolynomOrder [|5 .. 2 .. 60|] 
            (fun xData yData -> 
                let  windowSize = getWindowWidth yData
                FSharp.Stats.Signal.PeakDetection.SecondDerivative.getPeaks parameters.MinSNR parameters.PolynomOrder windowSize xData yData
                )
        | Domain.XicProcessing.Wavelet parameters ->
            (fun xData yData -> 
                FSharpStats'.Wavelet.identify parameters xData yData
                )

    let quantify (library: string) (swathAnalysisParams: Domain.SWATHAnalysisParams) (instrumentOutput: string) (outFolder: string) =
        let outFilePath = outFolder + "/" + (System.IO.Path.GetFileNameWithoutExtension instrumentOutput) + ".squant"
        let consensusLibrary =
            Seq.fromFileWithCsvSchema<LibraryEntry>(library, '\t', true,schemaMode = FSharpAux.IO.SchemaReader.Csv.SchemaModes.Fill)
            |> Seq.toArray

        let inReader = Core.MzIO.Reader.getReader instrumentOutput
        let tr = inReader.BeginTransaction()
        let runID = Core.MzIO.Reader.getDefaultRunID inReader
        let swathIdx = SwathIndexer.SwathIndexer.Create(inReader, runID)
        let identifyPeaks = initIdentifyPeaks swathAnalysisParams.XicProcessing
        let peptideList =
            match swathAnalysisParams.PeptideList with
            | Some x -> x
            | None ->
                consensusLibrary
                |> Array.map (fun x -> x.Sequence)
                |> Array.distinct
        let queries =
            peptideList
            |> Array.map (fun peptide -> peptide, createPeptideQuery peptide consensusLibrary swathAnalysisParams.MatchingTolerancePPM swathAnalysisParams.QueryOffsetRange)
        let quant =
            queries
            |> Array.choose (fun (peptide,(rt,query)) ->
                let profile = getRTProfiles swathIdx inReader query
                // Match to look if RTProfile was found. If none was found, quantification is skipped.
                match profile with
                | Some rtProfile ->
                    let rtProfiles =
                        rtProfile
                        |> Array2D.toJaggedArray
                        |> Array.map (fun peak2DArr ->
                            peak2DArr
                            |> Array.map (fun peak2D ->
                                peak2D.Rt, peak2D.Intensity
                            )
                        )
                    let quants = 
                        rtProfiles
                        |> Array.choose (fun retInt ->
                            retInt
                            |> Array.unzip
                            |> fun (ret, intensity) ->
                                identifyPeaks ret intensity
                            |> fun peaks ->
                                match peaks with
                                | [||] -> None
                                | _    -> Some (BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks rt)
                        )
                        |> Array.map (fun peak ->
                            BioFSharp.Mz.Quantification.HULQ.quantifyPeak peak
                            |> fun x -> x.Area
                        )
                    Some {|Peptide = peptide; Area = (quants |> Array.sum)|}
                | None -> None
            )
        tr.Dispose()
        if quant.Length >= 1 then
            FSharpAux.IO.SeqIO.Seq.CSV "\t" true false quant
            |> FSharpAux.IO.SeqIO.Seq.writeOrAppend outFilePath
        else
            ()