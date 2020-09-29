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

    type LibraryEntry =
        {
            Charge         : float
            Iontype        : string
            MassOverCharge : float
            Number         : int
            Intensity      : float
            PepSequenceID  : int
            ModSequenceID  : int
            PSMId          : string
            PrecursorMZ    : float
            ScanTime       : float
            Count          : float
            Version        : int
            Sequence       : string
            GlobalMod      : int
            PercolatorScore: float
        }

    type LibraryEntryQuant =
        {
            Charge         : float
            Iontype        : string
            MassOverCharge : float
            Number         : int
            Intensity      : float
            PepSequenceID  : int
            ModSequenceID  : int
            PSMId          : string
            PrecursorMZ    : float
            ScanTime       : float
            Count          : float
            Version        : int
            Sequence       : string
            GlobalMod      : int
            PercolatorScore: float
            Quant          : float
        }

    /// Contains the three dimensional information of a peak, consisting of one intensity and the corresponding m/z amd retention time value.
    type Peak2DInfo(intensity:float, mz:float, rt:float, fragmentInfo: LibraryEntry) =
    
        inherit Peak1D(intensity, mz)

        new () = Peak2DInfo(0., 0., 0.,
                    {
                        Charge         = -1.
                        Iontype        = ""
                        MassOverCharge = -1.
                        Number         = -1
                        Intensity      = -1.
                        PepSequenceID  = -1
                        ModSequenceID  = -1
                        PSMId          = ""
                        PrecursorMZ    = -1.
                        ScanTime       = -1.
                        Count          = -1.
                        Version        = -1
                        Sequence       = ""
                        GlobalMod      = -1
                        PercolatorScore= -1.
                    }
        )
    
        member this.Rt
            with get() = rt

        member this.FragmentInfo
            with get() = fragmentInfo
    
        override this.ToString() =
            String.Format("intensity={0}, mz={1}, rt={2}, fragment Info={3}", this.Intensity, this.Mz, this.Rt, this.FragmentInfo)

    type SwathIndexer.SwathIndexer with
        
        member this.GetRTProfiles2(dataReader:IMzIODataReader, query: SwathQuery, entries: LibraryEntry[], getLockMz: bool, spectrumSelector: seq<SwathIndexer.MSSwath> -> seq<SwathIndexer.MSSwath> list ,?mzRangeSelector: Peak1DArray * RangeQuery -> Peak1D) =

            let getClosestMz (peaks: Peak1DArray, mzRange: RangeQuery) =
                peaks.Peaks
                    .DefaultIfEmpty(new Peak1D(0., mzRange.LockValue))
                    .ItemAtMin(fun x -> Math.Abs(x.Mz - mzRange.LockValue))

            let mzRangeSelector = defaultArg mzRangeSelector getClosestMz

            let swathSpectra = 
                this.SwathList.SearchAllTargetMz(query.TargetMz)
                    |> spectrumSelector

            swathSpectra
            |> List.map (fun spec ->
                let swathSpectrum = spec.SelectMany(fun x -> x.SearchAllRt(query)).ToArray()

                if swathSpectrum.Length > 0 then

                    let profile = Array2D.create query.CountMS2Masses swathSpectrum.Length (new Peak2DInfo())

                    for specIdx = 0 to swathSpectrum.Length - 1 do

                        let swathSpec = swathSpectrum.[specIdx]
                        let pa = dataReader.ReadSpectrumPeaks(swathSpec.SpectrumID)

                        for ms2MassIndex = 0 to query.CountMS2Masses - 1 do
                        
                            let mzRange = query.Ms2Masses.[ms2MassIndex]
                            let p = mzRangeSelector(pa, mzRange)

                            if getLockMz then
                                profile.[ms2MassIndex, specIdx] <- new Peak2DInfo(p.Intensity, mzRange.LockValue, swathSpec.Rt, entries.[ms2MassIndex])
                            else
                                profile.[ms2MassIndex, specIdx] <- new Peak2DInfo(p.Intensity, p.Mz, swathSpec.Rt,  entries.[ms2MassIndex])

                    Some profile
                else
                    None
            )

    type PeptideInformation =
        {
            PepSequenceID  : int
            ModSequenceID  : int
            PrecursorMZ    : float
            Sequence       : string
            GlobalMod      : int
            PercolatorScore: float
        }

    let createPeptideInformation pepSeqID modSeqID precMz seq globalMod percolatorScore =
        {
            PepSequenceID   = pepSeqID
            ModSequenceID   = modSeqID
            PrecursorMZ     = precMz
            Sequence        = seq
            GlobalMod       = globalMod
            PercolatorScore = percolatorScore
        }

    type PeptideInformationQuant =
        {
            PepSequenceID  : int
            ModSequenceID  : int
            PrecursorMZ    : float
            StringSequence : string
            GlobalMod      : int
            PercolatorScore: float
            Quantification : float
        }

    let createPeptideInformationQuant pepSeqID modSeqID precMz seq globalMod percolatorScore quant=
        {
            PepSequenceID   = pepSeqID
            ModSequenceID   = modSeqID
            PrecursorMZ     = precMz
            StringSequence  = seq
            GlobalMod       = globalMod
            PercolatorScore = percolatorScore
            Quantification  = quant
        }


    type PeptideInformationQuantJoined =
        {
            PepSequenceID  : int
            PrecursorMZ_L  : float
            PrecursorMZ_H  : float 
            StringSequence : string
            PercolatorScore: float
            PercScoreL     : float
            PercScoreH     : float
            Quant_Light    : float
            Quant_Heavy    : float
        }

    let createPeptideInformationQuantJoined pepSeqID precMzL precMzH stringSeq percScore percScoreL percScoreH quantL quantH=
        {
            PepSequenceID  = pepSeqID
            PrecursorMZ_L  = precMzL
            PrecursorMZ_H  = precMzH
            StringSequence = stringSeq
            PercolatorScore= percScore
            PercScoreL     = percScoreL
            PercScoreH     = percScoreL
            Quant_Light    = quantL
            Quant_Heavy    = quantH
        }

    let parseSpectrumSelection (swathParam: Domain.SWATHAnalysisParams) =
        match swathParam.SpectrumSelectionF with
        |Domain.SpectrumSelection.All -> (fun x -> x |> Seq.fold (fun acc y -> seq[y]::acc)[])
        |Domain.SpectrumSelection.First -> (fun x -> [x.Take(1)])

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

    let createPeptideQuery (entries: LibraryEntry[]) matchingTolerance offsetRange =

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

    let getRTProfiles (swathIndexer: SwathIndexer.SwathIndexer) (reader: IMzIODataReader) (spectrumSelector: seq<SwathIndexer.MSSwath> -> seq<SwathIndexer.MSSwath> list) (entries: LibraryEntry[]) (swathQuery: SwathQuery) =
        swathIndexer.GetRTProfiles2(reader, swathQuery, entries, false, spectrumSelector ,getClosestMz)

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
            |> Array.minBy (fun (w,ac) -> (ac - noiseAutoCorr) |> abs)
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
            | Some x ->
                consensusLibrary
                |> Array.groupBy (fun x -> x.Sequence, x.GlobalMod, x.ModSequenceID, x.PepSequenceID)
                |> Array.map snd
                |> Array.filter (fun y -> x.Contains y.[0].Sequence)
            | None ->
                consensusLibrary
                |> Array.groupBy (fun x -> x.Sequence, x.GlobalMod, x.ModSequenceID, x.PepSequenceID)
                |> Array.map snd
        let queries =
            peptideList
            |> Array.map (fun entries -> entries, createPeptideQuery entries swathAnalysisParams.MatchingTolerancePPM swathAnalysisParams.QueryOffsetRange)
        let quant =
            queries
            |> Array.map (fun (entries,(rt,query)) ->
                let profiles = getRTProfiles swathIdx inReader (parseSpectrumSelection swathAnalysisParams) entries query
                profiles
                |> List.choose (fun profile ->
                // Match to look if RTProfile was found. If none was found, quantification is skipped.
                    match profile with
                    | Some rtProfile ->
                        let rtProfiles =
                            rtProfile
                            |> Array2D.toJaggedArray
                            |> Array.map (fun peak2DArr ->
                                peak2DArr
                                |> Array.map (fun peak2D ->
                                    peak2D.Rt, peak2D.Intensity, peak2D.FragmentInfo
                                )
                            )
                        let quants = 
                            rtProfiles
                            |> Array.choose (fun retInt ->
                                retInt
                                //process XIC
                                |> Array.mapi (fun i (rt,intensity,lib) ->
                                    if i = 0 || i = retInt.Length-1 || intensity > 0. then
                                        Some (rt,intensity,lib)
                                    else
                                        let rt',intensity',lib' = retInt.[i-1]
                                        if intensity' = 0. then
                                            Some (rt,intensity,lib)
                                        elif intensity' > (100. * (intensity+1.)) then
                                            None
                                        else
                                            Some (rt,intensity,lib)
                                  )
                                |> Array.choose id
                                |> Array.unzip3
                                |> fun (ret, intensity, libEntry) ->
                                    let checkSameLibEntry =
                                        libEntry
                                        |> Array.distinct
                                        |> fun x ->
                                            if x.Length > 1 then failwith "Too many Library Entries for one fragment"
                                            else
                                                x.[0]
                                    identifyPeaks ret intensity, checkSameLibEntry
                                |> fun (peaks,libEntry) ->
                                    match peaks with
                                    | [||] -> None
                                    | _    -> Some ((BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks rt), libEntry)
                            )
                            |> Array.map (fun (peak, libEntry) ->
                                BioFSharp.Mz.Quantification.HULQ.quantifyPeak peak
                                |> fun x -> x.Area, libEntry
                            )
                        Some (quants)
                    | None -> None
                )
            )
            |> List.concat
            |> Array.concat
            |> Array.map (fun (quant,x) ->
                {
                    Charge         = x.Charge
                    Iontype        = x.Iontype
                    MassOverCharge = x.MassOverCharge
                    Number         = x.Number
                    Intensity      = x.Intensity
                    PepSequenceID  = x.PepSequenceID
                    ModSequenceID  = x.ModSequenceID
                    PSMId          = x.PSMId
                    PrecursorMZ    = x.PrecursorMZ
                    ScanTime       = x.ScanTime
                    Count          = x.Count
                    Version        = x.Version
                    Sequence       = x.Sequence
                    GlobalMod      = x.GlobalMod
                    PercolatorScore= x.PercolatorScore
                    Quant          = quant
                }
            )
        tr.Dispose()
        if quant.Length >= 1 then
            FSharpAux.IO.SeqIO.Seq.CSV "\t" true false quant
            |> FSharpAux.IO.SeqIO.Seq.writeOrAppend outFilePath
        else
            ()