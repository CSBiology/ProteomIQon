namespace ProteomIQon

open System.IO
open System
open FSharpAux.Colors.Table.StatisticalGraphics24
open FSharpAux.IO
open FSharp.Stats
open FSharpAux.IO.SchemaReader
open FSharp.Plotly
open BioFSharp
open Microsoft
open Microsoft.ML
open Microsoft.ML.Data   
open Dto
open Dto.QuantificationResult
open Newtonsoft.Json
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
open BioFSharp.Mz.SearchDB

module SwathAnalysis =
    let quantifyPeptides p o dbConnection i ii = ()
    ///
    type SparsePeakArray = { 
        Length: int
        Data: System.Collections.Generic.IDictionary<int,float>
        }
   
    ///
    let dot (x:SparsePeakArray) (y:SparsePeakArray) =
        x.Data
        |> Seq.fold (fun (acc:float) xi -> 
            let present,yi = y.Data.TryGetValue xi.Key
            if present then acc + (yi * xi.Value) else acc
            ) 0.
           
    ///
    let getBinIdx' width offset x = int ((x / width) + offset)

    ///
    let peaksToNearestBinVector binWidth offset (pkarr:BioFSharp.Mz.PeakArray<_>) (minMassBoarder:float) (maxMassBoarder:float) = 
        let minMassBoarder = getBinIdx' binWidth offset minMassBoarder
        let maxMassBoarder = getBinIdx' binWidth offset maxMassBoarder
        let maxIndex = maxMassBoarder - minMassBoarder + 1        
        let keyValues = 
            pkarr 
            |> Array.choose (fun p ->  
                let index = (getBinIdx' binWidth offset p.Mz) - minMassBoarder
                if index < maxIndex-1 && index > -1 then Some (index, p.Intensity) else None
                )
            |> Array.groupBy fst 
            |> Array.map (fun (idx,data) -> idx, data |> Array.sumBy snd)
        { 
            Length = maxIndex
            Data = keyValues |> dict
        } 
    ///
    type Library =
        {
        FileName: string
        Data: ProteomIQon.Dto.PeptideIon[]
        }

    ///
    let readLibraryFrom path= 
        let tmp: ProteomIQon.Dto.PeptideIon[] = ProteomIQon.Json.ReadAndDeserialize path
        {
            FileName = Path.GetFileNameWithoutExtension path
            Data = tmp
        }


    type IonSpecies = {
        ModSequenceID                               : int
        Charge                                      : int
        GlobalMod                                   : int
        }
    
    let createIonSpecies modSequenceID charge globalMod = {
        ModSequenceID  = modSequenceID
        Charge         = charge
        GlobalMod      = globalMod
        }
    
    type QuantifiedFragment = {
        Fragment                : Dto.FragmentIon
        CorrelationWithTrace    : float
        Quant                   : float
        }
    
    let createQuantifiedFragment fragment correlationWithTrace quant = {
        Fragment                = fragment
        CorrelationWithTrace    = correlationWithTrace
        Quant                   = quant
        }

    type QuantifiedPepIon = {
        PeptideIon              : PeptideIon
        QuantifiedFragments     : QuantifiedFragment list
        FittedScanTime          : float
        }

    let createQuantifiedPepIon peptideIon quantifiedFragments fittedScanTime = {
        PeptideIon              = peptideIon
        QuantifiedFragments     = quantifiedFragments
        FittedScanTime          = fittedScanTime
        }
        


    ///
    type SwathIndexer.SwathIndexer with
            
        member this.GetRTProfilesFirstWnd(dataReader:IMzIODataReader, query: SwathQuery, getLockMz: bool, spectrumSelector: seq<SwathIndexer.MSSwath> -> seq<SwathIndexer.MSSwath> list ,?mzRangeSelector: Peak1DArray * RangeQuery -> Peak1D) =

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

                    let profile = Array2D.create query.CountMS2Masses swathSpectrum.Length (new Peak2D())

                    for specIdx = 0 to swathSpectrum.Length - 1 do

                        let swathSpec = swathSpectrum.[specIdx]
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
            )

    ///
    let getClosestMz (peaks: Peak1DArray, mzRange: RangeQuery) =
        let p1d = 
            peaks.Peaks
                .DefaultIfEmpty(new Peak1D(0., mzRange.LockValue))
                .ItemAtMin(fun x -> Math.Abs(x.Mz - mzRange.LockValue))
        if p1d.Mz > mzRange.HighValue || p1d.Mz < mzRange.LowValue then
            new Peak1D(0., mzRange.LockValue)
        else
            p1d

    ///
    let getRTProfiles (swathIndexer: SwathIndexer.SwathIndexer) (reader: IMzIODataReader) (*(spectrumSelector: seq<SwathIndexer.MSSwath> -> seq<SwathIndexer.MSSwath> list)*) (swathQuery: SwathQuery) =
        swathIndexer.GetRTProfilesFirstWnd(reader, swathQuery, false, (fun x -> [x.Take(1)]) ,getClosestMz)
    
    ///
    let getPlotFilePathFilePath outputDir (swathFileName:string) (plotName:string) (libraryFileName:string) =
        let swathFileName = swathFileName + "_" + (libraryFileName) + "_" + plotName 
        Path.Combine [|outputDir;swathFileName|]
    
    let filterePepions (processParams:Domain.ConsensusSpectralLibraryParams) (pepIons:PeptideIon[]) = 
        pepIons
        |> Array.map (fun pepIons -> {pepIons with Fragments = pepIons.Fragments |> List.filter (fun x -> x.Number > processParams.MinFragmentLadderIdx) })
        |> Array.filter (fun pepIons -> pepIons.Fragments.Length > processParams.MinFragmentCount)
        |> Array.filter (fun pepIons -> pepIons.StringSequence |> String.filter Char.IsUpper |> String.length >= processParams.MinPeptideLength)
            

    /// Retrieves the Stabw based on the fitted parameter values (HULQ output).
    let tryGetApexIntensity (fit:BioFSharp.Mz.Quantification.HULQ.QuantifiedPeak) = 
        if Array.length fit.EstimatedParams < 3 then 
            None
        else 
            Some fit.EstimatedParams.[0] 

    /// Retrieves the scan time based on the fitted parameter values (HULQ output).
    let tryGetScanTime (fit:BioFSharp.Mz.Quantification.HULQ.QuantifiedPeak) = 
        if Array.length fit.EstimatedParams < 3 then 
            None
        else 
            Some fit.EstimatedParams.[1] 

    /// Retrieves the Stabw based on the fitted parameter values (HULQ output).
    let tryGetStdev (fit:BioFSharp.Mz.Quantification.HULQ.QuantifiedPeak) = 
        if Array.length fit.EstimatedParams < 3 then 
            None
        else 
            Some fit.EstimatedParams.[2] 

    ///
    let checkApex fittedApex measuredApex = 
        fittedApex > 0.5 * measuredApex && fittedApex < 1.5 * measuredApex
    
    ///
    let checkStdev fittedstdev (peakXData:float[]) = 
        let minDist = 
            FSharp.Stats.Signal.Padding.HelperFunctions.getMinimumSpacing 
                (Array.map (fun x -> x,1.) peakXData) 
                FSharp.Stats.Signal.Padding.HelperFunctions.Float.getDiffFloat 
        let minMinusMax = 
            let min = Array.min peakXData
            let max = Array.max peakXData
            min - max
            |> abs
        fittedstdev > minDist && fittedstdev < minMinusMax
    
    let quantifyPepIon getPlotFilePathFilePath getRTProfiles (processParams:Domain.ConsensusSpectralLibraryParams) pepIon = 
        let targetFragmentMzs = 
            pepIon.Fragments
            |> List.map (fun f -> f.MeanFragMz,f.MeanRelativeIntensity_Frags)
            |> Array.ofSeq
            |> BioFSharp.Mz.PeakArray.zipMzInt
        let targetFragmentVector = 
            targetFragmentMzs 
            //|> fun x -> BioFSharp.Mz.PeakArray.peaksToNearestUnitDaltonBinVector x 100. 2000.
            |> fun (x) -> peaksToNearestBinVector processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset x (processParams.MS2ScanRange |> fst) (processParams.MS2ScanRange |> snd)
        let rtQuery   = MzIO.Processing.Query.createRangeQuery pepIon.ScanTime processParams.RtWindowWidth
        let mzQueries = 
            targetFragmentMzs 
            |> Array.map (fun pk -> MzIO.Processing.Query.createRangeQuery pk.Mz 0.05)                    
        let targetSwathQuery = MzIO.Processing.Query.createSwathQuery pepIon.PrecursorMZ rtQuery mzQueries
        let fragVecs  = 
            let tmp :Peak2D [,] option list  = getRTProfiles targetSwathQuery
            match tmp with 
            | (Some queryRes)::t ->
                let fragVecs = 
                    queryRes
                    //Out: mz inner: intensitäten zu rts 
                    |> JaggedArray.ofArray2D
                    |> JaggedArray.transpose
                    |> Array.map (fun rtGroup ->
                        rtGroup.[0].Rt, 
                        rtGroup
                        |> Array.map (fun p -> p.Mz,p.Intensity)
                        |> Array.ofSeq
                        |> BioFSharp.Mz.PeakArray.zipMzInt
                        |> fun (x) -> peaksToNearestBinVector processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset x (processParams.MS2ScanRange |> fst) (processParams.MS2ScanRange |> snd)          
                        )
                    |> Array.sortBy fst
                Some fragVecs
            | _ -> None
        let quantifiedCorrelationProfile = 
            match fragVecs with 
            | Some fragVecs -> 
                let correlationTrace = 
                    fragVecs
                    |> Array.map (fun (rt, measuredFrags) -> 
                        rt, dot targetFragmentVector measuredFrags
                        )
                let xData,yData = correlationTrace |> Array.unzip
                try
                let s = FSharpStats'.Wavelet.identify FSharpStats'.Wavelet.p xData yData
                if Array.isEmpty s then 
                    None
                elif s.Length = 1 then
                    let peakToQuantify =s.[0]
                    let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify 
                    Some (quantP,peakToQuantify,correlationTrace)
                else
                    let sorted = s |> Array.sortByDescending (fun x -> x.Apex.YVal)
                    let mainPeak,sndPeak = sorted.[0],sorted.[1]
                    let ratio = sndPeak.Apex.YVal / mainPeak.Apex.YVal
                    if ratio > processParams.MaxRatioMostAbundandVsSecondAbundandPeak then None else
                    let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy s pepIon.ScanTime
                    let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify 
                    Some (quantP,peakToQuantify,correlationTrace)
            
                with 
                | _ -> None
            | _ -> None

        match quantifiedCorrelationProfile with 
        | Some (quantifiedProfile,profile,correlationTrace) -> 
            let targetApex     = tryGetApexIntensity quantifiedProfile
            let targetScanTime = tryGetScanTime quantifiedProfile 
            let targetStdev    = tryGetStdev quantifiedProfile
            match targetApex, targetScanTime, targetStdev, fragVecs with 
            | Some apex, Some targetScanTime, Some stdev, Some fragVecs -> 
                if checkApex apex quantifiedProfile.MeasuredApexIntensity && checkStdev stdev profile.XData then
                    let profileMinRT = profile.XData |> Array.min
                    let profileMaxRT = profile.XData |> Array.max
                    let filteredfragVecs = 
                        fragVecs
                        |> Array.filter (fun (rt, mzData) -> 
                                rt >= profileMinRT && rt <= profileMaxRT
                            )
                    let fragCorrQuant =
                        pepIon.Fragments
                        |> List.map (fun frag -> 
                                let minMassBoarder = getBinIdx' processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset (processParams.MS2ScanRange |> fst)
                                let idx = (getBinIdx' processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset frag.MeanFragMz) - minMassBoarder 
                                let fragTrace = 
                                    filteredfragVecs 
                                    |> Array.map (fun (rt, mzData) -> 
                                        // can probably be checked earlier
                                        if mzData.Data.ContainsKey idx then 
                                            mzData.Data.[idx] 
                                        else 0.
                                        )
                                let corr = FSharp.Stats.Correlation.Seq.pearson profile.YData fragTrace
                                let max = fragTrace |> Array.max 
                                createQuantifiedFragment frag corr max
                            )
                        |> List.filter (fun quantFrag -> 
                            quantFrag.CorrelationWithTrace > 0.75
                            )
                    [
                        [
                        Chart.Point(correlationTrace)
                        Chart.Point([targetScanTime],[apex])
                        ]
                        |> Chart.Combine
                        Chart.Point(fragCorrQuant |> List.map (fun x -> targetScanTime,x.Quant))
                    ]
                    |> Chart.Combine
                    |> Chart.SaveHtmlAs (getPlotFilePathFilePath ((pepIon.StringSequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + pepIon.GlobalMod.ToString() + "Ch" + pepIon.Charge.ToString()) ("_") ) 
                    createQuantifiedPepIon pepIon fragCorrQuant targetScanTime
                    |> Some
                else None 

            | _ -> 
                None 
        | None -> None
    
    ///
    let labledQuantification (logger:NLog.Logger) (*(peptideLookUp:string -> int -> LookUpResult<AminoAcids.AminoAcid>)*) getPlotFilePathFilePath getRTProfiles (processParams:Domain.ConsensusSpectralLibraryParams) (pepIons:PeptideIon[]) = 
        let mutable c = 0
        pepIons
        |> Array.choose (fun pepIon -> 
                if c % 100 = 0 then logger.Trace (sprintf "peps quantified: %i" c)
                c <- c + 1
                //let unlabledPeptide = peptideLookUp pepIon.StringSequence 0
                //let labeledPeptide  = peptideLookUp pepIon.StringSequence 1
                //let targetPeptide = if pepIon.GlobalMod = 0 then unlabledPeptide else labeledPeptide 
                let quantifiedTargetIon = quantifyPepIon getPlotFilePathFilePath getRTProfiles processParams pepIon
                quantifiedTargetIon
            )


    ///
    let quantify (processParams:Domain.ConsensusSpectralLibraryParams) (outputDir:string) (targetSwathFile:string) (libraryFile:string) = 
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension libraryFile)
        logger.Trace (sprintf "Input directory containing library: %s" libraryFile)
        logger.Trace (sprintf "Input Swath file to quantify to: %s" targetSwathFile)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension targetSwathFile) + ".squant"
            Path.Combine [|outputDir;fileName|]
         
        logger.Trace "Reading and preparing .csl file for quantification"
        let libraryFile = readLibraryFrom libraryFile
        logger.Trace "Reading and preparing .csl file for quantification: finished"
        
        logger.Trace "Init connection to swath mass spectrum data."
        let inReader = Core.MzIO.Reader.getReader targetSwathFile :?> MzIO.MzSQL.MzSQL
        inReader.Connection.Open()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader
        let inTr = inReader.BeginTransaction()
        logger.Trace "Init connection to swath mass spectrum data.:finished"

        logger.Trace "Create Swath RetentionTime index"
        let retTimeIdxed = Query.getSwathIdx inReader inRunID
        let getRTProfiles = getRTProfiles retTimeIdxed inReader
        logger.Trace "Create Swath RetentionTime index:finished"

        logger.Trace "Quantify peptide ions in Swath file"
        let swathFileName = Path.GetFileNameWithoutExtension targetSwathFile
        let getPlotFilePathFilePath = getPlotFilePathFilePath outputDir swathFileName                    
        let quantifiedIons = labledQuantification logger getPlotFilePathFilePath getRTProfiles processParams libraryFile.Data
        logger.Trace "Quantify peptide ions in Swath file: finished"
        
        logger.Trace "Writing Output"
        Json.serializeAndWrite outFilePath quantifiedIons
        logger.Trace "Writing Output:finished"
        
        logger.Trace "Quantification: finished"
        
  
        
            


        
        
            
        
