namespace ProteomIQon

open System.IO
open System
open FSharpAux.Colors.Table.StatisticalGraphics24
open FSharpAux.IO
open FSharp.Stats
open FSharpAux.IO.SchemaReader
open Plotly.NET
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
open BioFSharp.Mz.Quantification
open SparsePeakArray'

module SwathAnalysis =

    /// Return Option
    let getPeakClosestPeakBy (peaks:FSharp.Stats.Signal.PeakDetection.IdentifiedPeak []) x =
        let isPartOfPeak x (p:FSharp.Stats.Signal.PeakDetection.IdentifiedPeak)  = 
            p.LeftEnd.XVal < x && p.RightEnd.XVal > x       
        match Array.tryFind (isPartOfPeak x) peaks with
        | Some p        ->  p
        | Option.None   -> Array.minBy (fun p -> abs(p.Apex.XVal - x)) peaks 

    ///// Return Option
    //let getHighestPeakBy (peaks:FSharp.Stats.Signal.PeakDetection.IdentifiedPeak []) x =
    //    let isPartOfPeak x (p:FSharp.Stats.Signal.PeakDetection.IdentifiedPeak)  = 
    //        p.LeftEnd.XVal < x && p.RightEnd.XVal > x       
    //    match Array.tryFind (isPartOfPeak x) peaks with
    //    | Some p        ->  p
    //    | Option.None   -> Array.minBy (fun p -> abs(p.Apex.XVal - x)) peaks 


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
        Quant_Max               : float
        Quant_Sum               : float
        QuantArea_Uni           : float  
        QuantArea_Trap          : float  
        QuantArea_Fit           : float  
        }
    
    let createQuantifiedFragment fragment correlationWithTrace quant_Max quant_Sum quantArea_Uni quantArea_Trap quantArea_Fit = {
        Fragment                = fragment
        CorrelationWithTrace    = correlationWithTrace
        Quant_Max               = quant_Max
        Quant_Sum               = quant_Sum
        QuantArea_Uni           = quantArea_Uni           
        QuantArea_Trap          = quantArea_Trap           
        QuantArea_Fit           = quantArea_Fit 
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
    
    let filterPepions (processParams:Domain.ConsensusSpectralLibraryParams) (pepIons:PeptideIon[]) = 
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
        (2.*fittedstdev) > minDist && fittedstdev < minMinusMax
    
    let identifyPeak pepIon (processParams:Domain.SWATHAnalysisParams) targetFragmentVector fragVecs = 
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
                Result.Error "No Peak detected in Correlation Trace"
            elif s.Length = 1 then
                let peakToQuantify =s.[0]
                let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify 
                Result.Ok (quantP,peakToQuantify,correlationTrace)
            else
                let sorted = s |> Array.sortByDescending (fun x -> x.Apex.YVal)
                let mainPeak,sndPeak = sorted.[0],sorted.[1]
                let ratio = sndPeak.Apex.YVal / mainPeak.Apex.YVal
                if ratio > processParams.MaxRatioMostAbundandVsSecondAbundandPeak then 
                    Result.Error (sprintf "sndPeak.Apex.YVal / mainPeak.Apex.YVal > %f" processParams.MaxRatioMostAbundandVsSecondAbundandPeak) 
                else
                    let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy s pepIon.ScanTime
                    let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify 
                    Result.Ok (quantP,peakToQuantify,correlationTrace)
            with 
            | _ -> Result.Error "PrecursorMz not covered by Swath windows."
        | _ -> Result.Error "PrecursorMz not covered by Swath windows."
    
    ///
    let quantifyPeakBy estApex estMeanX estStabw (x: float []) (y: float []) =
        //match estimateMoments p with 
        //| Some moments -> 
            let gaussian = BioFSharp.Mz.Quantification.HULQ.tryFitGaussian estApex estMeanX estStabw x y
            //let emg      = BioFSharp.Mz.Quantification.HULQ.tryFitEMG moments.ModeY moments.MeanX moments.Std moments.Tau p.XData p.YData
            match gaussian(*, emg*) with 
            | Some g(*, Some emg*) -> 
                let finalFit   = g // selectModel [|g;emg|]
                BioFSharp.Mz.Quantification.HULQ.createQuantifiedPeak  (Some finalFit.Model) finalFit.YPredicted finalFit.EstimatedParams finalFit.StandardErrorOfPrediction finalFit.Area estApex
            //| Some finalFit, Option.None  | Option.None, Some finalFit -> 
            //    createQuantifiedPeak  (Some finalFit.Model) finalFit.YPredicted finalFit.EstimatedParams finalFit.StandardErrorOfPrediction finalFit.Area  p.Apex.YVal
            | _                       ->
                let area = BioFSharp.Mz.Quantification.Integration.trapezEstAreaOf x y
                BioFSharp.Mz.Quantification.HULQ.createQuantifiedPeak Option.None [||] [||] nan area estApex
        //| Option.None   -> 
        //    let area = trapezEstAreaOf p.XData p.YData
        //    createQuantifiedPeak Option.None [||] [||] nan area p.Apex.YVal  

    let quantifyFragments diagCharts (logger:NLog.Logger) getPlotFilePathFilePath (processParams:Domain.SWATHAnalysisParams) pepIon (fragVecs:(float*SparsePeakArray)[] option) (quantifiedCorrelationProfile:(Mz.Quantification.HULQ.QuantifiedPeak*Signal.PeakDetection.IdentifiedPeak*((float*float) []))) = 
        let (quantifiedProfile,profile,correlationTrace) = quantifiedCorrelationProfile 
        let targetApex     = tryGetApexIntensity quantifiedProfile
        let targetScanTime = tryGetScanTime quantifiedProfile 
        let targetStdev    = tryGetStdev quantifiedProfile
        match targetApex, targetScanTime, targetStdev, fragVecs with 
        | Some apex, Some targetScanTime, Some stdev, Some fragVecs -> 
            if checkApex apex quantifiedProfile.MeasuredApexIntensity |> not then 
                logger.Trace (sprintf "checkApex failed with apex: %f quantifiedProfile.MeasuredApexIntensity: %f" apex quantifiedProfile.MeasuredApexIntensity)
                Result.Error ("checkApex failed.")
            elif checkStdev stdev profile.XData |> not then
                logger.Trace (sprintf "checkStdev failed with stdev: %f xData: %A" stdev profile.XData)
                Result.Error ("checkStdev failed.")
            else
                try
                let xCorrelationTrace,yCorrelationTrace = correlationTrace |> Array.unzip
                let profileMinRT = 
                    match xCorrelationTrace |> Array.tryFindBack (fun rt -> rt <= (targetScanTime - (2.*stdev))) with
                    | Some x -> x
                    | None   -> xCorrelationTrace |> Array.min
                let profileMaxRT = 
                    match xCorrelationTrace |> Array.tryFind (fun rt -> rt >= (targetScanTime + (2.*stdev))) with 
                    | Some x -> x
                    | None   -> xCorrelationTrace |> Array.max
                let correlationTraceToCorrelate = 
                    correlationTrace 
                    |> Array.filter (fun (rt,y) -> rt >= profileMinRT && rt <= profileMaxRT)
                    |> Array.map snd
                let filteredfragVecs = 
                    fragVecs
                    |> Array.filter (fun (rt, mzData) -> 
                            rt >= profileMinRT && rt <= profileMaxRT
                        )
                let fragCorrQuant,xFrags,yFrags =
                    try
                    pepIon.Fragments
                    |> List.map (fun frag -> 
                            let idx = initMzToBinIdx processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset frag.MeanFragMz
                            let fragTraceX,fragTraceY = 
                                filteredfragVecs 
                                |> Array.map (fun (rt, mzData) -> 
                                    // can probably be checked earlier
                                    if mzData.Data.ContainsKey idx then 
                                        rt,mzData.Data.[idx] 
                                    else rt,0.
                                    )
                                |> Array.unzip
                            let corr = FSharp.Stats.Correlation.Seq.pearson correlationTraceToCorrelate fragTraceY
                            let max = if fragTraceY |> Array.isEmpty then 0. else fragTraceY |> Array.max 
                            let sum = fragTraceY |> Array.sum       
                            let trapAreaUni = BioFSharp.Mz.Quantification.Integration.trapezEstAreaOfUniform fragTraceX fragTraceY
                            let trapArea = BioFSharp.Mz.Quantification.Integration.trapezEstAreaOf fragTraceX fragTraceY
                            let qp = quantifyPeakBy max targetScanTime stdev fragTraceX fragTraceY
                            createQuantifiedFragment frag corr max sum trapAreaUni trapArea qp.Area,fragTraceX,fragTraceY
                        )
                    |> List.filter (fun (quantFrag,_,_) -> 
                        quantFrag.CorrelationWithTrace > 0.75
                        )
                    |> List.unzip3
                    with 
                    | _ -> [],[],[]
                if diagCharts then 
                    [
                        [
                        Chart.Point(correlationTrace)
                        Chart.Point([targetScanTime],[apex])
                        ]
                        |> Chart.Combine
                        Chart.Point(fragCorrQuant |> List.map (fun x -> targetScanTime,x.Quant_Max))
                        |> Chart.withTraceName "Max"
                        Chart.Point(fragCorrQuant |> List.map (fun x -> targetScanTime,x.Quant_Sum))
                        |> Chart.withTraceName "Sum"
                        List.map2 (fun (x:float[]) (y:float[]) -> Chart.Line(x,y)) xFrags yFrags
                        |> Chart.Combine 
                    ]
                    |> Chart.Combine
                    |> Chart.SaveHtmlAs (getPlotFilePathFilePath ((pepIon.StringSequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + pepIon.GlobalMod.ToString() + "Ch" + pepIon.Charge.ToString()) ("_") ) 
                
                match fragCorrQuant with 
                | [] -> 
                    Result.Error "No Fragments with Correlation > 0.75"
                | _  ->
                    createQuantifiedPepIon pepIon fragCorrQuant targetScanTime
                    |> Result.Ok
                with 
                |_ -> Result.Error ("No Fragments.")
        | _ -> Result.Error "Fit did not return all Parameters"

    let quantifyPepIon diagCharts logger getPlotFilePathFilePath getRTProfiles (processParams:Domain.SWATHAnalysisParams) pepIon = 
        let targetFragmentMzs = 
            pepIon.Fragments
            |> List.map (fun f -> f.MeanFragMz,f.MeanRelativeIntensity_Frags)
            |> Array.ofSeq
            |> BioFSharp.Mz.PeakArray.zipMzInt
        let targetFragmentVector = 
            targetFragmentMzs 
            |> peaksToNearestBinVector processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset (processParams.MS2ScanRange |> fst) (processParams.MS2ScanRange |> snd)
        let rtQuery   = MzIO.Processing.Query.createRangeQuery pepIon.ScanTime processParams.RtWindowWidth
        let mzQueries = 
            targetFragmentMzs 
            |> Array.map (fun pk -> MzIO.Processing.Query.createRangeQuery pk.Mz processParams.FragMatchingBinWidth)                    
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
                        |> peaksToNearestBinVector processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset (processParams.MS2ScanRange |> fst) (processParams.MS2ScanRange |> snd)          
                        )
                    |> Array.sortBy fst
                Some fragVecs
            | _ -> None

        //let quantifiedCorrelationProfile = identifyPeak pepIon processParams targetFragmentVector fragVecs
        identifyPeak pepIon processParams targetFragmentVector fragVecs
        |> Result.bind (quantifyFragments diagCharts logger getPlotFilePathFilePath processParams pepIon fragVecs)           
     
    ///
    let labledQuantification diagCharts (logger:NLog.Logger) (*(peptideLookUp:string -> int -> LookUpResult<AminoAcids.AminoAcid>)*) getPlotFilePathFilePath getRTProfiles (processParams:Domain.SWATHAnalysisParams) (pepIons:PeptideIon[]) = 
        let results = 
            pepIons
            |> Array.mapi (fun i pepIon -> 
                    if i % 100 = 0 then logger.Trace (sprintf "peps quantified: %i" i)
                    //let unlabledPeptide = peptideLookUp pepIon.StringSequence 0
                    //let labeledPeptide  = peptideLookUp pepIon.StringSequence 1
                    //let targetPeptide = if pepIon.GlobalMod = 0 then unlabledPeptide else labeledPeptide 
                    let quantifiedTargetIon = quantifyPepIon diagCharts logger getPlotFilePathFilePath getRTProfiles processParams pepIon
                    quantifiedTargetIon
                )
        results
        |> Array.groupBy (fun x ->
            match x with 
            | Result.Ok x    -> "Result"
            | Result.Error s -> s
        )
        |> Array.map (fun (message,items) -> 
                message,items.Length
            )
        |> Chart.Column
        |> Chart.SaveHtmlAs (getPlotFilePathFilePath ("Errors") ("_") ) 
        results
        |> Array.choose (fun x ->
                match x with 
                | Result.Ok x -> Some x
                | _           -> None
            )


    ///
    let quantify diagCharts (processParams:Domain.SWATHAnalysisParams) (outputDir:string) (targetSwathFile:string) (libraryFile:string) = 
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
        let quantifiedIons = labledQuantification diagCharts logger getPlotFilePathFilePath getRTProfiles processParams libraryFile.Data
        logger.Trace "Quantify peptide ions in Swath file: finished"
        
        logger.Trace "Writing Output"
        Json.serializeAndWrite outFilePath quantifiedIons
        logger.Trace "Writing Output:finished"
        
        logger.Trace "Quantification: finished"
        
  
        
            


        
        
            
        
