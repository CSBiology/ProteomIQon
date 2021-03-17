namespace ProteomIQon

open System.IO
open System.Data.SQLite
open ProteomIQon.Core
open Core.MzIO
open Dto
open FSharp.Stats
open BioFSharp.Mz.Quantification
open BioFSharp.Mz
open FSharpAux.IO.SchemaReader
open Plotly.NET
open BioFSharp
open MzIO.Processing
open BioFSharp.Mz.SearchDB
open SearchDB'

module PSMBasedQuantification =

    type PeptideIon = 
        {
            Sequence             : string
            GlobalMod            : int
            Charge               : int
            ModSequenceID        : int
            PepSequenceID        : int
        }
            
    type averagePSM = {
        MeanPrecMz   : float
        MeanScanTime : float
        WeightedAvgScanTime:float
        MeanScore   : float
        X_Xic         : float []
        Y_Xic         : float []
        Y_Xic_uncorrected: float []
        }

    ///
    let createAveragePSM meanPrecMz meanScanTime weightedAvgScanTime meanScore xXic yXic yXic_uncorrected = {
        MeanPrecMz    = meanPrecMz
        MeanScanTime  = meanScanTime
        WeightedAvgScanTime= weightedAvgScanTime
        MeanScore = meanScore
        X_Xic         = xXic
        Y_Xic         = yXic
        Y_Xic_uncorrected= yXic_uncorrected
        }


    type PeakComparison = {
        Mz: float
        MeasuredIntensity: float
        MeasuredIntensityCorrected:float
        PredictedRelFrequency:float
        }


    type ClusterComparison = {
        PeakComparisons     : PeakComparison []
        KLDiv_UnCorrected   : float
        KLDiv_Corrected     : float
        }

    ///
    let getBaseLineCorrectionOffsetAt tarRT x_Xic y_Xic y_Xic_uncorrected =
        let (rt,y,yUncorr) = 
            Array.zip3 x_Xic y_Xic y_Xic_uncorrected 
            |> Array.minBy (fun (rt,y,yUncorr) -> abs (rt - tarRT))
        yUncorr - y
       
    ///
    let getClosestMs1 (ms1s: MzIO.Model.MassSpectrum []) scanTime = 
         ms1s
         |> Array.minBy (fun ms -> abs (MassSpectrum.getScanTime ms - scanTime))

    ///
    let getSpec (reader:MzIO.IO.IMzIODataReader) (ms1: MzIO.Model.MassSpectrum)  =
        Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(ms1.ID).Peaks)
        |> fun (mzData,intensityData) -> PeakArray.zip mzData intensityData
       
    /// Calculates the Kullback-Leibler divergence Dkl(p||q) from q (theory, model, description, or approximation of p) 
    /// to p (the "true" distribution of data, observations, or a precisely calculated theoretical distribution).
    let klDiv (p:float []) (q:float []) = 
        Array.fold2 (fun acc p q -> (System.Math.Log(p/q)*p) + acc ) 0. p q
     


    ///
    let substractBaseLine (logger: NLog.Logger) (baseLineParams:Domain.BaseLineCorrection) (yData:float []) =
        if yData.Length > 500 then
            yData
        else
            let baseLine = FSharp.Stats.Signal.Baseline.baselineAls' baseLineParams.MaxIterations baseLineParams.Lambda baseLineParams.P yData |> Array.ofSeq
            Array.map2 (fun y b ->
                           let c = y - b
                           if c < 0. then 0. else c
                       ) yData baseLine

    ///
    let initGetProcessedXIC logger (baseLineCorrection:Domain.BaseLineCorrection option) reader idx scanTimeWindow mzWindow_Da meanScanTime meanPrecMz =
        let rtQuery = Query.createRangeQuery meanScanTime scanTimeWindow
        let mzQuery = Query.createRangeQuery meanPrecMz mzWindow_Da
        let retData',itzData' =
            let tmp =
                Query.getXIC reader idx rtQuery mzQuery
                |> Array.map (fun p -> p.Rt , p.Intensity)
            tmp
            |> Array.mapi (fun i (rt,intensity) ->
                            if i = 0 || i = tmp.Length-1 || intensity > 0. then
                                Some (rt,intensity)
                            else
                                let rt',intensity' = tmp.[i-1]
                                if intensity' = 0. then
                                    Some (rt,intensity)
                                elif intensity' > (100. * (intensity+1.)) then
                                    None
                                else
                                    Some (rt,intensity)
                          )
            |> Array.choose id
            |> Array.unzip
        match baseLineCorrection with
        | Some baseLineParams ->
            retData', substractBaseLine logger baseLineParams itzData', itzData'
        | None ->
            retData',itzData',itzData'

    ///
    let initGetIsotopicEnvelope reader idx scanTimeWindow mzWindow_Da ch meanScanTime meanPrecMz =
        let rtQuery   = Query.createRangeQuery meanScanTime scanTimeWindow
        let mzQueries = 
            [|
                for i = 0 to 2 do 
                    let mz = meanPrecMz + (float i) * (Mass.Table.PMassInU / (ch|> float)) 
                    Query.createRangeQuery mz mzWindow_Da 
            |]
            |> Array.filter (fun x ->  abs (meanPrecMz-x.LockValue) < 1.)                 
        let retData',itzData' =
            let query = Query.getXICs reader idx rtQuery mzQueries 
            [|
                for i = 0 to query.[*,0].Length-1 do 
                let tmp = query.[i,*]
                yield (tmp.[0].Rt,tmp |> Array.sumBy (fun p -> p.Intensity))

            |]
            |> Array.unzip 
        retData',itzData'

    ///
    let weightedMean (weights:seq<'T>) (items:seq<'T>) =
        let sum,n = Seq.fold2 (fun (sum,n) w i -> w*i+sum,n + w ) (0.,0.) weights items
        sum / n
    ///
    let average getXic (psms:(PSMStatisticsResult*float) []) =
            let meanPrecMz   = psms |> Seq.meanBy (fun (psm,m) -> psm.PrecursorMZ)
            let meanScanTime = psms |> Seq.meanBy (fun (psm,m) -> psm.ScanTime)
            let meanScore = psms |> Seq.averageBy (fun (psm,m) -> psm.PercolatorScore)
            let psms' = 
                let tmp = Array.sortByDescending (fun (psm,m) -> m) psms
                if tmp.Length > 3 then tmp.[..2] else tmp 
            let weightedAvgScanTime =
                let scanTimes = psms' |> Array.map (fun (psm,m) -> psm.ScanTime)
                let weights = psms' |> Array.map snd
                weightedMean weights scanTimes
            let (retData,itzDataCorrected,ItzDataUncorrected) = getXic weightedAvgScanTime meanPrecMz
            createAveragePSM meanPrecMz meanScanTime weightedAvgScanTime meanScore retData itzDataCorrected ItzDataUncorrected


    type InferredPeak = {
        Model                       :HULQ.PeakModel option
        Area                        :float
        StandardErrorOfPrediction   :float
        MeasuredApexIntensity       :float
        EstimatedParams             :float[]
        X_Xic                       :float[]
        Y_Xic                       :float[]
        Y_Xic_uncorrected           :float[]
        xPeak                       :float[]
        yFitted                     :float[]
        }

    ///
    let quantifyInferredPeak getXic identifyPeaks targetMz targetScanTime inferredScanTime =
        try
            let (retData,itzData,uncorrectedItzData)   =
                getXic targetScanTime targetMz
            let peaks = identifyPeaks retData itzData
            if Array.isEmpty peaks then
                {
                    Model                       = None 
                    Area                        = nan
                    StandardErrorOfPrediction   = nan
                    MeasuredApexIntensity       = nan
                    EstimatedParams             = [||]
                    X_Xic                       = [||]
                    Y_Xic                       = [||]
                    Y_Xic_uncorrected           = [||]
                    xPeak                       = [||]
                    yFitted                     = [||]
                }
            else
                let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks inferredScanTime
                let quantP         = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify
                {
                    Model                     = quantP.Model 
                    Area                      = quantP.Area
                    StandardErrorOfPrediction = quantP.StandardErrorOfPrediction
                    MeasuredApexIntensity     = quantP.MeasuredApexIntensity
                    EstimatedParams           = quantP.EstimatedParams
                    X_Xic                     = retData
                    Y_Xic                     = itzData
                    Y_Xic_uncorrected         = uncorrectedItzData
                    xPeak                     = peakToQuantify.XData
                    yFitted                   = quantP.YPredicted
                }
        with
        | _ -> 
            {
                Model                       = None
                Area                        = nan
                StandardErrorOfPrediction   = nan
                MeasuredApexIntensity       = nan
                EstimatedParams             = [||]
                X_Xic                       = [||]
                Y_Xic                       = [||]
                Y_Xic_uncorrected           = [||]
                xPeak                       = [||]
                yFitted                     = [||]
            }

    /// Calculates the difference between the scan time used to retreave the peak and the fitted peak midpoint.
    let searchRTMinusFittedRtTarget searchRT (fit:HULQ.QuantifiedPeak) = 
        try
            searchRT - fit.EstimatedParams.[1] 
        with
        | _ -> nan

    /// Calculates the difference between the scan time used to retreave the peak and the fitted peak midpoint.
    let searchRTMinusFittedRtInferred searchRT fit = 
        try
            searchRT - fit.EstimatedParams.[1] 
        with
        | _ -> nan
    
    /// Aims to provide the best scan time estimate when quanitfying a inferred peak.
    let chooseScanTime maxDiff searchRTMinusFittedRT initialScanTime (quantP:HULQ.QuantifiedPeak) = 
        if Array.isEmpty quantP.EstimatedParams then
            initialScanTime 
        elif abs searchRTMinusFittedRT > maxDiff then
            initialScanTime 
        else
            quantP.EstimatedParams.[1]

    ///
    let saveChart sequence globalMod ch (xXic:float[]) (yXic:float[]) ms2s avgScanTime (xToQuantify:float[]) (ypToQuantify:float[]) (fitY:float[])
            (xXicInferred:float[]) (yXicinferred:float[]) (xInferred:float[]) (inferredFit:float[]) (*(xEnvelopeSum:float[]) (yEnvelopeSum:float[])*) (peaks:FSharp.Stats.Signal.PeakDetection.IdentifiedPeak []) (pattern:PeakComparison []) plotDirectory =
        let xic = 
            [
            Chart.Point(xXic, yXic)                     |> Chart.withTraceName "Target XIC"
            peaks
            |> Array.map (fun x -> Chart.Point(x.XData,x.YData)) 
            |> Chart.Combine
            Chart.Point(ms2s)                           |> Chart.withTraceName "MS2s with scores"
            Chart.Point([avgScanTime],[1.])             |> Chart.withTraceName "Weighted Mean of Ms2 scan times"
            Chart.Point((xToQuantify), (ypToQuantify))  |> Chart.withTraceName "Identified Target Peak"
            Chart.Line(xToQuantify,fitY)                |> Chart.withTraceName "Fit of target Peak"
            Chart.Point(xXicInferred, yXicinferred)     |> Chart.withTraceName "Inferred XIC"
            Chart.Line(xInferred,inferredFit)           |> Chart.withTraceName "Fit of inferred Peak"
            //Chart.Point(xEnvelopeSum, yEnvelopeSum)     |> Chart.withTraceName "Target Envelope Sum"
            ]
            |> Chart.Combine
        let pattern = 
            [
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensity))          |> Chart.withTraceName "Measured"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensityCorrected)) |> Chart.withTraceName "Measured Corrected"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.PredictedRelFrequency))      |> Chart.withTraceName "Predicted Relative Frequency"
            //Chart.Point(xEnvelopeSum, yEnvelopeSum)     |> Chart.withTraceName "Target Envelope Sum"
            ]
            |> Chart.Combine
        [xic;pattern]
        |> Chart.Stack(2, 0.1)
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i" sequence globalMod)
        |> Chart.withSize(2500.,800.)
        |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + globalMod.ToString() + "Ch" + ch.ToString())|])
            
           
    // Method is based on: https://doi.org/10.1021/ac0600196
    /// Estimates the autocorrelation at lag 1 of a blank signal (containing only noise). Subsequently, the signal of interest is smoothed
    /// several times by a savitzky golay filter using constant polynomial order and variing windowWidth. For each iteration, the deviation
    /// of the smoothed to the original signal is computed and the autocorrelation at lag 1 of this residual noise is computed. The function returns the optimized
    /// window width yielding a autocorrelation at lag 1 closest to the value computed for the blank signal.
    let optimizeWindowWidth polOrder (windowWidthToTest:int[]) noiseAutoCorr (signalOfInterest:float[]) =
        let signalOfInterest' = signalOfInterest |> vector
        //let noiseAutoCorr = Correlation.Vector.autoCorrelation 1 (blankSignal |> vector)
        let filterF w yData = FSharp.Stats.Signal.Filtering.savitzkyGolay w polOrder 0 0 yData |> vector
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
        
    ///
    let calcCorrelation (xValues:float []) (quantifiedPeak:HULQ.QuantifiedPeak) (inferredPeak:InferredPeak) = 
        let getValue (model:HULQ.PeakModel) estParams =
            match model with
            | HULQ.PeakModel.Gaussian m -> 
                m.GetFunctionValue (vector estParams)
            | HULQ.PeakModel.EMG m -> 
                m.GetFunctionValue (vector estParams)
        match quantifiedPeak.Model, inferredPeak.Model with 
        | Some q , Some i -> 
            let xValuesBW = 
                [|for i = 1 to xValues.Length-1 do abs(xValues.[i] - xValues.[i-1])|]
                |> Array.median
            let xValues = [|xValues.[0] .. xValuesBW/2. .. xValues.[xValues.Length-1]|]
            let fQ = getValue q quantifiedPeak.EstimatedParams
            let yQ = xValues |> Array.map fQ
            let fI = getValue i inferredPeak.EstimatedParams
            let yI = xValues |> Array.map fI
            FSharp.Stats.Correlation.Seq.pearson yQ yI
        | _ -> nan

    ///Predicts an isotopic distribution of the given formula at the given charge, normalized by the sum of probabilities, using the MIDAs algorithm
    let generateIsotopicDistributionOfFormulaBySum (charge:int) (seq:AminoAcids.AminoAcid list) =
        seq
        |> BioList.toFormula
        |> Formula.add Formula.Table.H2O
        |> IsotopicDistribution.MIDA.ofFormula IsotopicDistribution.MIDA.normalizeByProbSum 0.01 0.001 charge 
        |> Array.ofList

    ///
    let initComparePredictedAndMeasuredIsotopicCluster inReader ms1s ms1AccuracyEstimate x_Xic y_Xic y_Xic_uncorrected ch peptideSequence tarRt tarMz =    
        /// IsotopicCluster
        let targetIsotopicPattern_predicted = 
            generateIsotopicDistributionOfFormulaBySum ch peptideSequence
    
        let baseLineCorrectionF = getBaseLineCorrectionOffsetAt tarRt x_Xic y_Xic y_Xic_uncorrected
        let closestMS1 = getClosestMs1 ms1s tarRt
        let peaks' = 
            getSpec inReader closestMS1
            |> Array.filter (fun x -> x.Mz < tarMz + 1. && x.Mz > tarMz - 0.6)

        let recordedVsPredictedPattern = 
            targetIsotopicPattern_predicted
            |> Array.choose (fun (mz,relFreq) -> 
                match peaks' |> Array.tryFind (fun peak -> abs(peak.Mz - mz) < 4. * ms1AccuracyEstimate  ) with 
                | None -> None
                | Some peak -> Some (peak,relFreq)
                )
            |> Array.groupBy fst
            |> Array.map (fun ((peak),list) -> 
                peak.Mz,peak.Intensity,list |> Array.sumBy snd  (*,list |>*)
                )
            |> Array.map (fun (mz,measuredIntensity,predictedRelFrequency) -> 
                    {Mz=mz;MeasuredIntensity= measuredIntensity;MeasuredIntensityCorrected= measuredIntensity - baseLineCorrectionF;PredictedRelFrequency= predictedRelFrequency}
                )
            |> Array.filter (fun (isoP:PeakComparison) -> isoP.MeasuredIntensityCorrected > 0.)
        let recordedVsPredictedPatternNorm = 
            let sumMeasured = recordedVsPredictedPattern |> Array.sumBy (fun x -> x.MeasuredIntensity)
            let sumMeasuredCorr = recordedVsPredictedPattern |> Array.sumBy (fun x -> x.MeasuredIntensityCorrected)
            let sumMeasuredPred = recordedVsPredictedPattern |> Array.sumBy (fun x -> x.PredictedRelFrequency)
            recordedVsPredictedPattern
            |>Array.map (fun isoP -> 
                {isoP with 
                    MeasuredIntensity = isoP.MeasuredIntensity / sumMeasured; 
                    MeasuredIntensityCorrected = isoP.MeasuredIntensityCorrected / sumMeasuredCorr; 
                    PredictedRelFrequency = isoP.PredictedRelFrequency / sumMeasuredPred }
                )
        let klUnCorr = klDiv (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.MeasuredIntensity)) (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.PredictedRelFrequency))
        let klCorr   = klDiv (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.MeasuredIntensityCorrected)) (recordedVsPredictedPatternNorm |> Array.map (fun x -> x.PredictedRelFrequency))
        {
            PeakComparisons     = recordedVsPredictedPatternNorm
            KLDiv_UnCorrected   = klUnCorr
            KLDiv_Corrected     = klCorr
        }
       
    ///
    let quantifyPeptides (processParams:Domain.QuantificationParams) (outputDir:string) (cn:SQLiteConnection) (instrumentOutput:string) (scoredPSMs:string)  =

        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension scoredPSMs)
        logger.Trace (sprintf "Input file: %s" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        logger.Trace (sprintf "Now performing Quantification using: %s and %s, Results will be written to: %s" instrumentOutput scoredPSMs outputDir)
        // initialize Reader and Transaction
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".quant"
            Path.Combine [|outputDir;fileName|]
        logger.Trace (sprintf "outFilePath:%s" outFilePath)
        //
        let plotDirectory =
            let fileName = sprintf "%s_plots" (Path.GetFileNameWithoutExtension instrumentOutput)
            let path = Path.Combine [|outputDir;fileName|]
            if System.IO.Directory.Exists path then
                path
            else
                System.IO.Directory.CreateDirectory path |> ignore
                path
        logger.Trace (sprintf "plotDirectory:%s" plotDirectory)
        logger.Trace "Copy peptide DB into Memory"
        let memoryDB = SearchDB.copyDBIntoMemory cn
        logger.Trace "Copy peptide DB into Memory: finished"
        logger.Trace "Get peptide lookUp function"
        let dBParams     = getSDBParams memoryDB
        //let massLookUp = prepareSelectMassByModSequenceAndGlobalMod memoryDB
        let peptideLookUp = getThreadSafePeptideLookUpFromFileBySequenceAndGMod memoryDB dBParams
        let calcIonSeries aal = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal
        logger.Trace "Get peptide lookUp function: finished"
        // initialize Reader and Transaction
        logger.Trace "Init connection to mass spectrum data."
        let inReader = Core.MzIO.Reader.getReader instrumentOutput :?> MzIO.MzSQL.MzSQL
        inReader.Connection.Open()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader       
        let inTr = inReader.BeginTransaction()
        let countMatchedMasses (peptide: LookUpResult<AminoAcids.AminoAcid>)(psms: PSMStatisticsResult []) =
            let frag = 
                let ionSeries = (calcIonSeries peptide.BioSequence).TargetMasses
                [1. .. 2.]
                |> List.map (fun ch -> 
                    ionSeries 
                    |> List.map (fun x -> x.MainPeak.Mass)
                    |> List.map (fun x -> Mass.toMZ x ch)
                )
                |> List.concat                    
            psms
            |> Array.map (fun psm -> 
                let spec = inReader.ReadSpectrumPeaks psm.PSMId
                let sum = 
                    spec.Peaks 
                    |> Seq.filter (fun peak -> 
                        frag
                        |> List.exists (fun ion -> abs (ion - peak.Mz) <= (Mass.deltaMassByPpm 100. peak.Mz))
                    )
                    |> Seq.sumBy (fun x -> x.Intensity)
                psm,sum
            )
        logger.Trace "Create RetentionTime index"
        let retTimeIdxed = Query.getMS1RTIdx inReader inRunID
        logger.Trace "Create RetentionTime index:finished"
        
        logger.Trace "Read and sort ms1s"
        /// Returns a Sequence containing all MassSpectra of a single MS-run
        let massSpectra = inReader.ReadMassSpectra(inRunID)
        /// Returns a Array that contains all MS1s sorted by their scanTime
        let ms1SortedByScanTime =
            massSpectra
            |> Seq.filter (fun ms -> MassSpectrum.getMsLevel ms = 1)
            |> Seq.sortBy MassSpectrum.getScanTime
            |> Array.ofSeq
        logger.Trace "Read and sort ms1s:finished"
        
        logger.Trace "Read scored PSMs"
        ///
        let peptides =
            Csv.CsvReader<PSMStatisticsResult>(SchemaMode=Csv.Fill).ReadFile(scoredPSMs,'\t',false,1)
            |> Array.ofSeq
        logger.Trace "Read scored PSMs:finished"
 
        logger.Trace "Estimate ms1 mz accuracy"
        ///
        let ms1AccuracyEstimate = 
            peptides
            |> Seq.stDevBy (fun x -> abs(x.PrecursorMZ - Mass.toMZ x.TheoMass (float x.Charge)) )
        logger.Trace (sprintf "Estimate ms1 mz accuracy:finished, Accuracy: %f" ms1AccuracyEstimate) 
            
        logger.Trace "init lookup functions"        
        ///
        let comparePredictedAndMeasuredIsotopicCluster = 
            initComparePredictedAndMeasuredIsotopicCluster inReader ms1SortedByScanTime ms1AccuracyEstimate     
        
        ///
        let getXIC = 
            initGetProcessedXIC logger processParams.BaseLineCorrection inReader retTimeIdxed processParams.XicExtraction.ScanTimeWindow processParams.XicExtraction.MzWindow_Da     
        
        ///
        let identifyPeaks = 
            initIdentifyPeaks processParams.XicExtraction.XicProcessing
        
        logger.Trace "init lookup functions:finished"
        
        logger.Trace "init quantification functions"
        ///
        let labledQuantification (pepIon:PeptideIon) (psms:PSMStatisticsResult []) = 
            try
            let bestQValue,bestPepValue,prots = psms |> Array.minBy (fun x -> x.QValue) |> fun x -> x.QValue, x.PEPValue,x.ProteinNames
            let unlabledPeptide = peptideLookUp pepIon.Sequence 0
            let labeledPeptide  = peptideLookUp pepIon.Sequence 1
            let targetPeptide = if pepIon.GlobalMod = 0 then unlabledPeptide else labeledPeptide            
            let psmsWithMatchedSums = countMatchedMasses targetPeptide psms 
            let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
            let averagePSM = average getXIC psmsWithMatchedSums
            let avgMass = Mass.ofMZ (averagePSM.MeanPrecMz) (pepIon.Charge |> float)
            let peaks = identifyPeaks averagePSM.X_Xic averagePSM.Y_Xic
            if Array.isEmpty peaks then 
                Chart.Point(averagePSM.X_Xic, averagePSM.Y_Xic)
                |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i_noPeaks" pepIon.Sequence pepIon.GlobalMod)
                |> Chart.withSize(1500.,800.)
                |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((pepIon.Sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + pepIon.GlobalMod.ToString() + "Ch" + pepIon.Charge.ToString() + "_notQuantified")|])
                None
            else
            let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks averagePSM.WeightedAvgScanTime
            let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify                        
            let searchRTMinusFittedRT = searchRTMinusFittedRtTarget averagePSM.WeightedAvgScanTime quantP
            let inferredScanTime = chooseScanTime processParams.XicExtraction.ScanTimeWindow searchRTMinusFittedRT averagePSM.WeightedAvgScanTime quantP 
            let clusterComparison_Target = comparePredictedAndMeasuredIsotopicCluster averagePSM.X_Xic averagePSM.Y_Xic averagePSM.Y_Xic_uncorrected pepIon.Charge targetPeptide.BioSequence quantP.EstimatedParams.[1] averagePSM.MeanPrecMz            
            if pepIon.GlobalMod = 0 then
                let mz_Heavy = Mass.toMZ (labeledPeptide.Mass) (pepIon.Charge|> float)
                let inferred_Heavy = quantifyInferredPeak getXIC identifyPeaks mz_Heavy averagePSM.WeightedAvgScanTime inferredScanTime
                let searchRTMinusFittedRT_Heavy = searchRTMinusFittedRtInferred inferredScanTime inferred_Heavy
                let clusterComparison_Heavy = comparePredictedAndMeasuredIsotopicCluster inferred_Heavy.X_Xic inferred_Heavy.Y_Xic inferred_Heavy.Y_Xic_uncorrected pepIon.Charge labeledPeptide.BioSequence quantP.EstimatedParams.[1] mz_Heavy
                let corrLightHeavy  = calcCorrelation averagePSM.X_Xic quantP inferred_Heavy  
                //let chart = saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                //                    peakToQuantify.XData peakToQuantify.YData quantP.YPredicted inferred_Heavy.X_Xic inferred_Heavy.Y_Xic inferred_Heavy.xPeak inferred_Heavy.yFitted peaks clusterComparison_Target.PeakComparisons plotDirectory
                {
                StringSequence                              = pepIon.Sequence
                GlobalMod                                   = pepIon.GlobalMod
                Charge                                      = pepIon.Charge
                PepSequenceID                               = pepIon.PepSequenceID
                ModSequenceID                               = pepIon.ModSequenceID
                PrecursorMZ                                 = averagePSM.MeanPrecMz
                MeasuredMass                                = avgMass
                TheoMass                                    = unlabledPeptide.Mass
                AbsDeltaMass                                = abs(avgMass-unlabledPeptide.Mass)
                MeanPercolatorScore                         = averagePSM.MeanScore
                QValue                                      = bestQValue
                PEPValue                                    = bestPepValue
                ProteinNames                                = prots
                QuantMz_Light                               = averagePSM.MeanPrecMz
                Quant_Light                                 = quantP.Area
                MeasuredApex_Light                          = quantP.MeasuredApexIntensity
                Seo_Light                                   = quantP.StandardErrorOfPrediction
                Params_Light                                = quantP.EstimatedParams            
                Difference_SearchRT_FittedRT_Light          = searchRTMinusFittedRT
                KLDiv_Observed_Theoretical_Light            = clusterComparison_Target.KLDiv_UnCorrected
                KLDiv_CorrectedObserved_Theoretical_Light   = clusterComparison_Target.KLDiv_Corrected
                QuantMz_Heavy                               = mz_Heavy
                Quant_Heavy                                 = inferred_Heavy.Area
                MeasuredApex_Heavy                          = inferred_Heavy.MeasuredApexIntensity
                Seo_Heavy                                   = inferred_Heavy.StandardErrorOfPrediction
                Params_Heavy                                = inferred_Heavy.EstimatedParams 
                Difference_SearchRT_FittedRT_Heavy          = searchRTMinusFittedRT_Heavy
                KLDiv_Observed_Theoretical_Heavy            = clusterComparison_Heavy.KLDiv_UnCorrected
                KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterComparison_Heavy.KLDiv_Corrected
                Correlation_Light_Heavy                     = corrLightHeavy
                QuantificationSource                        = QuantificationSource.PSM
                IsotopicPatternMz_Light                     = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.Mz)
                IsotopicPatternIntensity_Observed_Light     = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                IsotopicPatternIntensity_Corrected_Light    = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                RtTrace_Light                               = averagePSM.X_Xic 
                IntensityTrace_Observed_Light               = averagePSM.Y_Xic_uncorrected
                IntensityTrace_Corrected_Light              = averagePSM.Y_Xic
                IsotopicPatternMz_Heavy                     = clusterComparison_Heavy.PeakComparisons |> Array.map (fun x -> x.Mz)
                IsotopicPatternIntensity_Observed_Heavy     = clusterComparison_Heavy.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                IsotopicPatternIntensity_Corrected_Heavy    = clusterComparison_Heavy.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                RtTrace_Heavy                               = inferred_Heavy.X_Xic 
                IntensityTrace_Observed_Heavy               = inferred_Heavy.Y_Xic_uncorrected
                IntensityTrace_Corrected_Heavy              = inferred_Heavy.Y_Xic
                }
                |> Option.Some
            else
                let mz_Light          = Mass.toMZ (unlabledPeptide.Mass) (pepIon.Charge|> float)
                let inferred_Light       = quantifyInferredPeak getXIC identifyPeaks mz_Light averagePSM.WeightedAvgScanTime inferredScanTime
                let searchRTMinusFittedRT_Light = searchRTMinusFittedRtInferred inferredScanTime inferred_Light
                let clusterComparison_Light = comparePredictedAndMeasuredIsotopicCluster inferred_Light.X_Xic inferred_Light.Y_Xic inferred_Light.Y_Xic_uncorrected pepIon.Charge unlabledPeptide.BioSequence quantP.EstimatedParams.[1] mz_Light
                let corrLightHeavy        = calcCorrelation averagePSM.X_Xic quantP inferred_Light
                //let chart = saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
                //                    peakToQuantify.XData peakToQuantify.YData quantP.YPredicted inferred_Light.X_Xic inferred_Light.Y_Xic inferred_Light.xPeak inferred_Light.yFitted (*envelopeSumX envelopeSumY*) peaks clusterComparison_Target.PeakComparisons plotDirectory
                {
                StringSequence                              = pepIon.Sequence
                GlobalMod                                   = pepIon.GlobalMod
                Charge                                      = pepIon.Charge
                PepSequenceID                               = pepIon.PepSequenceID
                ModSequenceID                               = pepIon.ModSequenceID
                PrecursorMZ                                 = averagePSM.MeanPrecMz
                MeasuredMass                                = avgMass
                TheoMass                                    = labeledPeptide.Mass
                AbsDeltaMass                                = abs(avgMass-labeledPeptide.Mass)
                MeanPercolatorScore                         = averagePSM.MeanScore
                QValue                                      = bestQValue
                PEPValue                                    = bestPepValue
                ProteinNames                                = prots
                QuantMz_Light                               = mz_Light
                Quant_Light                                 = inferred_Light.Area
                MeasuredApex_Light                          = inferred_Light.MeasuredApexIntensity
                Seo_Light                                   = inferred_Light.StandardErrorOfPrediction
                Params_Light                                = inferred_Light.EstimatedParams 
                Difference_SearchRT_FittedRT_Light          = searchRTMinusFittedRT_Light
                KLDiv_Observed_Theoretical_Light            = clusterComparison_Light.KLDiv_UnCorrected
                KLDiv_CorrectedObserved_Theoretical_Light   = clusterComparison_Light.KLDiv_Corrected
                QuantMz_Heavy                               = averagePSM.MeanPrecMz
                Quant_Heavy                                 = quantP.Area
                MeasuredApex_Heavy                          = quantP.MeasuredApexIntensity
                Seo_Heavy                                   = quantP.StandardErrorOfPrediction
                Params_Heavy                                = quantP.EstimatedParams 
                Difference_SearchRT_FittedRT_Heavy          = searchRTMinusFittedRT
                KLDiv_Observed_Theoretical_Heavy            = clusterComparison_Target.KLDiv_UnCorrected
                KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterComparison_Target.KLDiv_Corrected
                Correlation_Light_Heavy                     = corrLightHeavy
                QuantificationSource                        = QuantificationSource.PSM
                IsotopicPatternMz_Light                     = clusterComparison_Light.PeakComparisons |> Array.map (fun x -> x.Mz)
                IsotopicPatternIntensity_Observed_Light     = clusterComparison_Light.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                IsotopicPatternIntensity_Corrected_Light    = clusterComparison_Light.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                RtTrace_Light                               = inferred_Light.X_Xic 
                IntensityTrace_Observed_Light               = inferred_Light.Y_Xic_uncorrected
                IntensityTrace_Corrected_Light              = inferred_Light.Y_Xic
                IsotopicPatternMz_Heavy                     = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.Mz)
                IsotopicPatternIntensity_Observed_Heavy     = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                IsotopicPatternIntensity_Corrected_Heavy    = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                RtTrace_Heavy                               = averagePSM.X_Xic 
                IntensityTrace_Observed_Heavy               = averagePSM.Y_Xic_uncorrected
                IntensityTrace_Corrected_Heavy              = averagePSM.Y_Xic
                }
                |> Option.Some
            with
            | ex ->
                
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None

        let lableFreeQuantification (pepIon:PeptideIon)  (psms:PSMStatisticsResult []) = 
            try
            if pepIon.GlobalMod <> 0 then None 
            else
            let bestQValue,bestPepValue,prots = psms |> Array.minBy (fun x -> x.QValue) |> fun x -> x.QValue, x.PEPValue,x.ProteinNames
            let unlabledPeptide = peptideLookUp pepIon.Sequence 0
            let psmsWithMatchedSums = countMatchedMasses unlabledPeptide psms 
            let ms2s = psmsWithMatchedSums |> Array.map (fun (psm,m) -> psm.ScanTime,m)
            let averagePSM = average getXIC psmsWithMatchedSums
            let avgMass = Mass.ofMZ (averagePSM.MeanPrecMz) (pepIon.Charge |> float)
            
            let peaks = 
                try
                    identifyPeaks averagePSM.X_Xic averagePSM.Y_Xic 
                with 
                | ex ->
                    logger.Trace (sprintf "Quantfailed: %A" ex)
                    [||]
            if Array.isEmpty peaks then 
                Chart.Point(averagePSM.X_Xic, averagePSM.Y_Xic)
                |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i_noPeaks" pepIon.Sequence pepIon.GlobalMod)
                |> Chart.withSize(1500.,800.)
                |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((pepIon.Sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + pepIon.GlobalMod.ToString() + "Ch" + pepIon.Charge.ToString() + "_notQuantified")|])
                None
            else
            let peakToQuantify = BioFSharp.Mz.Quantification.HULQ.getPeakBy peaks averagePSM.WeightedAvgScanTime
            let quantP = BioFSharp.Mz.Quantification.HULQ.quantifyPeak peakToQuantify
            let searchRTMinusFittedRT = searchRTMinusFittedRtTarget averagePSM.WeightedAvgScanTime quantP 
            let clusterComparison_Target = comparePredictedAndMeasuredIsotopicCluster averagePSM.X_Xic averagePSM.Y_Xic averagePSM.Y_Xic_uncorrected pepIon.Charge unlabledPeptide.BioSequence quantP.EstimatedParams.[1] averagePSM.MeanPrecMz
            //let chart = saveChart pepIon.Sequence pepIon.GlobalMod pepIon.Charge averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime
            //                    peakToQuantify.XData peakToQuantify.YData quantP.YPredicted  [||] [||] [||] [||] peaks clusterComparison_Target.PeakComparisons plotDirectory
            {
            StringSequence                              = pepIon.Sequence
            GlobalMod                                   = pepIon.GlobalMod
            Charge                                      = pepIon.Charge
            PepSequenceID                               = pepIon.PepSequenceID
            ModSequenceID                               = pepIon.ModSequenceID
            PrecursorMZ                                 = averagePSM.MeanPrecMz
            MeasuredMass                                = avgMass
            TheoMass                                    = unlabledPeptide.Mass
            AbsDeltaMass                                = abs(avgMass-unlabledPeptide.Mass)
            MeanPercolatorScore                         = averagePSM.MeanScore
            QValue                                      = bestQValue
            PEPValue                                    = bestPepValue
            ProteinNames                                = prots
            QuantMz_Light                               = averagePSM.MeanPrecMz
            Quant_Light                                 = quantP.Area
            MeasuredApex_Light                          = quantP.MeasuredApexIntensity
            Seo_Light                                   = quantP.StandardErrorOfPrediction
            Params_Light                                = quantP.EstimatedParams 
            Difference_SearchRT_FittedRT_Light          = searchRTMinusFittedRT
            KLDiv_Observed_Theoretical_Light            = clusterComparison_Target.KLDiv_UnCorrected
            KLDiv_CorrectedObserved_Theoretical_Light   = clusterComparison_Target.KLDiv_Corrected
            QuantMz_Heavy                               = nan
            Quant_Heavy                                 = nan
            MeasuredApex_Heavy                          = nan
            Seo_Heavy                                   = nan
            Params_Heavy                                = [||]
            Difference_SearchRT_FittedRT_Heavy          = nan
            KLDiv_Observed_Theoretical_Heavy            = nan
            KLDiv_CorrectedObserved_Theoretical_Heavy   = nan
            Correlation_Light_Heavy                     = nan
            QuantificationSource                        = QuantificationSource.PSM
            IsotopicPatternMz_Light                     = clusterComparison_Target.PeakComparisons|> Array.map (fun x -> x.Mz)
            IsotopicPatternIntensity_Observed_Light     = clusterComparison_Target.PeakComparisons|> Array.map (fun x -> x.MeasuredIntensity)
            IsotopicPatternIntensity_Corrected_Light    = clusterComparison_Target.PeakComparisons|> Array.map (fun x -> x.MeasuredIntensityCorrected)
            RtTrace_Light                               = averagePSM.X_Xic 
            IntensityTrace_Observed_Light               = averagePSM.Y_Xic_uncorrected
            IntensityTrace_Corrected_Light              = averagePSM.Y_Xic
            IsotopicPatternMz_Heavy                     = [||]
            IsotopicPatternIntensity_Observed_Heavy     = [||]
            IsotopicPatternIntensity_Corrected_Heavy    = [||]
            RtTrace_Heavy                               = [||]
            IntensityTrace_Observed_Heavy               = [||]
            IntensityTrace_Corrected_Heavy              = [||]
            }
            |> Option.Some
            with
            | ex ->
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None
        logger.Trace "init quantification functions:finished"
        
        logger.Trace "executing quantification"
        peptides        
        |> Array.groupBy (fun x -> 
            {
                Sequence             = x.StringSequence     
                GlobalMod            = x.GlobalMod    
                Charge               = x.Charge       
                ModSequenceID        = x.ModSequenceID
                PepSequenceID        = x.PepSequenceID
            }
            )
        |> Array.mapi (fun i (pepIon,psms) -> 
            if i % 100 = 0 then logger.Trace (sprintf "%i peptides quantified" i)
            if processParams.PerformLabeledQuantification then 
                labledQuantification pepIon psms
            else
                lableFreeQuantification pepIon psms
            )
        |> Array.filter Option.isSome
        |> Array.map (fun x -> x.Value)
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()
        memoryDB.Dispose()
        logger.Trace "executing quantification:finished"
        