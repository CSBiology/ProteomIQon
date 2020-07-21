namespace ProteomIQon

open System.IO
open Argu
open System.Data.SQLite
open System
open ProteomIQon.Core
open Core.MzIO
open Dto
open FSharp.Stats
open BioFSharp.Mz.Quantification
open BioFSharp.Mz
open FSharpAux.IO.SchemaReader
open FSharp.Plotly
open BioFSharp
open MzIO.Processing
open BioFSharp.Mz.SearchDB
open System.Data

module AlignmentBasedQuantification =

    type PeptideIon = 
        {
            Sequence             : string
            GlobalMod            : int
            Charge               : int
            ModSequenceID        : int
            PepSequenceID        : int
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
                    X_Xic                       = retData
                    Y_Xic                       = itzData
                    Y_Xic_uncorrected           = uncorrectedItzData
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

    ///
    let saveChart sequence globalMod ch (xXic:float[]) (yXic:float[]) avgScanTime (xToQuantify:float[]) (ypToQuantify:float[]) (fitY:float[])
            (xXicInferred:float[]) (yXicinferred:float[]) (xInferred:float[]) (inferredFit:float[]) (*(xEnvelopeSum:float[]) (yEnvelopeSum:float[])*) (pattern:PeakComparison []) plotDirectory =
        let xic = 
            [
            Chart.Point(xXic, yXic)                     |> Chart.withTraceName "Target XIC"
            Chart.Point([avgScanTime],[1.])             |> Chart.withTraceName "Inferred Scan Time"
            Chart.Point((xToQuantify), (ypToQuantify))  |> Chart.withTraceName "Identified Target Peak"
            Chart.Line(xToQuantify,fitY)                |> Chart.withTraceName "Fit of target Peak"
            Chart.Point(xXicInferred, yXicinferred)     |> Chart.withTraceName "Inferred XIC"
            Chart.Line(xInferred,inferredFit)           |> Chart.withTraceName "Fit of inferred Peak"

            ]
            |> Chart.Combine
        let pattern = 
            [
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensity))          |> Chart.withTraceName "Measured"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.MeasuredIntensityCorrected)) |> Chart.withTraceName "Measured Corrected"
            Chart.Point(pattern |> Array.map (fun x -> x.Mz), pattern |> Array.map (fun x -> x.PredictedRelFrequency))      |> Chart.withTraceName "Predicted Relative Frequency"
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
        
    ///
    let searchRTMinusFittedRt searchRT fit = 
        try
            searchRT - fit.EstimatedParams.[1] 
        with
        | _ -> nan

    /// Aims to provide the best scan time estimate when quanitfying a inferred peak.
    let chooseScanTime maxDiff searchRTMinusFittedRT initialScanTime quantP = 
        if Array.isEmpty quantP.EstimatedParams then
            initialScanTime 
        elif abs searchRTMinusFittedRT > maxDiff then
            initialScanTime 
        else
            quantP.EstimatedParams.[1]

    ///
    let calcCorrelation (xValues:float []) (targetPeak:InferredPeak) (inferredPeak:InferredPeak) = 
        let getValue (model:HULQ.PeakModel) estParams =
            match model with
            | HULQ.PeakModel.Gaussian m -> 
                m.GetFunctionValue (vector estParams)
            | HULQ.PeakModel.EMG m -> 
                m.GetFunctionValue (vector estParams)
        match targetPeak.Model, inferredPeak.Model with 
        | Some q , Some i -> 
            let xValuesBW = 
                [|for i = 1 to xValues.Length-1 do abs(xValues.[i] - xValues.[i-1])|]
                |> Array.median
            let xValues = [|xValues.[0] .. xValuesBW/2. .. xValues.[xValues.Length-1]|]
            let fQ = getValue q targetPeak.EstimatedParams
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
                peak.Mz,peak.Intensity,list |> Array.sumBy snd 
                )
            |> Array.map (fun (mz,measuredIntensity,predictedRelFrequency) -> 
                    {Mz=mz;MeasuredIntensity = measuredIntensity;MeasuredIntensityCorrected= measuredIntensity - baseLineCorrectionF;PredictedRelFrequency= predictedRelFrequency}
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
    let quantifyPeptides (processParams:Domain.AlignmentBasedQuantificationParams) (outputDir:string) (cn:SQLiteConnection) (instrumentOutput:string) (scoredPSMs:string)  =

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
        let dBParams = SearchDB'.getSDBParams memoryDB
        let peptideLookUp = SearchDB'.getThreadSafePeptideLookUpFromFileBySequenceAndGMod memoryDB dBParams
        let calcIonSeries aal = Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal
        logger.Trace "Get peptide lookUp function: finished"
        // initialize Reader and Transaction
        logger.Trace "Init connection to mass spectrum data."
        let inReader = Core.MzIO.Reader.getReader instrumentOutput
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader
        let inTr = inReader.BeginTransaction()

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
            Csv.CsvReader<AlignmentResult>(SchemaMode=Csv.Fill).ReadFile(scoredPSMs,'\t',false,1)
            |> Array.ofSeq
        logger.Trace "Read scored PSMs:finished"

        logger.Trace "Estimate ms1 mz accuracy"
        ///
        let ms1AccuracyEstimate = 0.01
            //peptides
            //|> Seq.stDevBy (fun x -> abs(x.PrecursorMZ - Mass.toMZ x.TheoMass (float x.Charge)) )
        logger.Trace "Estimate ms1 mz accuracy:finished"           
        
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
        let labledQuantification (psms:AlignmentResult) = 
            try
            let unlabledPeptide = peptideLookUp psms.StringSequence 0
            let labeledPeptide  = peptideLookUp psms.StringSequence 1
            let targetPeptide = if psms.GlobalMod = 0 then unlabledPeptide else labeledPeptide            
            let targetMz = Mass.toMZ (targetPeptide.Mass) (psms.Charge|> float)
            let targetQuant = quantifyInferredPeak getXIC identifyPeaks targetMz psms.PredictedScanTime psms.PredictedScanTime
            if Array.isEmpty targetQuant.EstimatedParams then 
                Chart.Point(targetQuant.X_Xic, targetQuant.Y_Xic)
                |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i_noPeaks" psms.StringSequence psms.GlobalMod)
                |> Chart.withSize(1500.,800.)
                |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((psms.StringSequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + psms.GlobalMod.ToString() + "Ch" + psms.Charge.ToString() + "_notQuantified")|])
                None
            else
            let searchRTMinusFittedRT = searchRTMinusFittedRt psms.PredictedScanTime targetQuant
            let inferredScanTime = chooseScanTime processParams.XicExtraction.ScanTimeWindow searchRTMinusFittedRT psms.PredictedScanTime targetQuant
            let clusterComparison_Target = comparePredictedAndMeasuredIsotopicCluster targetQuant.X_Xic targetQuant.Y_Xic targetQuant.Y_Xic_uncorrected psms.Charge targetPeptide.BioSequence targetQuant.EstimatedParams.[1] targetMz
            if psms.GlobalMod = 0 then
                let mz_Heavy = Mass.toMZ (labeledPeptide.Mass) (psms.Charge|> float)
                let inferred_Heavy = quantifyInferredPeak getXIC identifyPeaks mz_Heavy psms.PredictedScanTime inferredScanTime
                let searchRTMinusFittedRT_Heavy = searchRTMinusFittedRt inferredScanTime inferred_Heavy
                let clusterComparison_Heavy = comparePredictedAndMeasuredIsotopicCluster inferred_Heavy.X_Xic inferred_Heavy.Y_Xic inferred_Heavy.Y_Xic_uncorrected psms.Charge labeledPeptide.BioSequence targetQuant.EstimatedParams.[1] mz_Heavy
                let corrLightHeavy  = calcCorrelation targetQuant.X_Xic targetQuant inferred_Heavy  
                let chart = saveChart psms.StringSequence psms.GlobalMod psms.Charge targetQuant.X_Xic targetQuant.Y_Xic psms.PredictedScanTime
                                    targetQuant.xPeak targetQuant.yFitted targetQuant.yFitted 
                                    inferred_Heavy.X_Xic inferred_Heavy.Y_Xic inferred_Heavy.xPeak inferred_Heavy.yFitted 
                                    clusterComparison_Target.PeakComparisons plotDirectory
                {
                StringSequence                              = psms.StringSequence
                GlobalMod                                   = psms.GlobalMod
                Charge                                      = psms.Charge
                PepSequenceID                               = psms.PepSequenceID
                ModSequenceID                               = psms.ModSequenceID
                PrecursorMZ                                 = psms.Mz
                MeasuredMass                                = targetPeptide.Mass
                TheoMass                                    = targetPeptide.Mass
                AbsDeltaMass                                = nan
                MeanPercolatorScore                         = 0.
                QValue                                      = 0.
                PEPValue                                    = 0.
                ProteinNames                                = psms.ProteinNames
                QuantMz_Light                               = psms.Mz
                Quant_Light                                 = targetQuant.Area
                MeasuredApex_Light                          = targetQuant.MeasuredApexIntensity
                Seo_Light                                   = targetQuant.StandardErrorOfPrediction
                Params_Light                                = targetQuant.EstimatedParams            
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
                QuantificationSource                        = QuantificationSource.Alignment
                IsotopicPatternMz_Light                     = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.Mz)
                IsotopicPatternIntensity_Observed_Light     = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                IsotopicPatternIntensity_Corrected_Light    = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                RtTrace_Light                               = targetQuant.X_Xic 
                IntensityTrace_Observed_Light               = targetQuant.Y_Xic_uncorrected
                IntensityTrace_Corrected_Light              = targetQuant.Y_Xic
                IsotopicPatternMz_Heavy                     = clusterComparison_Heavy.PeakComparisons |> Array.map (fun x -> x.Mz)
                IsotopicPatternIntensity_Observed_Heavy     = clusterComparison_Heavy.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                IsotopicPatternIntensity_Corrected_Heavy    = clusterComparison_Heavy.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                RtTrace_Heavy                               = inferred_Heavy.X_Xic 
                IntensityTrace_Observed_Heavy               = inferred_Heavy.Y_Xic_uncorrected
                IntensityTrace_Corrected_Heavy              = inferred_Heavy.Y_Xic
                }
                |> Option.Some
            else
                let mz_Light = Mass.toMZ (unlabledPeptide.Mass) (psms.Charge|> float)
                let inferred_Light = quantifyInferredPeak getXIC identifyPeaks mz_Light psms.PredictedScanTime inferredScanTime
                let searchRTMinusFittedRT_Light = searchRTMinusFittedRt inferredScanTime inferred_Light
                let clusterComparison_Light = comparePredictedAndMeasuredIsotopicCluster inferred_Light.X_Xic inferred_Light.Y_Xic inferred_Light.Y_Xic_uncorrected psms.Charge unlabledPeptide.BioSequence targetQuant.EstimatedParams.[1] mz_Light 
                let corrLightHeavy  = calcCorrelation targetQuant.X_Xic targetQuant inferred_Light  
                let chart = saveChart psms.StringSequence psms.GlobalMod psms.Charge targetQuant.X_Xic targetQuant.Y_Xic psms.PredictedScanTime
                                    targetQuant.xPeak targetQuant.yFitted targetQuant.yFitted 
                                    inferred_Light.X_Xic inferred_Light.Y_Xic inferred_Light.xPeak inferred_Light.yFitted 
                                    clusterComparison_Target.PeakComparisons plotDirectory
                {
                StringSequence                              = psms.StringSequence
                GlobalMod                                   = psms.GlobalMod
                Charge                                      = psms.Charge
                PepSequenceID                               = psms.PepSequenceID
                ModSequenceID                               = psms.ModSequenceID
                PrecursorMZ                                 = psms.Mz
                MeasuredMass                                = targetPeptide.Mass
                TheoMass                                    = targetPeptide.Mass
                AbsDeltaMass                                = nan
                MeanPercolatorScore                         = 0.
                QValue                                      = 0.
                PEPValue                                    = 0.
                ProteinNames                                = psms.ProteinNames
                QuantMz_Light                               = mz_Light
                Quant_Light                                 = inferred_Light.Area
                MeasuredApex_Light                          = inferred_Light.MeasuredApexIntensity
                Seo_Light                                   = inferred_Light.StandardErrorOfPrediction
                Params_Light                                = inferred_Light.EstimatedParams 
                Difference_SearchRT_FittedRT_Light          = searchRTMinusFittedRT_Light
                KLDiv_Observed_Theoretical_Light            = clusterComparison_Light.KLDiv_UnCorrected
                KLDiv_CorrectedObserved_Theoretical_Light   = clusterComparison_Light.KLDiv_Corrected
                QuantMz_Heavy                               = psms.Mz
                Quant_Heavy                                 = targetQuant.Area
                MeasuredApex_Heavy                          = targetQuant.MeasuredApexIntensity
                Seo_Heavy                                   = targetQuant.StandardErrorOfPrediction
                Params_Heavy                                = targetQuant.EstimatedParams 
                Difference_SearchRT_FittedRT_Heavy          = searchRTMinusFittedRT
                KLDiv_Observed_Theoretical_Heavy            = clusterComparison_Target.KLDiv_UnCorrected
                KLDiv_CorrectedObserved_Theoretical_Heavy   = clusterComparison_Target.KLDiv_Corrected
                Correlation_Light_Heavy                     = corrLightHeavy
                QuantificationSource                        = QuantificationSource.Alignment
                IsotopicPatternMz_Light                     = clusterComparison_Light.PeakComparisons |> Array.map (fun x -> x.Mz)
                IsotopicPatternIntensity_Observed_Light     = clusterComparison_Light.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                IsotopicPatternIntensity_Corrected_Light    = clusterComparison_Light.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                RtTrace_Light                               = inferred_Light.X_Xic 
                IntensityTrace_Observed_Light               = inferred_Light.Y_Xic_uncorrected
                IntensityTrace_Corrected_Light              = inferred_Light.Y_Xic
                IsotopicPatternMz_Heavy                     = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.Mz)
                IsotopicPatternIntensity_Observed_Heavy     = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensity)
                IsotopicPatternIntensity_Corrected_Heavy    = clusterComparison_Target.PeakComparisons |> Array.map (fun x -> x.MeasuredIntensityCorrected)
                RtTrace_Heavy                               = targetQuant.X_Xic 
                IntensityTrace_Observed_Heavy               = targetQuant.Y_Xic_uncorrected
                IntensityTrace_Corrected_Heavy              = targetQuant.Y_Xic

                }
                |> Option.Some
            with
            | ex ->
             
                logger.Trace (sprintf "Quantfailed: %A" ex)
                Option.None

        let lableFreeQuantification (psms:AlignmentResult) = 
            try
            if psms.GlobalMod <> 0 then 
                None
            else
                let unlabledPeptide = peptideLookUp psms.StringSequence 0
                let labeledPeptide  = peptideLookUp psms.StringSequence 1
                let targetPeptide = if psms.GlobalMod = 0 then unlabledPeptide else labeledPeptide            
                let targetMz = Mass.toMZ (targetPeptide.Mass) (psms.Charge|> float)
                let targetQuant = quantifyInferredPeak getXIC identifyPeaks targetMz psms.PredictedScanTime psms.PredictedScanTime
                if Array.isEmpty targetQuant.EstimatedParams then 
                    Chart.Point(targetQuant.X_Xic, targetQuant.Y_Xic)
                    |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i_noPeaks" psms.StringSequence psms.GlobalMod)
                    |> Chart.withSize(1500.,800.)
                    |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((psms.StringSequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + psms.GlobalMod.ToString() + "Ch" + psms.Charge.ToString() + "_notQuantified")|])
                    None
                else
                let searchRTMinusFittedRT = searchRTMinusFittedRt psms.PredictedScanTime targetQuant 
                let clusterComparison_Target = comparePredictedAndMeasuredIsotopicCluster targetQuant.X_Xic targetQuant.Y_Xic targetQuant.Y_Xic_uncorrected psms.Charge targetPeptide.BioSequence targetQuant.EstimatedParams.[1] targetMz
                let chart = saveChart psms.StringSequence psms.GlobalMod psms.Charge targetQuant.X_Xic targetQuant.Y_Xic psms.PredictedScanTime
                                    targetQuant.xPeak targetQuant.yFitted targetQuant.yFitted 
                                    [||] [||] [||] [||] 
                                    clusterComparison_Target.PeakComparisons plotDirectory
                {
                StringSequence                              = psms.StringSequence
                GlobalMod                                   = psms.GlobalMod
                Charge                                      = psms.Charge
                PepSequenceID                               = psms.PepSequenceID
                ModSequenceID                               = psms.ModSequenceID
                PrecursorMZ                                 = psms.Mz
                MeasuredMass                                = targetPeptide.Mass
                TheoMass                                    = targetPeptide.Mass
                AbsDeltaMass                                = nan
                MeanPercolatorScore                         = 0.
                QValue                                      = 0.
                PEPValue                                    = 0.
                ProteinNames                                = psms.ProteinNames
                QuantMz_Light                               = psms.Mz
                Quant_Light                                 = targetQuant.Area
                MeasuredApex_Light                          = targetQuant.MeasuredApexIntensity
                Seo_Light                                   = targetQuant.StandardErrorOfPrediction
                Params_Light                                = targetQuant.EstimatedParams            
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
                QuantificationSource                        = QuantificationSource.Alignment
                IsotopicPatternMz_Light                     = clusterComparison_Target.PeakComparisons|> Array.map (fun x -> x.Mz)
                IsotopicPatternIntensity_Observed_Light     = clusterComparison_Target.PeakComparisons|> Array.map (fun x -> x.MeasuredIntensity)
                IsotopicPatternIntensity_Corrected_Light    = clusterComparison_Target.PeakComparisons|> Array.map (fun x -> x.MeasuredIntensityCorrected)
                RtTrace_Light                               = targetQuant.X_Xic 
                IntensityTrace_Observed_Light               = targetQuant.Y_Xic_uncorrected
                IntensityTrace_Corrected_Light              = targetQuant.Y_Xic
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

        |> Array.choose (fun psms -> 
            if processParams.PerformLabeledQuantification then 
                labledQuantification psms
            else
                lableFreeQuantification psms
            )
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        logger.Trace "executing quantification:finished"
        