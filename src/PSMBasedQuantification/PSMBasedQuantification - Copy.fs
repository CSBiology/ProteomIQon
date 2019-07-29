namespace ProteomIQon

open System.IO
open Argu
open System.Data.SQLite
open System
open ProteomIQon.Core
open Core.MzLite
open Dto
open FSharp.Stats
open BioFSharp.Mz.Quantification
open BioFSharp.Mz
open FSharpAux.IO.SchemaReader
open FSharp.Plotly
open BioFSharp
open BioFSharp.Mz.Quantification

//open FSharp.pl

module PSMBasedQuantification_old = 

    ///
    let setIndexOnModSequenceAndGlobalMod (cn:SQLiteConnection) =
        let querystring = "CREATE INDEX IF NOT EXISTS SequenceAndGlobalModIndex ON ModSequence (Sequence,GlobalMod)"
        let cmd = new SQLiteCommand(querystring, cn)    
        cmd.ExecuteNonQuery()

    /// Prepares statement to select a ModSequence entry by Massrange (Between selected Mass -/+ the selected toleranceWidth)
    let prepareSelectMassByModSequenceAndGlobalMod (cn:SQLiteConnection) =
        let querystring = "SELECT RealMass FROM ModSequence WHERE Sequence=@sequence AND GlobalMod=@globalMod"
        let cmd = new SQLiteCommand(querystring, cn) 
        cmd.Parameters.Add("@sequence", Data.DbType.String) |> ignore
        cmd.Parameters.Add("@globalMod", Data.DbType.Int32) |> ignore
        (fun (sequence:string) (globalMod:int) ->
        cmd.Parameters.["@sequence"].Value  <- sequence
        cmd.Parameters.["@globalMod"].Value <- globalMod
        use reader = cmd.ExecuteReader()            
        match reader.Read() with 
        | true  -> Some (reader.GetDouble(0))
        | false -> Option.None)

    type PathPair =
        {
            Name : string
            MZInputPath : string
            PeptideInputPath: string
        }

    let createPathPair name mzliteInputPath peptideInputPath=
        {
            Name = name
            MZInputPath = mzliteInputPath
            PeptideInputPath = peptideInputPath
        }
   
    type averagePSM = {
        MeanPrecMz   : float 
        MeanScanTime : float
        WeightedAvgScanTime:float
        MeanScore   : float
        X_Xic         : float []
        Y_Xic         : float []
        }

    let createAveragePSM meanPrecMz meanScanTime weightedAvgScanTime meanScore xXic yXic = {
        MeanPrecMz    = meanPrecMz   
        MeanScanTime  = meanScanTime 
        WeightedAvgScanTime= weightedAvgScanTime
        MeanScore = meanScore
        X_Xic         = xXic         
        Y_Xic         = yXic         
        }

    let weightedMean (weights:seq<'T>) (items:seq<'T>) =
        let sum,n = Seq.fold2 (fun (sum,n) w i -> w*i+sum,n + w ) (0.,0.) weights items 
        sum / n 


    let substractBaseLine (baseLineParams:Domain.BaseLineCorrection) (yData:float []) = 
        if yData.Length > 300 then 
            printfn "xic Length > 300" 
            yData 
        else
        let baseLine = FSharp.Stats.Signal.Baseline.baselineAls baseLineParams.MaxIterations baseLineParams.Lambda baseLineParams.P yData |> Array.ofSeq
        Array.map2 (fun y b -> 
                       let c = y - b
                       if c < 0. then 0. else c
                   ) yData baseLine

    let initGetProcessedXIC (baseLineCorrection:Domain.BaseLineCorrection option) reader idx scanTimeWindow mzWindow_Da meanScanTime meanPrecMz =
        let rtQuery = Query.createRangeQuery meanScanTime scanTimeWindow
        let mzQuery = Query.createRangeQuery meanPrecMz mzWindow_Da
        let retData',itzData' =
            Query.getXIC reader idx rtQuery mzQuery
            |> Array.map (fun p -> p.Rt , p.Intensity)
            |> Array.unzip
        match baseLineCorrection with
        | Some baseLineParams -> 
            retData', substractBaseLine baseLineParams itzData' 
        | None -> 
            retData',itzData'


    let average getXic (psms:(PSMStatisticsResult*float) []) =
            let meanPrecMz          = psms |> Seq.meanBy (fun (x,scanTime) -> x.PrecursorMZ)
            let meanScanTime        = psms |> Seq.meanBy (fun (x,scanTime) -> scanTime)
            let (retData,itzData)   = getXic meanScanTime meanPrecMz
            let meanScore = psms |> Seq.averageBy (fun (x,scanTime) -> x.PercolatorScore)
            let weightedAvgScanTime =
                let scanTimes = 
                    psms 
                    |> Seq.map snd
                let weights =
                    scanTimes
                    |> Seq.map (FSharp.Stats.Signal.PeakDetection.idxOfClosestPeakBy retData itzData)
                    |> Seq.map (fun idx -> itzData.[idx])
                    |> Seq.map (fun x -> if x <= 0. then 1. else x)
                weightedMean weights scanTimes 
            createAveragePSM meanPrecMz meanScanTime weightedAvgScanTime meanScore retData itzData

    let quantifyBy minSNR polOrder windowSize getXic targetMz targetScanTime =
        let (retData,itzData)   =
            getXic targetScanTime targetMz
        let peaks          = Signal.PeakDetection.SecondDerivative.getPeaks minSNR polOrder windowSize retData itzData
        let peakToQuantify = BioFSharp.Mz.Quantification.MyQuant.getPeakBy peaks targetScanTime
        let quantP         = BioFSharp.Mz.Quantification.MyQuant.quantifyPeak peakToQuantify  
        quantP,retData,itzData, peakToQuantify.XData

    let saveChart windowWidth sequence globalMod ch (xXic:float[]) (yXic:float[]) ms2s avgScanTime (xToQuantify:float[]) (ypToQuantify:float[]) (fitY:float[]) 
            (xXicInferred:float[]) (yXicinferred:float[]) (xInferred:float[]) (inferredFit:float[]) plotDirectory =
        [
        Chart.Point(xXic, yXic)                     |> Chart.withTraceName "Target XIC"
        Chart.Point(ms2s)                           |> Chart.withTraceName "MS2s with scores"
        Chart.Point([avgScanTime],[1.])             |> Chart.withTraceName "Weighted Mean of Ms2 scan times"
        Chart.Point((xToQuantify), (ypToQuantify))  |> Chart.withTraceName "Identified Target Peak"
        Chart.Line(xToQuantify,fitY)                |> Chart.withTraceName "Fit of target Peak"
        Chart.Point(xXicInferred, yXicinferred)     |> Chart.withTraceName "Inferred XIC"
        Chart.Line(xInferred,inferredFit)           |> Chart.withTraceName "Fit of inferred Peak"

        ]
        |> Chart.Combine
        |> Chart.withTitle(sprintf "Sequence= %s,globalMod = %i, WindowWidth = %i" sequence globalMod windowWidth)
        |> Chart.withSize(1500.,800.)
        |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory; ((sequence |> String.filter (fun x -> x <> '*')) + "_GMod_" + globalMod.ToString() + "Ch" + ch.ToString())|])


    ///
    type Result = 
        {
        StringSequence               : string;
        GlobalMod                    : int;
        Charge                       : int;
        PrecursorMZ                  : float;
        MeasuredMass                 : float; 
        TheoMass                     : float;
        AbsDeltaMass                 : float;
        MeanPercolatorScore          : float;
        QValue                       : float;
        PEPValue                     : float;
        ProteinNames                 : string;
        N14QuantMz                   : float
        N14Quant                     : float
        N14Seo                       : float
        N14Params                    : string
        N15QuantMz                   : float
        N15Quant                     : float
        N15Seo                       : float
        N15Params                    : string
        N15Minus1QuantMz             : float
        N15Minus1Quant               : float
        N15Minus1Seo                 : float
        N15Minus1Params              : string
        }

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
    let quantifyPeptides (processParams:Domain.QuantificationParams) (outputDir:string) (cn:SQLiteConnection) (instrumentOutput:string) (scoredPSMs:string)  =
        printfn "Now performing Quantification using: %s and %s, Results will be written to: %s" instrumentOutput scoredPSMs outputDir

        // initialize Reader and Transaction
        let outFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".quant"
            Path.Combine [|outputDir;fileName|]
        printfn "outFilePath:%s" outFilePath

        //
        let plotDirectory = 
            let fileName = sprintf "%s_plots" (Path.GetFileNameWithoutExtension instrumentOutput) 
            let path = Path.Combine [|outputDir;fileName|]
            if System.IO.Directory.Exists path then
                path
            else 
                System.IO.Directory.CreateDirectory path |> ignore
                path
        printfn "plotDirectory:%s" plotDirectory
        
        printfn "Copy peptide DB into Memory"
        let memoryDB = SearchDB.copyDBIntoMemory cn 
        printfn "Copy peptide DB into Memory: finished"
        
        printfn  "Get peptide lookUp function"
        let massLookUp = prepareSelectMassByModSequenceAndGlobalMod memoryDB 
        printfn  "Get peptide lookUp function: finished"

        // initialize Reader and Transaction
        printfn "Init connection to mass spectrum data." 
        let inReader = Core.MzLite.Reader.getReader instrumentOutput  
        let inRunID  = Core.MzLite.Reader.getDefaultRunID inReader
        let inTr = inReader.BeginTransaction()

        printfn "Create RetentionTime index"
        let retTimeIdxed = Query.getMS1RTIdx inReader inRunID
        printfn "Create RetentionTime index:finished"
            
        // initialize Reader and Transaction
        printfn "Read scored PSMs."
        ///
        let peptides =
            Csv.CsvReader<PSMStatisticsResult>(SchemaMode=Csv.Fill).ReadFile(scoredPSMs,'\t',false,1)
            |> Array.ofSeq
        
        let getXIC = initGetProcessedXIC processParams.BaseLineCorrection inReader retTimeIdxed processParams.XicExtraction.ScanTimeWindow processParams.XicExtraction.MzWindow_Da 
        

        let noiseAutoCorrelationMedian = 
            match processParams.XicExtraction.WindowSize with 
            | Domain.WindowSize.Fixed _  -> None
            | Domain.WindowSize.Estimate ->
                printfn "Estimate noise autocorrelation"
                try
                let noiseAutoCorr = 
                    let gpeptides =
                        peptides
                        |> Array.groupBy (fun x -> x.StringSequence,x.Charge,x.GlobalMod)
                    let n = 
                        if gpeptides.Length < 260 then peptides.Length-1
                        else 250
                    gpeptides
                    |> Array.shuffleFisherYates
                    |> Array.take n
                    |> Array.choose (fun ((sequence,ch,globMod),psms) ->
                                    try
                                    printfn "sequence = %s,ch = %i,globMod = %i " sequence ch globMod
                                    let psmsWithScanTime = psms |> Array.map (fun x -> x, MassSpectrum.getScanTime (inReader.ReadMassSpectrum(x.PSMId)))
                                    let ms2s = psmsWithScanTime |> Array.map (fun (psm, scanTime) -> scanTime,psm.PercolatorScore)
                                    printfn "quantify target"
                                    let averagePSM = average getXIC psmsWithScanTime
                                    let avgMass = Mass.ofMZ (averagePSM.MeanPrecMz) (ch |> float)
                                    let peaks          = Signal.PeakDetection.SecondDerivative.getPeaks 0.1 2 11 averagePSM.X_Xic averagePSM.Y_Xic
                                    let NoNoise = peaks |> Array.map (fun x -> x.XData) |> Array.concat |> Set.ofArray
                                    let noiseArr = 
                                        Array.zip averagePSM.X_Xic averagePSM.Y_Xic
                                        |> Array.filter (fun (ret,intensity) -> NoNoise.Contains ret |> not)
                                        |> Array.map snd 
                                    
                                    let corr = Correlation.Vector.autoCorrelation 1 (noiseArr |> vector)                        
                                    if nan.Equals corr then Some 0. else Some corr 
                                    with
                                    | ex ->
                                        printfn "%A" ex
                                        None
                                 )
                let medianAutoCorr = Seq.median (noiseAutoCorr |> Array.filter (fun x -> nan.Equals(x) |> not) |> Array.sort)
                printfn "Estimate noise autocorrelation:finished"
                Chart.BoxPlot noiseAutoCorr
                |> Chart.SaveHtmlAs(Path.Combine[|plotDirectory;"NoiseAutoCorrelation"|])
                Some (medianAutoCorr)
                with
                | ex -> 
                    printfn "%A" ex
                    Some 0.5
                                    
        peptides 
        |> Array.groupBy (fun x -> x.StringSequence,x.Charge,x.GlobalMod)
        //|> fun x -> x.[3500..4500]
        |> Array.mapi (fun i ((sequence,ch,globMod),psms) ->
                        try
                        printfn "%i, sequence = %s,ch = %i,globMod =%i " i sequence ch globMod
                        let bestQValue,bestPepValue,prots = psms |> Array.minBy (fun x -> x.QValue) |> fun x -> x.QValue, x.PEPValue,x.ProteinNames
                        let psmsWithScanTime = psms |> Array.map (fun x -> x, MassSpectrum.getScanTime (inReader.ReadMassSpectrum(x.PSMId)))
                        let ms2s = psmsWithScanTime |> Array.map (fun (psm, scanTime) -> scanTime,psm.PercolatorScore)
                        printfn "quantify target"
                        let averagePSM = average getXIC psmsWithScanTime
                        let avgMass = Mass.ofMZ (averagePSM.MeanPrecMz) (ch |> float)
                        let windowWidth = 
                            match processParams.XicExtraction.WindowSize with 
                            | Domain.WindowSize.Fixed w  -> w
                            | Domain.WindowSize.Estimate ->
                                optimizeWindowWidth processParams.XicExtraction.PolynomOrder [|5 .. 2 .. 60|] noiseAutoCorrelationMedian.Value averagePSM.Y_Xic
                        let peaks          = Signal.PeakDetection.SecondDerivative.getPeaks processParams.XicExtraction.MinSNR processParams.XicExtraction.PolynomOrder windowWidth averagePSM.X_Xic averagePSM.Y_Xic
                        let peakToQuantify = BioFSharp.Mz.Quantification.MyQuant.getPeakBy peaks averagePSM.WeightedAvgScanTime
                        let quantP = BioFSharp.Mz.Quantification.MyQuant.quantifyPeak peakToQuantify 
                        let searchScanTime = 
                            if quantP.EstimatedParams |> Array.isEmpty then
                                averagePSM.WeightedAvgScanTime
                            else 
                                quantP.EstimatedParams.[1] 

                        let unlabledMass   = massLookUp sequence 0
                        let labeledMass    = massLookUp sequence 1
                        if globMod = 0 then 
                            let n15mz          = Mass.toMZ (labeledMass.Value) (ch|> float)
                            printfn "quantify inferred"
                            let n15Quant,rt,itz,rtP       = quantifyBy processParams.XicExtraction.MinSNR processParams.XicExtraction.PolynomOrder windowWidth  getXIC n15mz searchScanTime
                            let n15Minus1Mz    = n15mz - (Mass.Table.NMassInU / (ch|> float))
                            printfn "quantify n15Minus 1"
                            let n15Minus1Quant,_,_,_ = quantifyBy processParams.XicExtraction.MinSNR processParams.XicExtraction.PolynomOrder windowWidth getXIC n15Minus1Mz searchScanTime
                            let chart = saveChart windowWidth sequence globMod ch averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime 
                                                peakToQuantify.XData peakToQuantify.YData quantP.YPredicted rt itz rtP n15Quant.YPredicted plotDirectory
                    
                            {
                            StringSequence            = sequence   
                            GlobalMod                 = globMod  
                            Charge                    = ch   
                            PrecursorMZ               = averagePSM.MeanPrecMz  
                            MeasuredMass              = avgMass
                            TheoMass                  = unlabledMass.Value  
                            AbsDeltaMass              = abs(avgMass-unlabledMass.Value)  
                            MeanPercolatorScore       = averagePSM.MeanScore  
                            QValue                    = bestQValue
                            PEPValue                  = bestPepValue
                            ProteinNames              = prots  
                            N14QuantMz                = averagePSM.MeanPrecMz   
                            N14Quant                  = quantP.Area  
                            N14Seo                    = quantP.StandardErrorOfPrediction
                            N14Params                 = quantP.EstimatedParams |> Array.map (fun x -> x.ToString()) |> String.concat "; "
                            N15QuantMz                = n15mz  
                            N15Quant                  = n15Quant.Area   
                            N15Seo                    = n15Quant.StandardErrorOfPrediction 
                            N15Params                 = n15Quant.EstimatedParams |> Array.map (fun x -> x.ToString()) |> String.concat "; "
                            N15Minus1QuantMz          = n15Minus1Mz
                            N15Minus1Quant            = n15Minus1Quant.Area   
                            N15Minus1Seo              = n15Minus1Quant.StandardErrorOfPrediction 
                            N15Minus1Params           = n15Minus1Quant.EstimatedParams |> Array.map (fun x -> x.ToString()) |> String.concat "; "
                            } 
                            |> Option.Some

                        else
                            let n14mz          = Mass.toMZ (unlabledMass.Value) (ch|> float)
                            printfn "quantify inferred"
                            let n14Quant,rt,itz,rtP      = quantifyBy processParams.XicExtraction.MinSNR processParams.XicExtraction.PolynomOrder windowWidth getXIC n14mz searchScanTime
                            printfn "quantify n15Minus 1"
                            let n15Minus1Mz    = averagePSM.MeanPrecMz - (Mass.Table.NMassInU / (ch|> float))
                            let n15Minus1Quant,_,_,_ = quantifyBy processParams.XicExtraction.MinSNR processParams.XicExtraction.PolynomOrder windowWidth getXIC n15Minus1Mz searchScanTime
                            let chart = saveChart windowWidth sequence globMod ch averagePSM.X_Xic averagePSM.Y_Xic ms2s averagePSM.WeightedAvgScanTime 
                                                peakToQuantify.XData peakToQuantify.YData quantP.YPredicted rt itz rtP n14Quant.YPredicted plotDirectory
                            
                            {
                            StringSequence            = sequence   
                            GlobalMod                 = globMod  
                            Charge                    = ch   
                            PrecursorMZ               = averagePSM.MeanPrecMz  
                            MeasuredMass              = avgMass
                            TheoMass                  = labeledMass.Value  
                            AbsDeltaMass              = abs(avgMass-labeledMass.Value)  
                            MeanPercolatorScore       = averagePSM.MeanScore  
                            QValue                    = bestQValue
                            PEPValue                  = bestPepValue
                            ProteinNames              = prots  
                            N14QuantMz                = n14mz  
                            N14Quant                  = n14Quant.Area   
                            N14Seo                    = n14Quant.StandardErrorOfPrediction 
                            N14Params                 = n14Quant.EstimatedParams |> Array.fold (fun acc x -> acc + " " + x.ToString() + ";") ""
                            N15QuantMz                = averagePSM.MeanPrecMz   
                            N15Quant                  = quantP.Area  
                            N15Seo                    = quantP.StandardErrorOfPrediction
                            N15Params                 = quantP.EstimatedParams |> Array.fold (fun acc x -> acc + " " + x.ToString() + ";") ""    
                            N15Minus1QuantMz          = n15Minus1Mz
                            N15Minus1Quant            = n15Minus1Quant.Area   
                            N15Minus1Seo              = n15Minus1Quant.StandardErrorOfPrediction 
                            N15Minus1Params           = n15Minus1Quant.EstimatedParams |> Array.fold (fun acc x -> acc + " " + x.ToString() + ";") ""
                            } 
                            |> Option.Some
                        with
                        | ex -> 
                            printfn "%A" ex
                            Option.None
                       )
        |> Array.filter Option.isSome
        |> Array.map (fun x -> x.Value)
        |> FSharpAux.IO.SeqIO.Seq.toCSV "\t" true
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)