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
open PeakQuantification
open PSMQuantification

module PSMBasedQuantification =
    module Query = 
        open MzIO.Commons.Arrays
        open MzIO.Processing.MzIOLinq
         
        /// The retention time index and the scan times of the MS1 spectra, both sorted by scan time, from one pass
        /// over the spectrum descriptions with the streaming parser (SpectrumIndex.readMs1Index).
        let readMs1Index (reader: MzIO.MzSQL.MzSQL) (runID: string) (warn: string -> unit) =
            let ms1SortedByScanTime = SpectrumIndex.readMs1Index reader runID warn ignore
            let rtIndex =
                ms1SortedByScanTime
                |> Array.map (fun (scanTime, id) -> RtIndexEntry(scanTime, id))
                |> MzIOArray.ToMzIOArray
            rtIndex, ms1SortedByScanTime


    ///
    let lightQualityFilter lowerBorder upperBorder (quantResults:QuantificationResult[]) =
        let medianApexIntensities = 
            quantResults
            |> Array.map (fun x -> x.MeasuredApex_Light)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        let medianQuantIntensities = 
            quantResults
            |> Array.map (fun x -> x.Quant_Light)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        quantResults
        |> Array.filter (fun x -> 
            let qualR = 
                let apexNorm = x.MeasuredApex_Light / medianApexIntensities
                let quantNorm =  x.Quant_Light / medianQuantIntensities
                quantNorm / apexNorm
                |> FSharp.Stats.Ops.log2
            (qualR > lowerBorder && qualR < upperBorder) || (nan.Equals(qualR) )
            ) 

    ///
    let heavyQualityFilter lowerBorder upperBorder (quantResults:QuantificationResult[]) =
        let medianApexIntensities = 
            quantResults
            |> Array.map (fun x -> x.MeasuredApex_Heavy)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        let medianQuantIntensities = 
            quantResults
            |> Array.map (fun x -> x.Quant_Heavy)
            |> Array.filter (fun x -> nan.Equals x |> not)
            |> Array.median
        quantResults
        |> Array.filter (fun x -> 
            let qualR = 
                let apexNorm = x.MeasuredApex_Heavy / medianApexIntensities
                let quantNorm =  x.Quant_Heavy / medianQuantIntensities
                quantNorm / apexNorm
                |> FSharp.Stats.Ops.log2
            (qualR > lowerBorder && qualR < upperBorder) || (nan.Equals(qualR) )
            )
    
    ///
    let initGetProcessedXIC (baseLineCorrection:Domain.BaseLineCorrection option) readSortedPeaks idx scanTimeWindow mzWindow_Da meanScanTime meanPrecMz =
        let rtQuery = Query.createRangeQuery meanScanTime scanTimeWindow
        let rt, intensity = Xic.rtProfile readSortedPeaks idx rtQuery meanPrecMz mzWindow_Da nan nan
        Xic.processXic baseLineCorrection rt intensity

    ///
    let quantifyPeptides diagCharts zipCharts (processParams:Domain.QuantificationParams) (outputDir:string) (cn:SQLiteConnection) (instrumentOutput:string) (scoredPSMs:string)  =
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
        logger.Trace "Read scored PSMs"
        ///
        let qpsms =
            Csv.CsvReader<PSMStatisticsResult>(SchemaMode=Csv.Fill).ReadFile(scoredPSMs,'\t',false,1)
            |> Array.ofSeq
        logger.Trace "Read scored PSMs:finished"
        if Array.isEmpty qpsms then
            logger.Warn (sprintf "%s contains no PSMs, %s gets no quantifications." scoredPSMs outFilePath)
            Seq.empty<QuantificationResult>
            |> SeqIO'.csv "\t" true false
            |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        else
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
        let dBParams     = getSDBParamsByCn memoryDB
        //let massLookUp = prepareSelectMassByModSequenceAndGlobalMod memoryDB
        let peptideLookUp = getThreadSafePeptideLookUpFromFileBySequenceAndGMod memoryDB dBParams
        let isotopePatterns = new IsotopePatterns.IsotopePatternLookUp(logger, memoryDB)
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
                |> List.collect (fun ch -> 
                    ionSeries 
                    |> List.map (fun x -> Mass.toMZ x.MainPeak.Mass ch)
                )                  
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
        logger.Trace "Create RetentionTime index and read the ms1s"
        /// retTimeIdxed and ms1SortedByScanTime: the MS1s sorted by their scanTime
        let retTimeIdxed, ms1SortedByScanTime = Query.readMs1Index inReader inRunID (fun message -> logger.Warn message)
        logger.Trace (sprintf "Create RetentionTime index and read the ms1s:finished, %i MS1s" ms1SortedByScanTime.Length)
        
        let calibration = calibratePsms logger diagCharts plotDirectory qpsms

        logger.Trace "init lookup functions"        
        // The ions are quantified in the order of their scan time and their XICs read the MS1 spectra of a window of
        // 2 x ScanTimeWindow, so the spectra of twice the widest window are kept instead of all MS1 spectra of the file.
        let peakCache =
            let cacheSize, widestWindow = Xic.peakCacheSize (ms1SortedByScanTime |> Array.map fst) processParams.XicExtraction.ScanTimeWindow
            logger.Trace (sprintf "MS1 peaks cache: %i spectra (%i MS1s in the widest XIC window)" cacheSize widestWindow)
            Xic.PeakCache(Xic.readSortedPeaks inReader, cacheSize)
        let readSortedPeaksWithMem = peakCache.Get
        let comparePredictedAndMeasuredIsotopicCluster = initComparePredictedAndMeasuredIsotopicCluster isotopePatterns (fun ms1ID lowMz highMz -> Xic.getPeaksInMzRange (readSortedPeaksWithMem ms1ID) lowMz highMz) (SpectrumIndex.initClosestMs1 ms1SortedByScanTime) (calibration.HalfWindow / MzCalibration.mzWindowSdFactor)
        
        let mzWindow = MzCalibration.mzWindow processParams.XicExtraction.MzWindow_Da calibration
        logger.Trace (sprintf "XIC %s" (MzCalibration.describeWindow processParams.XicExtraction.MzWindow_Da mzWindow (calibration.HalfWindow = MzCalibration.fallbackMzWindow)))

        ///
        let getXIC = initGetProcessedXIC processParams.BaseLineCorrection readSortedPeaksWithMem retTimeIdxed processParams.XicExtraction.ScanTimeWindow mzWindow
        // PSMBasedQuantification has no ion mobility.
        let getXICWithIM scanTime mz (_: float) = getXIC scanTime mz
        
        match processParams.XicExtraction.XicProcessing with
        | Domain.XicProcessing.Gabor3D _ -> failwith "Gabor3D peak detection is only available in PSMBasedQuantificationTIMs."
        | _ -> ()
        ///
        let identifyPeaks = initIdentifyPeaks processParams.XicExtraction.XicProcessing
        logger.Trace "init lookup functions:finished"
        
        logger.Trace "init quantification functions"
        let quantificationFunctions : QuantificationFunctions =
            {
                Logger = logger
                DiagCharts = diagCharts
                PlotDirectory = plotDirectory
                ScanTimeWindow = processParams.XicExtraction.ScanTimeWindow
                CentrePSMs = processParams.CentrePSMs
                PeptideLookUp = peptideLookUp
                CountMatchedMasses = countMatchedMasses
                CorrectionAt = calibration.CorrectionAt
                IonMobilityOf = (fun _ -> nan)
                GetXic = getXICWithIM
                IdentifyPeaks = identifyPeaks
                ComparePredictedAndMeasuredIsotopicCluster = comparePredictedAndMeasuredIsotopicCluster
            }
        logger.Trace "init quantification functions:finished"
        
        logger.Trace "executing quantification"
        let quantResults = 
            qpsms 
            |> Array.groupBy (fun x -> 
                {
                    Sequence             = x.StringSequence     
                    GlobalMod            = x.GlobalMod    
                    Charge               = x.Charge       
                    ModSequenceID        = x.ModSequenceID
                    PepSequenceID        = x.PepSequenceID
                }
                )
            |> Array.map (fun (pepIon,psms) -> 
                    match processParams.XicExtraction.TopKPSMs with 
                    | Some x when psms.Length > x -> 
                        pepIon,
                        psms
                        |> Array.sortByDescending (fun x -> x.ModelScore)
                        |> Array.take x
                    | _ -> 
                        pepIon,
                        psms                  
                )
            |> Xic.mapiOrderedBy (fun (_, psms: PSMStatisticsResult []) -> if psms.Length = 0 then infinity else psms |> Array.averageBy (fun x -> x.ScanTime)) (fun i (pepIon,psms) -> 
                if i % 100 = 0 then logger.Trace (sprintf "%i peptides quantified" i)
                
                match processParams.PerformLabeledQuantification with 
                |Domain.Labeling.N15Labeling | Domain.Labeling.N15LabelingOnly | Domain.Labeling.Labelshift-> quantifyLabeled quantificationFunctions pepIon psms
                |Domain.Labeling.Unlabeled -> quantifyLabelFree quantificationFunctions pepIon psms
                )
            |> Array.choose id
                    
        let limit = processParams.QualityFilterLog2Limit
        let filteredResults = 
            match processParams.PerformLabeledQuantification with 
            |Domain.Labeling.N15Labeling ->
                quantResults
                |> heavyQualityFilter (-limit) limit
                |> lightQualityFilter (-limit) limit
            |Domain.Labeling.N15LabelingOnly ->
                quantResults
                |> heavyQualityFilter (-limit) limit
            |Domain.Labeling.Unlabeled ->
                quantResults
                |> lightQualityFilter (-limit) limit
            |Domain.Labeling.Labelshift ->
                quantResults
        filteredResults 
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        if zipCharts then
            plotDirectory
            |> Zipping.zipDirectory "*.html" logger
            |> fun zipped ->
                match zipped with
                | Ok byteArr ->
                    byteArr
                    |> Zipping.saveZippedDirectory outputDir logger ((Path.GetFileNameWithoutExtension instrumentOutput) + "_plots")
                    |> fun saved ->
                        match saved with
                        | Ok save -> 
                            save
                            Directory.Delete (plotDirectory, true)
                        | Error ex -> logger.Trace (sprintf "Error saving zipped directory: %A" ex)
                | Error ex -> logger.Trace (sprintf "Error zipping directory: %A" ex)
        logger.Trace (sprintf "Isotope patterns: %i read from the peptide data base, %i computed" isotopePatterns.Read isotopePatterns.Computed)
        if peakCache.Reads > 0 then logger.Trace (sprintf "%i MS1 spectra were read for %i MS1s of the file (cache of %i spectra)" peakCache.Reads ms1SortedByScanTime.Length peakCache.Capacity)
        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()
        (isotopePatterns :> System.IDisposable).Dispose()
        memoryDB.Dispose()
        logger.Trace "executing quantification:finished"
        