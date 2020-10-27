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

module ConsensusSpectralLibrary =
    
//let getBinIdx' width offset x = int ((x / width) + offset)

//let minMassBoarder =  getBinIdx' 0.01 0.0 100.
//let maxMassBoarder =  getBinIdx' 0.01 0.0 1501.
//let maxIndex = maxMassBoarder - minMassBoarder + 1
//let index = (getBinIdx' 0.01 0.0 1501.0) - minMassBoarder
    type SparsePeakArray = { 
        Length: int
        Data: System.Collections.Generic.IDictionary<int,float>
        }
    
    let dot (x:SparsePeakArray) (y:SparsePeakArray) =
        x.Data
        |> Seq.fold (fun (acc:float) xi -> 
            let present,yi = y.Data.TryGetValue xi.Key
            if present then acc + (yi * xi.Value) else acc
            ) 0.
            
    let getBinIdx' width offset x = int ((x / width) + offset)

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


    ///// Bins peaks to their nearest 1 Da bin. Filters out peaks where the mz < minMassBoarder & > maxMassBoarder
    //let peaksToNearestBinVector binWidth offset (pkarr:BioFSharp.Mz.PeakArray<_>) (minMassBoarder:float) (maxMassBoarder:float) =
    //    let minMassBoarder = getBinIdx' binWidth offset minMassBoarder
    //    let maxMassBoarder = getBinIdx' binWidth offset maxMassBoarder
    //    let maxIndex = maxMassBoarder - minMassBoarder + 1        
    //    let vector = Vector.create (maxIndex-1) 0.
    //    pkarr 
    //    |> Array.iter (fun p ->  
    //        let index = (getBinIdx' binWidth offset p.Mz) - minMassBoarder
    //        if index < maxIndex-1 && index > -1 then
    //            vector.[index] <- max vector.[index] p.Intensity
    //            )
    //    vector

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

    let getClosestMz (peaks: Peak1DArray, mzRange: RangeQuery) =
        let p1d = 
            peaks.Peaks
                .DefaultIfEmpty(new Peak1D(0., mzRange.LockValue))
                .ItemAtMin(fun x -> Math.Abs(x.Mz - mzRange.LockValue))
        if p1d.Mz > mzRange.HighValue || p1d.Mz < mzRange.LowValue then
            new Peak1D(0., mzRange.LockValue)
        else
            p1d

    let getRTProfiles (swathIndexer: SwathIndexer.SwathIndexer) (reader: IMzIODataReader) (*(spectrumSelector: seq<SwathIndexer.MSSwath> -> seq<SwathIndexer.MSSwath> list)*) (swathQuery: SwathQuery) =
        swathIndexer.GetRTProfilesFirstWnd(reader, swathQuery, false, (fun x -> [x.Take(1)]) ,getClosestMz)

    ///
    let readLibraryFrom path= 
        let b: ProteomIQon.Dto.PeptideIon[] = ProteomIQon.Json.ReadAndDeserialize path
        b

    ///
    let createSpectralLibrary (logger:NLog.Logger) (processParams:Domain.ConsensusSpectralLibraryParams) (outputDir:string) (libraryFiles:string) (targetSwathFile:string) = 
        logger.Trace (sprintf "Input directory containing library files: %s" libraryFiles)
        logger.Trace (sprintf "Input Swath file to align to: %s" targetSwathFile)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        
        let getPlotFilePathFilePath (plotName:string) (fileName:string) =
            let fileName = (Path.GetFileNameWithoutExtension fileName) + "_" + plotName 
            Path.Combine [|outputDir;fileName|]
        
        logger.Trace "Init align function"
        let ctx = new ML.MLContext()
        let rnd = new System.Random()
        let align = (*initAlign *) ctx
        logger.Trace "Init align function: finished"
         
        logger.Trace "Reading and preparing .sl files for alignment"
        let libraryFilePaths = System.IO.Directory.GetFiles (libraryFiles, "*.sl")
        let libraryFiles = libraryFilePaths |> Array.map readLibraryFrom 
        logger.Trace "Reading and preparing .sl files for alignment: finished"
        
        logger.Trace "Init connection to swath mass spectrum data."
        let inReader = Core.MzIO.Reader.getReader targetSwathFile
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader
        let inTr = inReader.BeginTransaction()
        logger.Trace "Init connection to swath mass spectrum data.:finished"

        logger.Trace "Create Swath RetentionTime index"
        let retTimeIdxed = Query.getSwathIdx inReader inRunID
        logger.Trace "Create Swath RetentionTime index:finished"

        logger.Trace "Identify features in Swath file"
        let getBinIdx width scantime = int ((scantime / width))                          
        let swathFileFeatures = 
            libraryFiles
            |> Array.map (fun library ->
                    let libBinned = 
                        library
                        |> Array.map (fun pepIons -> {pepIons with Fragments = pepIons.Fragments |> List.filter (fun x -> x.Number > processParams.MinFragmentLadderIdx) })
                        |> Array.filter (fun pepIons -> pepIons.Fragments.Length > processParams.MinFragmentCount)
                        |> Array.filter (fun pepIons -> pepIons.StringSequence |> String.filter Char.IsUpper |> String.length >= processParams.MinPeptideLength)
                        |> Array.groupBy (fun ion -> getBinIdx processParams.BinningWindowWidth ion.ScanTime)
                        |> Array.sortBy fst
                    Chart.Column(libBinned |> Array.map (fun (x,y) -> x,y.Length))
                    |> Chart.Show
                    let filtered =
                        libBinned 
                        |> Array.map (fun (binIdx,pepIons) -> 
                            let itemsSorted = pepIons |> Array.sortByDescending (fun x -> x.MeasuredApex)
                            let frac = ((itemsSorted.Length |> float) * processParams.FractionOfMostAbundandIonsPerBin) |> int
                            binIdx,
                            itemsSorted
                            |> Array.take (Math.Min(itemsSorted.Length-1,frac))
                            |> Array.sortBy (fun x -> x.ScanTime)
                            )                   
                    //[
                    Chart.Column(libBinned |> Array.map (fun (x,y) -> x,y.Length))
                    |> Chart.withTraceName "libBinned"
                    |> Chart.Show
                    Chart.Column(filtered |> Array.map (fun (x,y) -> x,y.Length))
                    |> Chart.withTraceName (sprintf "libBinned_filtered: top %f" processParams.FractionOfMostAbundandIonsPerBin)
                    |> Chart.Show
                    //]
                    //|> Chart.Combine
                    //|> Chart.Show
                    [
                    Chart.Column(libBinned |> Array.map (fun (x,y) -> x,y.Length))
                    |> Chart.withTraceName "libBinned"
                    Chart.Column(filtered |> Array.map (fun (x,y) -> x,y.Length))
                    |> Chart.withTraceName (sprintf "libBinned_filtered: top %f" processParams.FractionOfMostAbundandIonsPerBin)
                    
                    ]
                    |> Chart.Combine
                    |> Chart.Show
                    let estimatedScanTimesInSwathFile = 
                        filtered
                        |> Array.map snd
                        |> Array.concat
                        |> Array.choose (fun pepIon -> 
                                let fragmentMzs = 
                                    pepIon.Fragments
                                    |> List.map (fun f -> f.MeanFragMz,f.MeanRelativeIntensity_Frags)
                                    |> Array.ofSeq
                                    |> BioFSharp.Mz.PeakArray.zipMzInt
                                    
                                let fragmentVector = 
                                    fragmentMzs 
                                    //|> fun x -> BioFSharp.Mz.PeakArray.peaksToNearestUnitDaltonBinVector x 100. 2000.
                                    |> fun (x) -> peaksToNearestBinVector processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset x (processParams.MS2ScanRange |> fst) (processParams.MS2ScanRange |> snd)
                                let rtQuery   = MzIO.Processing.Query.createRangeQuery pepIon.ScanTime processParams.RtWindowWidth
                                let mzQueries = 
                                    fragmentMzs 
                                    |> Array.map (fun pk -> MzIO.Processing.Query.createRangeQuery pk.Mz 0.05)
                                
                                let swathQuery = MzIO.Processing.Query.createSwathQuery pepIon.PrecursorMZ rtQuery mzQueries
                                let profiles   = getRTProfiles retTimeIdxed inReader swathQuery
                                let correlationTrace = 
                                    match profiles with 
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
                                                //|> fun x -> BioFSharp.Mz.PeakArray.peaksToNearestUnitDaltonBinVector x 100. 2000.
                                                |> fun (x) -> peaksToNearestBinVector processParams.FragMatchingBinWidth processParams.FragMatchingBinOffset x (processParams.MS2ScanRange |> fst) (processParams.MS2ScanRange |> snd)
                                
                                                )
                                            |> Array.sortBy fst
                                        let correlationTrace = 
                                            fragVecs
                                            |> Array.map (fun (rt, measuredFrags) -> 
                                                rt, dot fragmentVector measuredFrags
                                                )
                                        [
                                        Chart.Point correlationTrace
                                        Chart.Point [pepIon.ScanTime,correlationTrace |> Array.maxBy snd |> snd]
                                        ]
                                        |> Chart.Combine
                                        |> Chart.Show
                                        Some correlationTrace
                                    | _ -> None 
                                correlationTrace
                            )
                    estimatedScanTimesInSwathFile
                )

            
        logger.Trace "Determining Alignment order: finished"
        
        logger.Trace "Plotting file distances"
        //let chart = 
        //    alignmentFilesOrdered
        //    |> Array.map (fun (target,sources) ->
        //            Chart.Point(sources |> Array.mapi (fun i x -> (snd x).FileName, fst x))
        //            |> Chart.withTraceName target.FileName
        //            |> Chart.withX_AxisStyle("FileNames")
        //            |> Chart.withY_AxisStyle("Median absolute difference of peptide ion scan times")
        //            |> Chart.withSize(1000.,1000.)
        //            |> Chart.SaveHtmlAs(getPlotFilePathFilePath "differences" target.FileName)
        //        )
        logger.Trace "Plotting file distances: finished"