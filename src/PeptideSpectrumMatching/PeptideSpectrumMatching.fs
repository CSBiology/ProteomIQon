namespace ProteomIQon

open System
open BioFSharp
open BioFSharp.Mz.SearchDB
open Domain
open Core
open System.IO
open BioFSharp.Mz
open MzIO
open MzIO.IO
open MzIO.Processing
open FSharpAux.IO
module PeptideSpectrumMatching =

    let tryGet tolerance (xTar:float) (data:#seq<(float*float)>) = 
        let xTmp,yTmp = data |> Seq.minBy (fun (x,y) -> abs(x-xTar))
        if abs(xTmp-xTar) < tolerance then Some(xTmp,yTmp) else None

    open System.IO
    open System.Data
    open BioFSharp.Mz.TheoreticalSpectra
    open BioFSharp.Mz.ChargeState

    let getPrecursorCharge (chParams:ChargeState.ChargeDetermParams) rnd inRunID (inReader: IMzIODataReader) =
        /// Returns a Sequence containing all MassSpectra of a single MS-run
        let massSpectra = inReader.ReadMassSpectra(inRunID)
        /// Returns a Array that contains all MS1s sorted by their scanTime
        let ms1SortedByScanTime =
            massSpectra
            |> Seq.filter (fun ms -> MassSpectrum.getMsLevel ms = 1)
            |> Seq.sortBy  MassSpectrum.getScanTime
            |> Array.ofSeq
        let ms2SortedByScanTime =
            massSpectra
            |> Seq.filter (fun ms -> MassSpectrum.getMsLevel ms = 2)
            |> Seq.sortBy  MassSpectrum.getScanTime
            |> Array.ofSeq
        let ms2AssignedToMS1 =
            ms2SortedByScanTime
            |> Array.map (fun ms2 ->
                            ms1SortedByScanTime
                            |> Array.findBack (fun ms1 ->
                                                    let ms1ScanTime = MassSpectrum.getScanTime ms1
                                                    let ms2ScanTime = MassSpectrum.getScanTime ms2
                                                    ms1ScanTime <= ms2ScanTime),ms2
                          )
            |> Array.groupBy fst
            |> Array.map (fun (ms1,ms1And2) -> ms1,ms1And2 |> Array.map snd)
        /// Returns a List of the type ChargeDetInterimResult. Each Item contains a ms1*ms2 tuple and a List of putative Chargestates
        let ms2PossibleChargestates =
            ms2AssignedToMS1
            |> Array.mapi (fun i (ms1,ms2s) ->
                            //printfn "%i" i
                            let mzdata,intensityData = Peaks.unzipIMzliteArray (inReader.ReadSpectrumPeaks(ms1.ID).Peaks)
                            ms2s
                            |> Array.filter (fun ms2 -> inReader.ReadSpectrumPeaks(ms2.ID).Peaks |> Seq.isEmpty = false)
                            |> Array.map (fun ms2 ->
                                                let assignedCharges = ChargeState.putativePrecursorChargeStatesBy chParams mzdata intensityData ms1.ID ms2.ID (MassSpectrum.getPrecursorMZ ms2)
                                                match assignedCharges with
                                                | [] ->
                                                    [
                                                        for i = 2 to 3 do
                                                            let precursorMz = (MassSpectrum.getPrecursorMZ ms2)
                                                            let mass = Mass.ofMZ precursorMz (float i)
                                                            let score = getScore 10 1 100.
                                                            yield createAssignedCharge ms1.ID ms2.ID precursorMz i mass 100. score [0.] 0 0. (Set[]) (Some 1.)
                                                    ]
                                                | _ -> assignedCharges
                                            )
                            )
            |> Array.filter (fun x ->  Array.isEmpty x |> not)
            |> Array.concat
            |> Array.filter (fun (assignedCharges)  -> assignedCharges <> [])
            |> List.ofArray
        /// Returns the standard deviation of the peak positions determined by using the best scored subsets of each attempt of mapping a ms2 to a chargestate
        let peakPosStdDev =
            ms2PossibleChargestates
            |> List.filter (fun assignedCharges -> assignedCharges.Head.PositionMetricPValue.IsNone)
            |> List.map (fun (assignedCharges) -> assignedCharges.Head )
            |> ChargeState.peakPosStdDevBy
        /// Returns a function that generates a usergiven amount of rnd spectra and calculates their mzDeviation. The
        /// input needed for this function is the nrOfPeaksInASubset and the putative chargestate
        let init = ChargeState.initMzDevOfRndSpec rnd {chParams with ExpectedMaximumCharge=8} peakPosStdDev//peakPosStdDev
        /// Returns a List of the type ChargeDetInterimResult. Each Item contains a ms1*ms2 tuple and a List of putative chargestates.
        /// Each putative chargestate is hypothesis tested and posesses a PValue <= 0.05
        let positionMetricScoredCharges =
            ms2PossibleChargestates
            |> List.mapi (fun i (assignedCharges) ->
                                let ms1ID = assignedCharges.Head.PrecursorSpecID
                                let ms2ID = assignedCharges.Head.ProductSpecID
                                let items =
                                    assignedCharges
                                    |> List.mapi (fun  i putativeCharge ->
                                                    match putativeCharge.PositionMetricPValue with
                                                    | Some x ->
                                                        putativeCharge
                                                    | None   ->
                                                        let pValue = ChargeState.empiricalPValueOfSim init (putativeCharge.SubSetLength,float putativeCharge.PrecCharge ) putativeCharge.MZChargeDev
                                                        {putativeCharge with PositionMetricPValue = Some pValue}
                                                 )
                                    |> List.filter (fun testIt -> testIt.PositionMetricPValue.Value <= 0.05)
                                    |> ChargeState.removeSubSetsOfBestHit
                                    |> (fun charges ->
                                            match charges with
                                            | [] -> [
                                                        for i = 2 to 3 do
                                                            let precMz = MassSpectrum.getPrecursorMZ (inReader.ReadMassSpectrum ms2ID)
                                                            let mass = Mass.ofMZ precMz (float i)
                                                            let score = ChargeState.getScore 10 1 100.
                                                            yield ChargeState.createAssignedCharge ms1ID ms2ID precMz i mass 100. score [0.] 0 0. (Set[]) None
                                                    ]
                                            | _  -> charges
                                        )
                                items
                            )

        positionMetricScoredCharges,peakPosStdDev

    let psm (processParams:PeptideSpectrumMatchingParams) peakPosStdDev lookUpF calcIonSeries (reader: IMzIODataReader) (outFilePath: string) (ms2sAndAssignedCharges: AssignedCharge list list) =
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension outFilePath)
        let resultWriter = new System.IO.StreamWriter(outFilePath, true)
        let (ms2IDAssignedCharge) =
                ms2sAndAssignedCharges
                |> List.mapi (fun i (assignedCharges) ->
                                assignedCharges |> List.map (fun assCh -> i,assCh.ProductSpecID, assCh)
                            )
                |> List.concat
                |> List.groupBy (fun (ascendingID,ms2Id,assCH) -> assCH.PrecCharge)
                |> List.map (fun (ch,ms2IDassCHL) ->
                                ch, ms2IDassCHL |> List.sortBy (fun (ascendingID,ms2ID,assCh) -> assCh.PutMass)
                            )
                |> List.sortBy fst

        ms2IDAssignedCharge
        |> List.iter (fun (ch, ms2IdAssCh) -> logger.Trace (sprintf "%i spectra with charge %i" ms2IdAssCh.Length ch))

        ms2IDAssignedCharge
        |> List.iter (fun (ch,ms2IdAssCh) ->
                        logger.Trace (sprintf "%i spectra with charge %i processed" ms2IdAssCh.Length ch)
                        ms2IdAssCh
                        |> List.iteri (fun i (ascendingID,ms2Id,assCh) ->
                                            if i%10000 = 0 then logger.Trace (sprintf  "%i" i)
                                            try
                                            let scanTime = MassSpectrum.getScanTime (reader.ReadMassSpectrum(ms2Id))
                                            let recSpec =
                                                Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(ms2Id).Peaks)
                                                |> fun (mzData,intensityData) -> PeakArray.zip mzData intensityData
                                            let floorToClosest10 x = 
                                                Math.Floor(x / 10.) * 10.            
                                            let ceilToClosest10 x =
                                                Math.Ceiling(x / 10.) * 10.
                                            let scanRange = //processParams.MS2ScanRange
                                                let low = Math.Max(0.,System.Math.Round((recSpec |> Array.minBy (fun x -> x.Mz)).Mz,0) |> floorToClosest10)
                                                let top = Math.Round((recSpec |> Array.maxBy (fun x -> x.Mz)).Mz,0)  |> ceilToClosest10
                                                low,top

                                            let getLookUpsResults precMz = 
                                                let lowerMass,upperMass =
                                                    let massWithH2O = BioFSharp.Mass.ofMZ precMz (ch|> float) 
                                                    Mass.rangePpm processParams.LookUpPPM massWithH2O
                                                let lookUpResults :SearchDB.LookUpResult<AminoAcids.AminoAcid> list =
                                                    lookUpF lowerMass upperMass
                                                lookUpResults
                                            let precMz' = 
                                                let ms1Spec = 
                                                    Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(assCh.PrecursorSpecID).Peaks)
                                                    ||> Seq.zip
                                                let mz = assCh.PrecursorMZ
                                                let monoIso = tryGet (3. * peakPosStdDev) mz ms1Spec
                                                let mzMinusOne = mz - (Mass.Table.PMassInU / (float ch)) 
                                                let minusOne = tryGet (3. * peakPosStdDev) mzMinusOne ms1Spec
                                                match monoIso,minusOne with
                                                //| Some (x,picked),Some (x2,minusOne) -> 
                                                //    if (picked / minusOne) > 0.8 then Some x2 else None
                                                //| Some (x,picked),None -> 
                                                //    Some mzMinusOne
                                                
                                                | _ -> 
                                                    //None
                                                    Some mzMinusOne
                                            let lookUpResults = 
                                                match precMz' with 
                                                | Some pMz -> 
                                                    (getLookUpsResults pMz)@(getLookUpsResults assCh.PrecursorMZ)
                                                    |> List.distinctBy (fun x -> x.ModSequenceID)
                                                | None -> 
                                                    (getLookUpsResults assCh.PrecursorMZ)
                                                
                                            let theoSpecsN =
                                                lookUpResults
                                                |> List.map (fun lookUpResult ->
                                                                let ionSeries = calcIonSeries lookUpResult.BioSequence
                                                                lookUpResult,ionSeries
                                                            )
                                       
                                            let theoSpecs = theoSpecsN(*@theoSpecsM1@theoSpecsP1*)

                                            let sequestTheoreticalSpecs = SequestLike.getTheoSpecs scanRange assCh.PrecCharge theoSpecs
                                            let sequestLikeScored =
                                                try
                                                SequestLike.calcSequestScore scanRange recSpec scanTime assCh.PrecCharge
                                                    assCh.PrecursorMZ sequestTheoreticalSpecs ms2Id
                                                with 
                                                | _ -> 
                                                    let recSpec =
                                                        Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(ms2Id).Peaks)
                                                        |> fun (mzData,intensityData) -> 
                                                            Array.zip mzData intensityData
                                                    logger.Trace (sprintf "spec with id: %s MinScan: %f MaxScan: %f , theoSpecLength: %A fails with: %A " ms2Id (fst scanRange) (snd scanRange) (sequestTheoreticalSpecs |> List.map (fun x -> x.DecoyTheoSpec.Length,x.TheoSpec.Length )) recSpec ); []
                                            
                                            
                                            let bestTargetSequest =
                                                sequestLikeScored
                                                |> List.filter (fun (x:SearchEngineResult.SearchEngineResult<float>) -> x.IsTarget)
                                                |> List.truncate 10
                                                |> List.map (fun x -> (x.ModSequenceID,x.GlobalMod),x)
                                                |> Map.ofList

                                            let bestDecoySequest =
                                                sequestLikeScored
                                                |> List.filter (fun x -> not x.IsTarget)
                                                |> List.truncate 10
                                                |> List.map (fun x -> (x.ModSequenceID,x.GlobalMod),x)
                                                |> Map.ofList

                                            let andromedaTheorticalSpecs =
                                                theoSpecs
                                                |> List.filter (fun (lookUpResult,fragments) ->
                                                                    bestTargetSequest |> Map.containsKey (lookUpResult.ModSequenceID,lookUpResult.GlobalMod) ||
                                                                    bestDecoySequest  |> Map.containsKey (lookUpResult.ModSequenceID,lookUpResult.GlobalMod)
                                                               )

                                                |> XScoring.getTheoSpecs scanRange assCh.PrecCharge    
                                                
                                            let andromedaLikeScored,xtandemScored = 
                                                XScoring.calcAndromedaAndXTandemScore processParams.AndromedaParams.PMinPMax scanRange processParams.AndromedaParams.MatchingIonTolerancePPM
                                                    recSpec scanTime assCh.PrecCharge assCh.PrecursorMZ andromedaTheorticalSpecs ms2Id

                                            let result = 
                                                List.map2 (fun (androRes:SearchEngineResult.SearchEngineResult<float>) (xTandemRes:SearchEngineResult.SearchEngineResult<float>) -> 
                                                    // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
                                                    //let pSMId = androRes.SpectrumID.Replace(' ', '-') + "_" + ascendingID.ToString() + "_" + ch.ToString() + "_" + i.ToString()
                                                    let label = if androRes.IsTarget then 1 else -1
                                                    let scanNr = ascendingID
                                                    let absDeltaMass = (androRes.TheoMass-androRes.MeasuredMass) |> abs
                                                    match androRes.IsTarget with
                                                    | true ->
                                                        match Map.tryFind (androRes.ModSequenceID,androRes.GlobalMod) bestTargetSequest with
                                                        | Some sequestScore ->
                                                            let res i : Dto.PeptideSpectrumMatchingResult = 
                                                                let pSMId i = androRes.SpectrumID.Replace(' ', '-') + "_" + ascendingID.ToString() + "_" + ch.ToString() + "_" + i.ToString()
                                                                {                                                                        
                                                                    PSMId                        = pSMId i
                                                                    GlobalMod                    = androRes.GlobalMod
                                                                    PepSequenceID                = androRes.PepSequenceID
                                                                    ModSequenceID                = androRes.ModSequenceID
                                                                    Label                        = label
                                                                    ScanNr                       = scanNr
                                                                    ScanTime                     = androRes.ScanTime
                                                                    Charge                       = ch
                                                                    PrecursorMZ                  = androRes.PrecursorMZ
                                                                    TheoMass                     = androRes.TheoMass
                                                                    AbsDeltaMass                 = absDeltaMass
                                                                    PeptideLength                = androRes.PeptideLength
                                                                    MissCleavages                = -1
                                                                    SequestScore                 = sequestScore.Score
                                                                    SequestNormDeltaBestToRest   = sequestScore.NormDeltaBestToRest
                                                                    SequestNormDeltaNext         = sequestScore.NormDeltaNext
                                                                    AndroScore                   = androRes.Score
                                                                    AndroNormDeltaBestToRest     = androRes.NormDeltaBestToRest
                                                                    AndroNormDeltaNext           = androRes.NormDeltaNext
                                                                    XtandemScore                 = xTandemRes.Score  
                                                                    XtandemNormDeltaBestToRest   = xTandemRes.NormDeltaBestToRest  
                                                                    XtandemNormDeltaNext         = xTandemRes.NormDeltaNext  
                                                                    StringSequence               = androRes.StringSequence
                                                                    ProteinNames                 = "PlaceHolder"
                                                                }
                                                            Some (label,res)
                                                        | None -> None 
                                                    | false ->
                                                        match Map.tryFind (androRes.ModSequenceID,androRes.GlobalMod) bestDecoySequest with
                                                        | Some sequestScore ->
                                                            let res i : Dto.PeptideSpectrumMatchingResult = 
                                                                let pSMId i = androRes.SpectrumID.Replace(' ', '-') + "_" + ascendingID.ToString() + "_" + ch.ToString() + "_" + i.ToString()
                                                                            
                                                                {                                                                        
                                                                    PSMId                        = pSMId i
                                                                    GlobalMod                    = androRes.GlobalMod
                                                                    PepSequenceID                = androRes.PepSequenceID
                                                                    ModSequenceID                = androRes.ModSequenceID
                                                                    Label                        = label
                                                                    ScanNr                       = scanNr
                                                                    ScanTime                     = androRes.ScanTime
                                                                    Charge                       = ch
                                                                    PrecursorMZ                  = androRes.PrecursorMZ
                                                                    TheoMass                     = androRes.TheoMass
                                                                    AbsDeltaMass                 = absDeltaMass
                                                                    PeptideLength                = androRes.PeptideLength
                                                                    MissCleavages                = -1
                                                                    SequestScore                 = sequestScore.Score
                                                                    SequestNormDeltaBestToRest   = sequestScore.NormDeltaBestToRest
                                                                    SequestNormDeltaNext         = sequestScore.NormDeltaNext
                                                                    AndroScore                   = androRes.Score
                                                                    AndroNormDeltaBestToRest     = androRes.NormDeltaBestToRest
                                                                    AndroNormDeltaNext           = androRes.NormDeltaNext
                                                                    XtandemScore                 = xTandemRes.Score  
                                                                    XtandemNormDeltaBestToRest   = xTandemRes.NormDeltaBestToRest  
                                                                    XtandemNormDeltaNext         = xTandemRes.NormDeltaNext 
                                                                    StringSequence               = androRes.StringSequence
                                                                    ProteinNames                 = "PlaceHolder"
                                                                }
                                                            Some (label,res)
                                                        | None -> None 
                                                ) andromedaLikeScored xtandemScored
                                                |> List.choose id
                                                |> List.groupBy (fun x -> fst x)
                                                |> List.sortBy fst
                                                |> List.map (fun (x,y) ->
                                                                y
                                                                |> List.map snd
                                                                |> List.mapi (fun i x -> x i)
                                                            )
                                                |> List.concat

                                            result
                                            |> SeqIO.Seq.toCSV "\t" false
                                            |> Seq.iter resultWriter.WriteLine

                                            with

                                            | _ as ex ->
                                                    logger.Trace (sprintf "spec with id: %s fails with: %A" ms2Id ex )
                                        )

                        )
        resultWriter.Dispose()
    //    Logger.printTimenWithPre pre "Finished PSM"

    let scoreSpectra (processParams:PeptideSpectrumMatchingParams) (outputDir:string) (d:string)  (instrumentOutput:string) =

        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension instrumentOutput)

        logger.Trace (sprintf "Input file: %s" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)

        // initialize Reader and Transaction
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".psm"
            Path.Combine [|outputDir;fileName|]
        logger.Trace (sprintf "Result file path: %s" outFilePath)

        logger.Trace "Copy peptide DB into Memory."
        let cn =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        let memoryDB = SearchDB.copyDBIntoMemory cn
        let memoryDBTr = memoryDB.BeginTransaction()
        cn.Dispose()
        logger.Trace "Copy peptide DB into Memory: finished."

        logger.Trace "Prepare processing functions."

        let chargeParams = processParams.ChargeStateDeterminationParams
        logger.Trace (sprintf "Charge parameters: %A" chargeParams)
        let dBParams     = SearchDB'.getSDBParams memoryDB
        logger.Trace (sprintf "DB parameters: %A" dBParams)

        let calcIonSeries aal  =
            Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal

        let rnd = new System.Random()

        let dbLookUp = SearchDB.getThreadSafePeptideLookUpFromFileBy memoryDB dBParams
        logger.Trace "Finished preparing processing functions."

        // initialize Reader and Transaction
        logger.Trace "Init connection to input data base."
        let inReader = Core.MzIO.Reader.getReader instrumentOutput :?> MzIO.MzSQL.MzSQL
        inReader.Connection.Open()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader
        logger.Trace (sprintf "Run ID: %s" inRunID)
        let inTr = inReader.BeginTransaction()

   
        logger.Trace "Starting charge state determination."
        let ms2sAndAssignedCharges,peakPosStdDev = getPrecursorCharge chargeParams rnd inRunID inReader
        logger.Trace "Finished charge state determination."
        
        logger.Trace "Starting peptide spectrum matching."
        psm processParams peakPosStdDev dbLookUp  calcIonSeries inReader outFilePath ms2sAndAssignedCharges
        logger.Trace "Finished peptide spectrum matching."

        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()
        memoryDB.Dispose()
        logger.Trace "Done."