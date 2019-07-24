namespace ProteomIQon

open System.Data.SQLite
open BioFSharp
open BioFSharp.Mz.SearchDB
open Domain
open Core 
open Logary
open System.IO
open BioFSharp.Mz
open MzLite 
open MzLite.Model
open MzLite.Binary
open Core.MzLite.Reader
open Core.MzLite.Peaks
open Core.MzLite
open FSharpAux.IO

module PeptideSpectrumMatching = 

    open System.IO
    open System.Data
    open MzLite.IO
    open BioFSharp.Mz.TheoreticalSpectra
    open BioFSharp.Mz.ChargeState

    /// Returns SearchDbParams of a existing database by filePath
    let getSDBParamsBy (cn :SQLiteConnection)= 
        let cn = 
            match cn.State with 
            | ConnectionState.Open ->
                cn
            | ConnectionState.Closed -> 
                cn.Open()
                cn
            | _ as x -> failwith "Data base is busy."
        match Db.SQLiteQuery.selectSearchDbParams cn with 
        | Some (iD,name,fo,fp,pr,minmscl,maxmscl,mass,minpL,maxpL,isoL,mMode,fMods,vMods,vThr) -> 
            createSearchDbParams 
                name fo fp id (Digestion.Table.getProteaseBy pr) minmscl maxmscl mass minpL maxpL 
                    (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchInfoIsotopic list>(isoL)) (Newtonsoft.Json.JsonConvert.DeserializeObject<MassMode>(mMode)) (massFBy (Newtonsoft.Json.JsonConvert.DeserializeObject<MassMode>(mMode))) 
                        (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchModification list>(fMods)) (Newtonsoft.Json.JsonConvert.DeserializeObject<SearchModification list>(vMods)) vThr
        | None ->
            failwith "This database does not contain any SearchParameters. It is not recommended to work with this file."

    let getPrecursorCharge (chParams:ChargeState.ChargeDetermParams) rnd inRunID (inReader: IMzLiteDataReader) =
        /// Returns a Sequence containing all MassSpectra of a single MS-run
        let massSpectra = inReader.ReadMassSpectra(inRunID)
        /// Returns a Array that contains all MS1s sorted by their scanTime
        let ms1SortedByScanTime =
            massSpectra
            |> Seq.filter (fun ms -> MassSpectrum.getMsLevel ms = 1)
            |> Seq.sortBy  Core.MzLite.MassSpectrum.getScanTime
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
                                                        for i = chParams.ExpectedMinimalCharge to chParams.ExpectedMaximumCharge do
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
        let init = ChargeState.initMzDevOfRndSpec rnd chParams peakPosStdDev//peakPosStdDev
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
                                                        for i = chParams.ExpectedMinimalCharge to chParams.ExpectedMaximumCharge do
                                                            let precMz = MassSpectrum.getPrecursorMZ (inReader.ReadMassSpectrum ms2ID)
                                                            let mass = Mass.ofMZ precMz (float i)
                                                            let score = ChargeState.getScore 10 1 100.
                                                            yield ChargeState.createAssignedCharge ms1ID ms2ID precMz i mass 100. score [0.] 0 0. (Set[]) None
                                                    ]
                                            | _  -> charges
                                        )
                                items
                            )

        positionMetricScoredCharges

    let psm (processParams:PeptideSpectrumMatchingParams) lookUpF calcIonSeries (reader: IMzLiteDataReader) outFilePath (ms2sAndAssignedCharges: AssignedCharge list list) =
        let resultWriter = new System.IO.StreamWriter(outFilePath, true)
        let (ms2IDAssignedCharge) =
                ms2sAndAssignedCharges
                |> List.mapi (fun i (assignedCharges) ->
                                assignedCharges |> List.map (fun assCh -> i,assCh.ProductSpecID, assCh)
                            )
                |> List.concat
                |> List.groupBy (fun (ascendingID,ms2Id,assCH) -> assCH.PrecCharge)
                |> List.map (fun (ch,ms2IDassCHL) ->
                                ch, ms2IDassCHL |> List.sortBy (fun (ascendingID,ms2ID,assCh) -> assCh.PutMass )
                            )
                |> List.sortBy fst
        ms2IDAssignedCharge 
        |> List.iter (fun (ch, ms2IdAssCh) -> printfn "%i spectra with charge %i" ms2IdAssCh.Length ch )

        ms2IDAssignedCharge
        |> List.iter (fun (ch,ms2IdAssCh) ->
                        printfn "%i spectra with charge %i processed" ms2IdAssCh.Length ch
                        ms2IdAssCh
                        |> List.iteri (fun i (ascendingID,ms2Id,assCh) ->
                                            if i%10000 = 0 then printfn "%i" i
                                            try                                                
                                            let scanTime = MassSpectrum.getScanTime (reader.ReadMassSpectrum(ms2Id))
                                            let recSpec =
                                                Peaks.unzipIMzliteArray (reader.ReadSpectrumPeaks(ms2Id).Peaks)
                                                |> fun (mzData,intensityData) -> PeakArray.zipMzInt mzData intensityData
                                            let lowerMass,upperMass =
                                                let massWithH2O = assCh.PutMass
                                                Mass.rangePpm processParams.LookUpPPM massWithH2O
                                            let lookUpResults :SearchDB.LookUpResult<AminoAcids.AminoAcid> list = 
                                                lookUpF lowerMass upperMass
                                            let theoSpecs = 
                                                lookUpResults
                                                |> List.map (fun lookUpResult -> 
                                                                let ionSeries = calcIonSeries lookUpResult.BioSequence
                                                                lookUpResult,ionSeries
                                                            )   
                                            let sequestTheoreticalSpecs = SequestLike.getTheoSpecs processParams.MS2ScanRange assCh.PrecCharge theoSpecs   
                                            let sequestLikeScored = 
                                                SequestLike.calcSequestScore processParams.MS2ScanRange recSpec scanTime assCh.PrecCharge 
                                                    assCh.PrecursorMZ sequestTheoreticalSpecs ms2Id

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
                                                |> AndromedaLike.getTheoSpecs processParams.MS2ScanRange assCh.PrecCharge    
                                                
                                            let andromedaLikeScored = 
                                                AndromedaLike.calcAndromedaScore processParams.AndromedaParams.PMinPMax processParams.MS2ScanRange processParams.AndromedaParams.MatchingIonTolerancePPM
                                                    recSpec scanTime assCh.PrecCharge assCh.PrecursorMZ andromedaTheorticalSpecs ms2Id

                                            let result = 
                                                andromedaLikeScored
                                                |> List.choose (fun androRes -> 
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
                                                                                //ScanTime                     = androRes.ScanTime
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
                                                                                //ScanTime                     = androRes.ScanTime
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
                                                                                StringSequence               = androRes.StringSequence
                                                                                ProteinNames                 = "PlaceHolder"
                                                                            }
                                                                        Some (label,res)
                                                                    | None -> None 
                                                            )
                                                |> List.groupBy (fun x -> fst x)
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

                                                    printfn "fail: \n %s" ex.Message
                                                    ()
                                        )


                        )
        resultWriter.Dispose()
    //    Logger.printTimenWithPre pre "Finished PSM"

    let scoreSpectra (processParams:PeptideSpectrumMatchingParams) (outputDir:string) (cn:SQLiteConnection) (instrumentOutput:string) =
        printfn "Now performing peptide spectrum matching: %s Results will be written to: %s" instrumentOutput outputDir

        // initialize Reader and Transaction
        let outFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".psm"
            Path.Combine [|outputDir;fileName|]

        printfn "Copy peptide DB into Memory"
        let memoryDB = SearchDB.copyDBIntoMemory cn 
        printfn "Copy peptide DB into Memory: finished"

        printfn "Prepare processing functions."
        
        let chargeParams = processParams.ChargeStateDeterminationParams
        let dBParams     = getSDBParamsBy memoryDB      

        let calcIonSeries aal  =
            Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal

        let rnd = new System.Random()

        let dbLookUp = SearchDB.getThreadSafePeptideLookUpFromFileBy memoryDB dBParams
        printfn "Finished preparing processing functions."
                
        // initialize Reader and Transaction
        printfn "Init connection to input data base." 
        let inReader = Core.MzLite.Reader.getReader instrumentOutput  
        let inRunID  = Core.MzLite.Reader.getDefaultRunID inReader
        let inTr = inReader.BeginTransaction()                    
     
            
        // Charge state determination 
        printfn "Starting charge state determination."
        let ms2sAndAssignedCharges = getPrecursorCharge chargeParams rnd inRunID inReader
        printfn "Finished charge state determination."
        
        printfn "Starting peptide spectrum matching."
        psm processParams dbLookUp calcIonSeries inReader outFilePath ms2sAndAssignedCharges
        printfn "Finished peptide spectrum matching."
        
        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()
        printfn "Done."
