namespace ProteomIQon

open ProteomIQon.Domain
open ProteomIQon.Core
open System.IO
open System.Collections.Generic
open BioFSharp.Mz
open ProteomIQon.Core.MzIO.Processing
open ProteomIQon.Dto
open ProteomIQon.Domain
open Deedle


module MsFraggerToPSM =

    let psms =
        Frame.ReadCsv("C:\Users\jonat\OneDrive\Desktop\psm.tsv", true, separators = "\t")
        |> Frame.mapRows (fun k row ->
            let assignedModifications =
                row.GetAs<string>("Assigned Modifications").Split(',')
                |> Array.sort
                |> Array.filter (fun s -> s <> "")
            let stringSequence = 
                row.GetAs<string>("Peptide").ToCharArray()
                |> Array.map string
            let modStringSequence =
                assignedModifications
                |> Array.fold (fun (acc: string []) (modString: string) ->
                    let modInfo = modString.Split('(') |> Array.head
                    let modPosition, modification =
                        if modInfo.StartsWith("N-term") then 0, "nTermAcetylation"
                        else
                            if modInfo.Length = 2 then
                                int modInfo.[0..0] - 1, modInfo.[1..]
                            else
                                int modInfo.[0..1] - 1, modInfo.[2..]
                    match modification with
                    | "M" ->
                        acc.[modPosition] <- "[ox]M"
                        acc
                    | "C" ->
                        acc.[modPosition] <- "[ca]C"
                        acc
                    | "nTermAcetylation" ->
                        acc.[0] <- "[ac]" + acc.[0]
                        acc
            
                ) stringSequence
                |> String.concat ""
            {
                PSMId = row.GetAs<string>("Spectrum")
                PepSequenceID = 0
                ModSequenceID = 0
                ScanTime = row.GetAs<float>("Retention") / 60.
                Charge = row.GetAs<int>("Charge")
                PrecursorMZ = row.GetAs<float>("Calibrated Observed M/Z")
                TheoMass = row.GetAs<float>("Calculated Peptide Mass")
                AbsDeltaMass = row.GetAs<float>("Delta Mass") |> abs
                IonMobility = row.GetAs<float>("Ion Mobility")
                PeptideLength = row.GetAs<int>("Peptide Length")
                MissCleavages = row.GetAs<int>("Number of Missed Cleavages")
                Expectscore = row.GetAs<float>("Expectation")
                Hyperscore = row.GetAs<float>("Hyperscore")
                Probability = row.GetAs<float>("Probability")
                StringSequence = modStringSequence
                ProteinNames = row.GetAs<string>("Protein")
                GlobalMod = 0
            }
        )
        |> Series.values



    let cn = SearchDB.getDBConnection @"C:\Users\jonat\OneDrive\Desktop\ChlamyTest.db"
    let memoryDB = SearchDB.copyDBIntoMemory cn 
    let pepDBTr = memoryDB.BeginTransaction()

    /// Prepares statement to select a ModSequence entry by Sequence
    let prepareSelectModsequenceBySequence (cn:SQLiteConnection) (tr) =
        let querystring = "SELECT * FROM ModSequence WHERE Sequence=@sequence AND GlobalMod=@globalMod"
        let cmd = new SQLiteCommand(querystring, cn, tr) 
        cmd.Parameters.Add("@sequence", System.Data.DbType.String) |> ignore
        cmd.Parameters.Add("@globalMod", System.Data.DbType.Int32) |> ignore
        (fun (sequence:string) (globalMod:int) ->        
            cmd.Parameters.["@sequence"].Value <- sequence
            cmd.Parameters.["@globalMod"].Value <- globalMod
            use reader = cmd.ExecuteReader()            
            match reader.Read() with
            | true -> 
                {|
                    ID = reader.GetInt32(0);
                    PepSequenceID = reader.GetInt32(1);
                    RealMass = reader.GetDouble(2);
                    RoundedMass = reader.GetInt64(3);
                    Sequence = reader.GetString(4);
                    GlobalMod = reader.GetInt32(5)
                |}
                |> Some
            | false -> None
        )

    let initModSeqLookup (cn:SQLiteConnection) tr =
        let selectModsequenceBySequence = prepareSelectModsequenceBySequence cn tr
        selectModsequenceBySequence

    let modSeqLookup = initModSeqLookup memoryDB pepDBTr

    let binBy (projection: 'a -> float) bandwidth (data: seq<'a>) =
        if bandwidth = 0. then raise (System.DivideByZeroException("Bandwidth cannot be 0."))
        let halfBw = bandwidth / 2.0
        let decBandwidth = decimal bandwidth
        let tmp = 
            data
            |> Seq.groupBy (fun x -> (decimal (projection x) / decBandwidth) |> float |> floor) 
            |> Seq.map (fun (k,values) -> 
                let count = (Seq.length(values)) |> float
                if k < 0. then
                    ((k  * bandwidth) + halfBw, values)   
                else
                    ((k + 1.) * bandwidth) - halfBw, values)
            |> Seq.sortBy fst
        tmp

    let idPSMs = 
        psms
        |> Seq.map (fun v -> 
            let res = modSeqLookup v.StringSequence 0
            if res.IsSome then
                {
                    v with
                        ModSequenceID = res.Value.ID
                        PepSequenceID = res.Value.PepSequenceID
                }
                |> Some
            else None
        )
        |> Seq.choose id

    let binnedPSMs =
        binBy (fun (x: PSMStatisticsResultFragpipe) -> x.IonMobility) 0.05 idPSMs
        |> Seq.toArray

    let mzliteFiles =
        System.IO.Directory.GetFiles(@"C:\Users\jonat\OneDrive\Desktop\Testfiles", "*.mzlite")

    let findClosest (targetT,targetMz) (numbers: (float*float*string) []) =
        numbers
        |> Array.minBy (fun (t,mz,_) -> abs ((t) - targetT) + abs ((mz) - targetMz))

    let inReader = Core.MzIO.Reader.getReader @"D:\Testfiles\binned_spectra_1.000.mzlite"
    Core.MzIO.Reader.openConnection inReader
    let inTr = inReader.BeginTransaction()
    let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader


    let ms = 
        inReader.ReadMassSpectra(inRunID)
        |> Seq.toArray
        |> Array.filter (fun x -> MassSpectrum.getMsLevel x = 2)
    1
    let lookupMap =
        ms
        |> Array.map (fun x -> MassSpectrum.getScanTime x, MassSpectrum.getPrecursorMZ x, x.ID)
    lookupMap |>Array.last
    let psmsWithSpecID =
        //binnedPSMs.[0..1]
        idPSMs
        |> Array.ofSeq
        |> Array.mapi (fun i (psms) ->
            //psms
            //|> Array.ofSeq
            //|> Array.mapi (fun i psm ->
            if i%1000 = 0 then printfn "%A" i 
            //    let closest = findClosest psm.ScanTime ms
            //    {
            //        psm with
            //            PSMId = (closest |> snd)
            //    }
            //)
            let closest = findClosest (psms.ScanTime,psms.PrecursorMZ) lookupMap
            {
                psms with
                    PSMId = (closest |> fun (_,_,a) -> a)
            }
        )
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (@"D:\Testfiles\binned_spectra_1.000.mzlite" + ".qpsm")