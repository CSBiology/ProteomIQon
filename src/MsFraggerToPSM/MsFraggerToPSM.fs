namespace ProteomIQon

open System.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute
open Deedle
open System.Data
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open System.Data.SQLite
open ProteomIQon
open MzIO.Processing
open ProteomIQon.Dto
open ProteomIQon.Domain


module MsFraggerToPSM =

    let readFragegrPsms (path: string) =
        Frame.ReadCsv(path, true, separators = "\t")
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
                    | m ->
                        printfn "Modification %s is not yet implemented" m
                        Array.append acc [|"UnknownMod"|]
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
        |> Seq.filter (fun x -> x. StringSequence.Contains "UnknownMod" |> not)

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

    let findClosest (targetT,targetMz) (numbers: (float*float*string) []) =
        numbers
        |> Array.minBy (fun (t,mz,_) -> abs ((t) - targetT) + abs ((mz) - targetMz))

    let convertToPSM (inputFileMsFragger: string) (inputFileMzLite: string) (outputDir: string) (cn: SQLite.SQLiteConnection) =

        let psms = readFragegrPsms inputFileMsFragger

        let outputPath = Path.Combine(outputDir, Path.GetFileNameWithoutExtension(inputFileMsFragger) + ".qpsm")

        let memoryDB = SearchDB.copyDBIntoMemory cn
        let pepDBTr = memoryDB.BeginTransaction()

        let modSeqLookup = initModSeqLookup memoryDB pepDBTr

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

        let inReader = Core.MzIO.Reader.getReader inputFileMzLite
        Core.MzIO.Reader.openConnection inReader
        let inTr = inReader.BeginTransaction()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader


        let ms = 
            inReader.ReadMassSpectra(inRunID)
            |> Seq.toArray
            |> Array.filter (fun x -> MassSpectrum.getMsLevel x = 2)
        let lookupMap =
            ms
            |> Array.map (fun x -> MassSpectrum.getScanTime x, MassSpectrum.getPrecursorMZ x, x.ID)

        let psmsWithSpecID =
            idPSMs
            |> Array.ofSeq
            |> Array.mapi (fun i (psms) ->
                let closest = findClosest (psms.ScanTime,psms.PrecursorMZ) lookupMap
                {
                    psms with
                        PSMId = (closest |> fun (_,_,a) -> a)
                }
            )

        psmsWithSpecID
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outputPath)