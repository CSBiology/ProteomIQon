namespace ProteomIQon

open FSharpAux.IO

module ConsensusSpectralLibrary =

    type IonInformation =
        {
            Charge         : float
            Iontype        : string
            MassOverCharge : float
            Number         : int
            Intensity      : float
            PepSequenceID  : int
            ModSequenceID  : int
            PSMId          : string
            PrecursorMZ    : float
            ScanTime       : float
            Sequence       : string
            GlobalMod      : int
            PercolatorScore: float
        }

    type ConsensIonInformation =
        {
            Charge         : float
            Iontype        : string
            MassOverCharge : float
            Number         : int
            Intensity      : float
            PepSequenceID  : int
            ModSequenceID  : int
            PSMId          : string
            PrecursorMZ    : float
            ScanTime       : float
            Count          : int
            Version        : int
            Sequence       : string
            GlobalMod      : int
            PercolatorScore: float
        }
        static member create charge ionType massOverCharge number intensity pepSeqID modSeqID psmID precMZ scanTime count version sequence globalMod percolatorScore =
            {
                Charge         = charge
                Iontype        = ionType
                MassOverCharge = massOverCharge
                Number         = number
                Intensity      = intensity
                PepSequenceID  = pepSeqID
                ModSequenceID  = modSeqID
                PSMId          = psmID
                PrecursorMZ    = precMZ
                ScanTime       = scanTime
                Count          = count
                Version        = version
                Sequence       = sequence
                GlobalMod      = globalMod
                PercolatorScore= percolatorScore
            }
    
        //static member createFromFile charge (ionType: string) massOverCharge number intensity pepSeqID modSeqID psmID precMZ scanTime count version sequence globalMod percolatorScore =
        //    {
        //        Charge         = charge
        //        Iontype        =
        //            ionType.Replace(" ","").Split(',')
        //            |> Array.sort
        //            |> String.concat ","
        //        MassOverCharge = massOverCharge
        //        Number         = number
        //        Intensity      = intensity
        //        PepSequenceID  = pepSeqID
        //        ModSequenceID  = modSeqID
        //        PSMId          = psmID
        //        PrecursorMZ    = precMZ
        //        ScanTime       = scanTime
        //        Count          = count
        //        Version        = version
        //        Sequence       = sequence
        //        GlobalMod      = globalMod
        //        PercolatorScore= percolatorScore
        //    }
  
    let binning (difference: float) (data: IonInformation[]) =
        let doubleDiff = difference * 2.
        data
        |> Array.groupBy (fun x -> floor (x.ScanTime / doubleDiff))
        |> Array.map snd

    let checkConflict (field: string) (arr: 'a[]) =
        arr
        |> Array.distinct
        |> fun x ->
            if x.Length > 1 then
                failwith (sprintf "Version conflict %s" field)
            else
                x |> Array.head

    let buildConsens (paths: string []) tolerance (outPath: string)=

        let entries =
            paths
            |> Array.collect (fun path ->
                Seq.fromFileWithCsvSchema<IonInformation>(path, '\t', true,schemaMode = FSharpAux.IO.SchemaReader.Csv.SchemaModes.Fill)
                |> Seq.toArray
            )
            |> Array.map (fun entry ->
                {
                    entry with
                        Iontype = 
                            entry.Iontype.Replace(" ","").Split(',')
                            |> Array.sort
                            |> String.concat ","
                }
            )
        let grouped =
            entries
            |> Array.groupBy (fun entry -> entry.Charge, entry.Iontype, entry.Number, entry.ModSequenceID, entry.GlobalMod)
            |> Array.map snd
        let binned =
            grouped
            |> Array.map (fun groupedEntries ->
                groupedEntries
                |> binning tolerance
            )
        let consens =
            binned
            |> Array.map (fun group -> 
                group
                |> Array.mapi (fun i version ->
                    let charge =
                        version
                        |> Array.map (fun entry -> entry.Charge)
                        |> (checkConflict "Charge")
                    let ionType =
                        version
                        |> Array.map (fun entry -> entry.Iontype)
                        |> (checkConflict "Iontype")
                    let massCharge =
                        version
                        |> Array.map (fun entry -> entry.MassOverCharge)
                        |> (checkConflict "m/z")
                    let number =
                        version
                        |> Array.map (fun entry -> entry.Number)
                        |> (checkConflict "Number")
                    let intensity =
                        version
                        |> Array.averageBy (fun entry -> entry.Intensity)
                    let pepSeqID =
                        version
                        |> Array.map (fun entry -> entry.PepSequenceID)
                        |> (checkConflict "PepSequenceID")
                    let modSeqID =
                        version
                        |> Array.map (fun entry -> entry.ModSequenceID)
                        |> (checkConflict "ModSequenceID")
                    let psmId =
                        version
                        |> Array.map (fun entry -> entry.PSMId)
                        |> Array.distinct
                        |> String.concat ";"
                    let precursorMz =
                        version
                        |> Array.map (fun entry -> entry.PrecursorMZ)
                        |> (checkConflict "PrecursorMZ")
                    let scanTime =
                        version
                        |> Array.averageBy (fun entry -> entry.ScanTime)
                    let sequence =
                        version
                        |> Array.map (fun entry -> entry.Sequence)
                        |> (checkConflict "Sequence")
                    let globalMod =
                        version
                        |> Array.map (fun entry -> entry.GlobalMod)
                        |> (checkConflict "GlobalMod")
                    let percolatorScore =
                        version
                        |> Array.averageBy (fun entry -> entry.PercolatorScore)
                    ConsensIonInformation.create
                        charge
                        ionType
                        massCharge
                        number
                        intensity
                        pepSeqID
                        modSeqID
                        psmId
                        precursorMz
                        scanTime
                        version.Length
                        i
                        sequence
                        globalMod
                        percolatorScore
                )
            )
            |> Array.concat
        consens
        |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
        |> Seq.write (sprintf @"%s/Consens.csl"outPath)
        //let head =
        //    Seq.fromFileWithCsvSchema<IonInformation>(paths.[0], '\t', true,schemaMode = FSharpAux.IO.SchemaReader.Csv.SchemaModes.Fill)
        //    |> Seq.toArray
        //    |> Array.map (fun ionInfo ->
        //        ConsensIonInformation.createFromFile
        //            ionInfo.Charge
        //            ionInfo.Iontype
        //            ionInfo.MassOverCharge
        //            ionInfo.Number
        //            ionInfo.Intensity
        //            ionInfo.PepSequenceID
        //            ionInfo.ModSequenceID
        //            ionInfo.PSMId
        //            ionInfo.PrecursorMZ
        //            ionInfo.ScanTime
        //            1.
        //            0
        //            ionInfo.Sequence
        //            ionInfo.GlobalMod
        //            ionInfo.PercolatorScore
        //    )
        //let rec loop i (acc: ConsensIonInformation []) =
        //    if i < paths.Length then
        //        let currentFile =
        //            Seq.fromFileWithCsvSchema<IonInformation>(paths.[i], '\t', true,schemaMode = FSharpAux.IO.SchemaReader.Csv.SchemaModes.Fill)
        //            |> Seq.toArray
        //            |> Array.map (fun ionInfo ->
        //                ConsensIonInformation.createFromFile
        //                    ionInfo.Charge
        //                    ionInfo.Iontype
        //                    ionInfo.MassOverCharge
        //                    ionInfo.Number
        //                    ionInfo.Intensity
        //                    ionInfo.PepSequenceID
        //                    ionInfo.ModSequenceID
        //                    ionInfo.PSMId
        //                    ionInfo.PrecursorMZ
        //                    ionInfo.ScanTime
        //                    1.
        //                    0
        //                    ionInfo.Sequence
        //                    ionInfo.GlobalMod
        //                    ionInfo.PercolatorScore
        //            )
        //        let total =
        //            Array.append acc currentFile
        //            |> Array.groupBy (fun info -> info.Charge, info.Iontype, info.Number, info.ModSequenceID, info.GlobalMod)
        //            |> Array.map snd
        //        let newAcc =
        //            total
        //            |> Array.map (fun sameIons ->
        //                if sameIons.Length = 1 then
        //                    [|sameIons.[0]|]
        //                else
        //                    //bin by rt with binsize = allowed difference
        //                    let binned =
        //                        sameIons
        //                        |> binning tolerance
        //                    binned
        //                    |> Array.mapi (fun i bin ->
        //                        bin.[1 ..]
        //                        |> Array.fold (fun acc ionInfo ->
        //                            {
        //                                acc with
        //                                    PSMId           = 
        //                                        if acc.PSMId <> ionInfo.PSMId then
        //                                            acc.PSMId + ";" + ionInfo.PSMId
        //                                        else
        //                                            acc.PSMId
        //                                    MassOverCharge  = (acc.MassOverCharge * acc.Count + ionInfo.MassOverCharge * ionInfo.Count) / (acc.Count + ionInfo.Count)
        //                                    Intensity       = (acc.Intensity * acc.Count + ionInfo.Intensity * ionInfo.Count) / (acc.Count + ionInfo.Count)
        //                                    ScanTime        = (acc.ScanTime * acc.Count + ionInfo.ScanTime * ionInfo.Count) / (acc.Count + ionInfo.Count)
        //                                    Count           = acc.Count + ionInfo.Count
        //                                    Version         = i
        //                                    PercolatorScore = (acc.PercolatorScore * acc.Count + ionInfo.PercolatorScore * ionInfo.Count) / (acc.Count + ionInfo.Count)
        //                            }
        //                        )bin.[0]
        //                    )
        //            )
        //            |> Array.concat
        //        loop (i+1) newAcc
        //    else
        //        acc
        //loop 1 head
        //|> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
        //|> Seq.write (sprintf @"%s/Consens.csl"outPath)
