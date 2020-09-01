namespace ProteomIQon

open FSharpAux.IO

module ConsensusSpectralLibrary =

    type IonInformation =
        {
            Charge        : float
            Iontype       : string
            MassOverCharge: float
            Number        : int
            Intensity     : float
            ModSequenceID : int
            PSMId         : string
            PrecursorMZ   : float
            ScanTime      : float
            Sequence      : string
            GlobalMod     : int
        }

    type ConsensIonInformation =
        {
            Charge        : float
            Iontype       : string
            MassOverCharge: float
            Number        : int
            Intensity     : float
            ModSequenceID : int
            PSMId         : string
            PrecursorMZ   : float
            ScanTime      : float
            Count         : float
            Version       : int
            Sequence      : string
            GlobalMod     : int
        }
        static member create charge ionType massOverCharge number intensity modSeqID psmID precMZ scanTime count version sequence globalMod=
            {
                Charge         = charge
                Iontype        = ionType
                MassOverCharge = massOverCharge
                Number         = number
                Intensity      = intensity
                ModSequenceID  = modSeqID
                PSMId          = psmID
                PrecursorMZ    = precMZ
                ScanTime       = scanTime
                Count          = count
                Version        = version
                Sequence       = sequence
                GlobalMod      = globalMod
            }
    
        static member createFromFile charge (ionType: string) massOverCharge number intensity modSeqID psmID precMZ scanTime count version sequence globalMod=
            {
                Charge         = charge
                Iontype        =
                    ionType.Replace(" ","").Split(',')
                    |> Array.sort
                    |> String.concat ","
                MassOverCharge = massOverCharge
                Number         = number
                Intensity      = intensity
                ModSequenceID  = modSeqID
                PSMId          = psmID
                PrecursorMZ    = precMZ
                ScanTime       = scanTime
                Count          = count
                Version        = version
                Sequence       = sequence
                GlobalMod      = globalMod
            }
  
    let binning (difference: float) (data: ConsensIonInformation[]) =
        let doubleDiff = difference * 2.
        data
        |> Array.groupBy (fun x -> floor (x.ScanTime / doubleDiff))
        |> Array.map snd


    let buildConsens (paths: string []) tolerance (outPath: string)=
        let head =
            Seq.fromFileWithCsvSchema<IonInformation>(paths.[0], '\t', true,schemaMode = FSharpAux.IO.SchemaReader.Csv.SchemaModes.Fill)
            |> Seq.toArray
            |> Array.map (fun ionInfo ->
                ConsensIonInformation.createFromFile
                    ionInfo.Charge
                    ionInfo.Iontype
                    ionInfo.MassOverCharge
                    ionInfo.Number
                    ionInfo.Intensity
                    ionInfo.ModSequenceID
                    ionInfo.PSMId
                    ionInfo.PrecursorMZ
                    ionInfo.ScanTime
                    1.
                    0
                    ionInfo.Sequence
                    ionInfo.GlobalMod
            )
        let rec loop i (acc: ConsensIonInformation []) =
            let currentFile =
                Seq.fromFileWithCsvSchema<IonInformation>(paths.[i], '\t', true,schemaMode = FSharpAux.IO.SchemaReader.Csv.SchemaModes.Fill)
                |> Seq.toArray
                |> Array.map (fun ionInfo ->
                    ConsensIonInformation.createFromFile
                        ionInfo.Charge
                        ionInfo.Iontype
                        ionInfo.MassOverCharge
                        ionInfo.Number
                        ionInfo.Intensity
                        ionInfo.ModSequenceID
                        ionInfo.PSMId
                        ionInfo.PrecursorMZ
                        ionInfo.ScanTime
                        1.
                        0
                        ionInfo.Sequence
                        ionInfo.GlobalMod
                )
            let total =
                Array.append acc currentFile
                |> Array.groupBy (fun info -> info.Charge, info.Iontype, info.Number, info.ModSequenceID, info.GlobalMod)
                |> Array.map snd
            total
            |> Array.map (fun sameIons ->
                if sameIons.Length = 1 then
                    [|sameIons.[0]|]
                else
                    //bin by rt with binsize = allowed difference
                    let binned =
                        sameIons
                        |> binning tolerance
                    binned
                    |> Array.mapi (fun i bin ->
                        bin.[1 ..]
                        |> Array.fold (fun acc ionInfo ->
                            {
                                acc with
                                    MassOverCharge = (acc.MassOverCharge * acc.Count + ionInfo.MassOverCharge * ionInfo.Count) / (acc.Count + ionInfo.Count)
                                    Intensity      = (acc.Intensity * acc.Count + ionInfo.Intensity * ionInfo.Count) / (acc.Count + ionInfo.Count)
                                    ScanTime       = (acc.ScanTime * acc.Count + ionInfo.ScanTime * ionInfo.Count) / (acc.Count + ionInfo.Count)
                                    Count          = acc.Count + ionInfo.Count
                                    Version        = i
                            }
                        )bin.[0]
                    )
            )
        loop 1 head
        |> Array.concat
        |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
        |> Seq.write (sprintf @"%s/Consens.csl"outPath)
