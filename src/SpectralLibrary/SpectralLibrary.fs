namespace ProteomIQon

open ProteomIQon.Dto
open FSharpAux.IO
open BioFSharp.Mz
open BioFSharp
open MzIO
open BioFSharp.Mz.SearchDB
open System.Data.SQLite
open System.Data
open ProteomIQon

module SpectralLibrary =

    /// Holds information about ion and in which spectrum it is found.
    type IonInformation =
        {
            Charge          : float
            Iontype         : Ions.IonTypeFlag
            Number          : int
            MassOverCharge  : float
            Intensity       : float
            RelIntensitySpec: float
            RelIntensityFrag: float
            RelIntensityRun : float
            PepSequenceID   : int
            ModSequenceID   : int
            PSMId           : string
            PrecursorMZ     : float
            ScanTime        : float
            Sequence        : string
            GlobalMod       : int
            PercolatorScore : float
            MzDelta         : float
        }

    let createIonInformation charge iontype number mOverZ intensity relIntSpec relIntFrag relIntRun pepSeqID modSeqID psmID precMZ scanTime sequence globalMod percolatorScore mzDelta=
        {
            Charge           = charge
            Iontype          = iontype
            Number           = number
            MassOverCharge   = mOverZ
            Intensity        = intensity
            RelIntensitySpec = relIntSpec
            RelIntensityFrag = relIntFrag
            RelIntensityRun  = relIntRun
            PepSequenceID    = pepSeqID
            ModSequenceID    = modSeqID
            PSMId            = psmID
            PrecursorMZ      = precMZ
            ScanTime         = scanTime
            Sequence         = sequence
            GlobalMod        = globalMod
            PercolatorScore  = percolatorScore
            MzDelta          = mzDelta
        }

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

    /// Prepares statement to select a ModSequence entry by ModSequenceID
    let prepareSelectModsequenceByModSequenceID (cn:SQLiteConnection) =
        let querystring = "SELECT * FROM ModSequence WHERE ID=@id"
        let cmd = new SQLiteCommand(querystring, cn)
        cmd.Parameters.Add("@id", DbType.Int64) |> ignore
        (fun (id:int)  ->
            cmd.Parameters.["@id"].Value <- id

            use reader = cmd.ExecuteReader()
            match reader.Read() with
            | true ->  (reader.GetInt32(0), reader.GetInt32(1),reader.GetDouble(2), reader.GetInt64(3), reader.GetString(4), reader.GetInt32(5))
            | false -> -1,-1,nan,-1L,"",-1
        )

    /// Returns a LookUpResult list
    let getThreadSafePeptideLookUpFromFileBy (cn:SQLiteConnection) sdbParams =
        let parseAAString = initOfModAminoAcidString sdbParams.IsotopicMod (sdbParams.FixedMods@sdbParams.VariableMods)
        let selectModsequenceByID = prepareSelectModsequenceByModSequenceID cn
        (fun id ->
                selectModsequenceByID id
                |> (createLookUpResultBy parseAAString)
        )

    let createSpectralLibrary (outDir: string) (spectralLibraryParams: Domain.SpectralLibraryParams) (cn: SQLiteConnection) (instrumentOutputAndScoredPSMs: string*string)  =


        let instrumentOutput, scoredPSMs = instrumentOutputAndScoredPSMs
        let logger = Logging.createLogger (System.IO.Path.GetFileNameWithoutExtension scoredPSMs)
        logger.Trace "Copy DB into memory"
        let memoryDB = SearchDB.copyDBIntoMemory cn
        let dBParams = getSDBParamsBy memoryDB

        logger.Trace (sprintf "Creating Library using the following files:\n%s\n%s" instrumentOutput scoredPSMs)

        let outFile = sprintf @"%s\%s.sl" outDir (System.IO.Path.GetFileNameWithoutExtension scoredPSMs)

        let peptideLookUp = getThreadSafePeptideLookUpFromFileBy memoryDB dBParams

        let calcIonSeries aal  =
            Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal

        let inReader = Core.MzIO.Reader.getReader instrumentOutput
        //let inRunID = Core.MzIO.Reader.getDefaultRunID inReader
        logger.Trace "Reading instrument output"
        let inTr = inReader.BeginTransaction()
        logger.Trace "Reading psm file"
        let psmFile =
            Seq.fromFileWithCsvSchema<PSMStatisticsResult>(scoredPSMs, '\t', false,schemaMode = FSharpAux.IO.SchemaReader.Csv.SchemaModes.Fill, skipLines = 1)
            |> Seq.toArray
        logger.Trace "Creating library"
        let assignIntensitiesToMasses (psms: PSMStatisticsResult []) (chargeList: float list) (matchingTolerance: float) =
            psms
            |> Seq.collect (fun psm ->
                let sequence = peptideLookUp psm.ModSequenceID
                let frag =
                    let ionSeries = (calcIonSeries sequence.BioSequence).TargetMasses
                    ProteomIQon.Fragmentation'.ladderElement ionSeries chargeList
                    |> List.map (fun frag -> frag.MainPeak::frag.DependentPeaks)
                    |> List.concat
                let spec = inReader.ReadSpectrumPeaks psm.PSMId
                let assigned =
                    let maxIntSpec =
                        spec.Peaks
                        |> Seq.maxBy (fun (x: MzIO.Binary.Peak1D) -> x.Intensity)
                        |> fun x -> x.Intensity
                    spec.Peaks
                    |> Seq.collect (fun (peak: MzIO.Binary.Peak1D) ->
                        frag
                        |> Seq.choose (fun ion ->
                            let deltaMass = abs (ion.MassOverCharge - peak.Mz)
                            if deltaMass <= (Mass.deltaMassByPpm matchingTolerance peak.Mz) then
                                Some 
                                    (createIonInformation
                                        ion.Charge
                                        ion.Iontype
                                        ion.Number
                                        ion.MassOverCharge
                                        peak.Intensity
                                        (peak.Intensity/maxIntSpec)
                                        nan
                                        nan
                                        psm.PepSequenceID
                                        psm.ModSequenceID
                                        psm.PSMId
                                        psm.TheoMass
                                        psm.ScanTime
                                        sequence.StringSequence
                                        psm.GlobalMod
                                        psm.PercolatorScore
                                        deltaMass
                                    )
                            else
                                None
                        )
                    )
                    |> fun fragments ->
                        let maxIntFrag =
                            fragments
                            |> Seq.maxBy (fun x -> x.Intensity)
                            |> fun x -> x.Intensity
                        fragments
                        |> Seq.map (fun fragment ->
                            {fragment with RelIntensityFrag = fragment.Intensity/maxIntFrag}
                        )
                assigned
            )
            |> fun run ->
                let maxIntRun =
                    run
                    |> Seq.maxBy (fun x -> x.Intensity)
                    |> fun x -> x.Intensity
                run
                |> Seq.map (fun fragment ->
                    {fragment with RelIntensityRun = fragment.Intensity/maxIntRun}
                )

        let ionInformations =
            assignIntensitiesToMasses psmFile spectralLibraryParams.ChargeList spectralLibraryParams.MatchingTolerancePPM

        ionInformations
        |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
        |> Seq.write outFile