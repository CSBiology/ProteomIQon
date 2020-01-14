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
            Charge        : float
            Iontype       : Ions.IonTypeFlag
            MassOverCharge: float
            Number        : int
            Intensity     : float
            ModSequenceID : int
            PSMId         : string
        }
    
    let createIonInformation charge iontype mOverZ number intensity modSeqID psmID =
        {
            Charge         = charge
            Iontype        = iontype
            MassOverCharge = mOverZ
            Number         = number
            Intensity      = intensity
            ModSequenceID  = modSeqID
            PSMId          = psmID
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

    let createSpectralLibrary (chargeList: float list) (matchingTolerance: float) (cn: SQLiteConnection) (instrumentOutput: string) (scoredPSMs: string)=

        let cn = SearchDB.getDBConnection @"C:\Users\jonat\source\repos\davedata\PSMStatisticsTest\Chlamy15N.db"
        let memoryDB = SearchDB.copyDBIntoMemory cn
        let dBParams = getSDBParamsBy memoryDB

        let peptideLookUp = getThreadSafePeptideLookUpFromFileBy memoryDB dBParams

        let calcIonSeries aal  =
            Fragmentation.Series.fragmentMasses Fragmentation.Series.bOfBioList Fragmentation.Series.yOfBioList dBParams.MassFunction aal

        let inReader = Core.MzIO.Reader.getReader instrumentOutput
        let inRunID = Core.MzIO.Reader.getDefaultRunID inReader
        inReader.BeginTransaction()

        let psmFile =
            Seq.fromFileWithCsvSchema<PSMStatisticsResult>(scoredPSMs, '\t', false,schemaMode = FSharpAux.IO.SchemaReader.Csv.SchemaModes.Fill, skipLines = 1)
            |> Seq.toArray

        let assignIntensitiesToMasses (psms: PSMStatisticsResult []) (chargeList: float list) (matchingTolerance: float) =
            psms
            |> Array.map (fun psm -> 
                let sequence = peptideLookUp psm.ModSequenceID
                let frag = 
                    let ionSeries = (calcIonSeries sequence.BioSequence).TargetMasses
                    ProteomIQon.Fragmentation'.ladderElement ionSeries chargeList
                    |> List.map (fun frag -> frag.MainPeak::frag.DependentPeaks)
                    |> List.concat
                let spec = inReader.ReadSpectrumPeaks psm.PSMId
                let assigned = 
                    spec.Peaks
                    |> Seq.collect (fun (peak: MzIO.Binary.Peak1D) ->
                        frag
                        |> List.choose (fun ion ->
                            if (abs (ion.MassOverCharge - peak.Mz)) <= (Mass.deltaMassByPpm matchingTolerance peak.Mz) then
                                Some (createIonInformation ion.Charge ion.Iontype ion.MassOverCharge ion.Number peak.Intensity psm.ModSequenceID psm.PSMId)
                            else
                                None
                        )
                    )
                assigned
            )

        assignIntensitiesToMasses psmFile chargeList matchingTolerance