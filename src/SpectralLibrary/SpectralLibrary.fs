namespace ProteomIQon

open ProteomIQon.Dto
open ProteomIQon.Dto.QuantificationResult
open FSharpAux.IO
open FSharpAux.IO.SchemaReader
open BioFSharp.Mz
open BioFSharp
open FSharp.Stats
open MzIO
open BioFSharp.Mz.SearchDB
open System.Data.SQLite
open System.Data
open ProteomIQon
open Dto.QuantificationResult
open Dto
open System.IO
open System
open FSharpAux.Colors.Table.StatisticalGraphics24
open FSharp.Stats
open FSharpAux.IO.SchemaReader
open FSharp.Plotly
open BioFSharp
open Microsoft



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
            MzDelta         : float
        }

    let createIonInformation charge iontype number mOverZ intensity relIntSpec mzDelta=
        {
            Charge           = charge
            Iontype          = iontype
            Number           = number
            MassOverCharge   = mOverZ
            Intensity        = intensity
            RelIntensitySpec = relIntSpec
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
    let getThreadSafePeptideLookUpFromFileBy (cn:SQLiteConnection) (sdbParams: SearchDbParams) =
        let parseAAString = initOfModAminoAcidString sdbParams.IsotopicMod (sdbParams.FixedMods@sdbParams.VariableMods)
        let selectModsequenceByID = prepareSelectModsequenceByModSequenceID cn
        (fun id ->
                selectModsequenceByID id
                |> (createLookUpResultBy parseAAString)
        )

    ///
    let getQuantifiedPeptides (quantFilePath:string) = 
        ///
        let peptides =
            Csv.CsvReader<QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(quantFilePath,'\t',false,1)
            |> Array.ofSeq
        let filteredPeptides =
            peptides
            |> Array.filter (fun (qp: QuantificationResult) -> 
                // Filter for peptides where fit assures a good estimation of the ScanTime
                (qp.GlobalMod = 0 && qp.Params_Light |> Array.isEmpty |> not) || (qp.GlobalMod = 1 && qp.Params_Heavy |> Array.isEmpty |> not)
                )
            |> Array.filter (fun qp -> tryTargetGetScanTime qp |> Option.isSome)
            |> Array.filter (fun qp -> (getTargetScanTimeDifference qp |> abs) / (getTargetStabw qp) < 2. )
        filteredPeptides 

    let createPSMLookupMap file =
        let psms =
            Csv.CsvReader<Dto.PSMStatisticsResult>(SchemaMode=Csv.Fill).ReadFile(file,'\t',false,1)
            |> Array.ofSeq
        let psmLookupMap =
            psms
            |> Array.groupBy (fun x -> x.ModSequenceID, x.Charge)
            |> Map.ofArray
        psmLookupMap

    let createSpectralLibrary (outDir: string) (spectralLibraryParams: Domain.SpectralLibraryParams) (cn: SQLiteConnection) (instrumentOutput, (scoredPSMs:string), quantFiles) =

        let logger = Logging.createLogger (System.IO.Path.GetFileNameWithoutExtension scoredPSMs)
        logger.Trace "Copy DB into memory"
        let memoryDB = SearchDB.copyDBIntoMemory cn
        let dBParams = getSDBParamsBy memoryDB

        let quantifiedPeptides =
            getQuantifiedPeptides quantFiles
        let psmLookupMap =
            createPSMLookupMap scoredPSMs

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
        let createPeptideIons (quant: QuantificationResult []) (psmLookup: Map<(int*int),PSMStatisticsResult[]>) (matchingTolerance: float) =
            quant
            |> Array.map (fun qr ->
                let psms =
                    psmLookup.[qr.ModSequenceID,qr.Charge]
                let ionInfoForQuant: IonInformation seq=
                    psms
                    |> Seq.collect (fun psm ->
                        let sequence = peptideLookUp psm.ModSequenceID
                        let frag =
                            let ionSeries = (calcIonSeries sequence.BioSequence).TargetMasses
                            ProteomIQon.Fragmentation'.ladderElement ionSeries [1. .. (float qr.Charge)]
                            |> List.map (fun frag -> frag.MainPeak::frag.DependentPeaks)
                            |> List.concat
                        let spec = inReader.ReadSpectrumPeaks psm.PSMId
                        let maxIntSpec =
                            spec.Peaks
                            |> Seq.maxBy (fun (x: MzIO.Binary.Peak1D) -> x.Intensity)
                            |> fun x -> x.Intensity
                        let ionInfos =
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
                                                peak.Mz
                                                peak.Intensity
                                                (peak.Intensity/maxIntSpec)
                                                deltaMass
                                            )
                                    else
                                        None
                                )
                            )
                        ionInfos
                    )
                let groupedIonInfo =
                    let totalIons = ionInfoForQuant |> Seq.length
                    let maxIntFrags = ionInfoForQuant |> Seq.maxBy (fun x -> x.Intensity) |> fun x -> x.Intensity
                    ionInfoForQuant
                    |> List.ofSeq
                    |> List.groupBy (fun x -> x.Iontype, x.Number, x.Charge)
                    |> List.map snd
                    |> List.map (fun ions ->
                        {
                            Charge                            = ions.[0].Charge
                            Iontype                           = ions.[0].Iontype.ToString()
                            Number                            = ions.[0].Number
                            // Brauch man glaub ich nicht weil das in der Quantinfo steckt
                            GlobalMod                         = qr.GlobalMod
                            CountAbsolute                     = ions.Length
                            CountFraction                     = (float ions.Length)/(float totalIons)
                            MeanFragMz                        = ions |> List.averageBy (fun x -> x.MassOverCharge)
                            CvMeanFragMz                      = ions |> Seq.cvBy (fun x -> x.MassOverCharge)
                            /// calculated fragment m/z - measured fragment m/z
                            MeanMzDelta                       = ions |> List.averageBy (fun x -> x.MzDelta)
                            MaxIntensity                      = ions |> List.maxBy (fun x -> x.Intensity) |> fun x -> x.Intensity
                            MinIntensity                      = ions |> List.minBy (fun x -> x.Intensity) |> fun x -> x.Intensity
                            MeanIntensity                     = ions |> List.averageBy (fun x -> x.Intensity)
                            // meinst du hier mit total das ganze spec?
                            /// mean of intensitiy relative to the spectrum the fragment was measured in
                            MeanRelativeIntensity_Total       = ions |> List.averageBy (fun x -> x.RelIntensitySpec)
                            /// cv of intensitiy relative to the spectrum the fragment was measured in
                            CVRelativeIntensity_Total         = ions |> Seq.cvBy (fun x -> x.RelIntensitySpec)
                            /// mean of intensitiy relative to the fragments of the precursor
                            MeanRelativeIntensity_Frags       = 
                                ions
                                |> List.map (fun x ->
                                    x.Intensity/maxIntFrags
                                )
                                |> List.average
                            /// cv of intensitiy relative to the fragments of the precursor
                            CVRelativeIntensity_Frags         =
                                ions
                                |> List.map (fun x ->
                                    x.Intensity/maxIntFrags
                                )
                                |> Seq.cv
                        }
                    )
                {
                    StringSequence                              = qr.StringSequence
                    GlobalMod                                   = qr.GlobalMod
                    Charge                                      = qr.Charge
                    PepSequenceID                               = qr.PepSequenceID
                    ModSequenceID                               = qr.ModSequenceID
                    PrecursorMZ                                 = qr.PrecursorMZ
                    MeasuredMass                                = qr.MeasuredMass
                    TheoMass                                    = qr.TheoMass
                    AbsDeltaMass                                = qr.AbsDeltaMass
                    MeanPercolatorScore                         = qr.MeanPercolatorScore
                    QValue                                      = qr.QValue
                    PEPValue                                    = qr.PEPValue
                    ProteinNames                                = qr.ProteinNames
                    Quant                                       = 
                        match qr.GlobalMod with
                        | 0 -> qr.Quant_Light
                        | 1 -> qr.QuantMz_Heavy
                    RelativeQuant                               = nan
                    MeasuredApex                                =
                        match qr.GlobalMod with
                        | 0 -> qr.MeasuredApex_Light
                        | 1 -> qr.MeasuredApex_Heavy
                    RelativeMeasuredApex                        = nan
                    Seo                                         =
                        match qr.GlobalMod with
                        | 0 -> qr.Seo_Light
                        | 1 -> qr.Seo_Heavy
                    ScanTime                                    = QuantificationResult.getTargetScanTime qr
                    ElutionWidth                                = QuantificationResult.getTargetStabw qr
                    Fragments                                   = groupedIonInfo
                }
            )
            |> fun peptideIons ->
                let maxQuant =
                    peptideIons
                    |> Array.maxBy (fun x -> x.Quant)
                    |> fun x -> x.Quant
                let maxMeasuredApex =
                    peptideIons
                    |> Array.maxBy (fun x -> x.MeasuredApex)
                    |> fun x -> x.MeasuredApex
                peptideIons
                |> Array.map (fun peptideIon ->
                    {
                        peptideIon with
                            RelativeQuant = peptideIon.Quant/maxQuant
                            RelativeMeasuredApex = peptideIon.MeasuredApex/maxMeasuredApex
                    }
                )
        let peptideIons = createPeptideIons quantifiedPeptides psmLookupMap spectralLibraryParams.MatchingTolerancePPM
        Json.serializeAndWrite outFile peptideIons
            
            //psms
            //|> Seq.collect (fun psm ->
            //    let sequence = peptideLookUp psm.ModSequenceID
            //    let frag =
            //        let ionSeries = (calcIonSeries sequence.BioSequence).TargetMasses
            //        ProteomIQon.Fragmentation'.ladderElement ionSeries chargeList
            //        |> List.map (fun frag -> frag.MainPeak::frag.DependentPeaks)
            //        |> List.concat
            //    let spec = inReader.ReadSpectrumPeaks psm.PSMId
            //    let assigned =
            //        let maxIntSpec =
            //            spec.Peaks
            //            |> Seq.maxBy (fun (x: MzIO.Binary.Peak1D) -> x.Intensity)
            //            |> fun x -> x.Intensity
            //        spec.Peaks
            //        |> Seq.collect (fun (peak: MzIO.Binary.Peak1D) ->
            //            frag
            //            |> Seq.choose (fun ion ->
            //                let deltaMass = abs (ion.MassOverCharge - peak.Mz)
            //                if deltaMass <= (Mass.deltaMassByPpm matchingTolerance peak.Mz) then
            //                    Some 
            //                        (createIonInformation
            //                            ion.Charge
            //                            ion.Iontype
            //                            ion.Number
            //                            ion.MassOverCharge
            //                            peak.Intensity
            //                            (peak.Intensity/maxIntSpec)
            //                            nan
            //                            nan
            //                            psm.PepSequenceID
            //                            psm.ModSequenceID
            //                            psm.PSMId
            //                            psm.TheoMass
            //                            psm.ScanTime
            //                            sequence.StringSequence
            //                            psm.GlobalMod
            //                            psm.PercolatorScore
            //                            deltaMass
            //                        )
            //                else
            //                    None
            //            )
            //        )
            //        |> fun fragments ->
            //            let maxIntFrag =
            //                fragments
            //                |> Seq.maxBy (fun x -> x.Intensity)
            //                |> fun x -> x.Intensity
            //            fragments
            //            |> Seq.map (fun fragment ->
            //                {fragment with RelIntensityFrag = fragment.Intensity/maxIntFrag}
            //            )
            //    assigned
            //)
            //|> fun run ->
            //    let maxIntRun =
            //        run
            //        |> Seq.maxBy (fun x -> x.Intensity)
            //        |> fun x -> x.Intensity
            //    run
            //    |> Seq.map (fun fragment ->
            //        {fragment with RelIntensityRun = fragment.Intensity/maxIntRun}
            //    )

        //let ionInformations =
        //    assignIntensitiesToMasses psmFile spectralLibraryParams.ChargeList spectralLibraryParams.MatchingTolerancePPM

        //ionInformations
        //|> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
        //|> Seq.write outFile