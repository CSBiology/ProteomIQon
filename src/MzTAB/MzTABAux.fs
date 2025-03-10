namespace ProteomIQon

open System
open MzIO
open FSharpAux
open FSharpAux.IO
open Plotly.NET
open Deedle
open Domain

module MzTABAux =

    type TableSort =
        {
            Protein: string
            Experiment: string
            DistinctPeptideCount: float
            Quant_Heavy: float option
            Quant_Light: float
            QValue: float
            StudySubject: float
            Subject_SEM: float option
            Subject_StDev: float option
        }

    let createTableSort protein experiment distPepCount quantHeavy quantLight qVal studySubject subjectSEM subjectStDev =
        {
            Protein              = protein
            Experiment           = experiment
            DistinctPeptideCount = distPepCount
            Quant_Heavy          = quantHeavy
            Quant_Light          = quantLight
            QValue               = qVal
            StudySubject         = studySubject
            Subject_SEM          = subjectSEM
            Subject_StDev        = subjectStDev
        }

    type InferredProteinClassItemOut =
        {
            ProteinGroup    : string
            PeptideSequence : string
            Class           : string
            TargetScore     : float
            DecoyScore      : float
            QValue          : float
        }

    type PSMStatisticsResult = {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_
        PSMId                        : string
        GlobalMod                    : int
        PepSequenceID                : int
        ModSequenceID                : int
        Label                        : int
        // ascending ms2 id (file specific)
        ScanNr                       : int
        ScanTime                     : float
        Charge                       : int
        PrecursorMZ                  : float
        TheoMass                     : float
        AbsDeltaMass                 : float
        PeptideLength                : int
        MissCleavages                : int
        SequestScore                 : float
        SequestNormDeltaBestToRest   : float
        SequestNormDeltaNext         : float
        AndroScore                   : float
        AndroNormDeltaBestToRest     : float
        AndroNormDeltaNext           : float
        XtandemScore                 : float
        XtandemNormDeltaBestToRest   : float
        XtandemNormDeltaNext         : float
        PercolatorScore              : float
        QValue                       : float
        PEPValue                     : float
        StringSequence               : string
        ProteinNames                 : string
        }

    type ProteinSection =
        {
            accession                                 : string[]
            description                               : string
            taxid                                     : int
            species                                   : string
            database                                  : string
            database_version                          : string
            search_engine                             : string
            best_search_engine_score                  : float option
            search_engine_score_ms_run                : (int*float option)[][]
            reliability                               : int
            num_psms_ms_run                           : (int*float option)[]
            num_peptides_distinct_ms_run              : (int*float option)[]
            num_peptides_unique_ms_run                : (int*float option)[]
            ambiguity_members                         : string[]
            modifications                             : string
            uri                                       : string
            go_terms                                  : string[]
            protein_coverage                          : float
            protein_abundance_assay                   : (int*float option)[]
            protein_abundance_study_variable          : (int*float option)[]
            protein_abundance_stdev_study_variable    : (int*float option)[]
            protein_abundance_std_error_study_variable: (int*float option)[]
        }

    type PeptideSection =
        {
            sequence                                  : string
            accession                                 : string[][]
            unique                                    : int
            database                                  : string
            database_version                          : string
            search_engine                             : string
            best_search_engine_score                  : (int*float) []
            search_engine_score_ms_run                : (int*float option)[][]
            reliability                               : int
            modifications                             : string
            retention_time                            : float
            retention_time_window                     : float*float
            charge                                    : int
            mass_to_charge                            : float
            uri                                       : string
            spectra_ref                               : string
            peptide_abundance_assay                   : (int*float option)[]
            peptide_abundance_study_variable          : (int*float option)[]
            peptide_abundance_stdev_study_variable    : (int*float option)[]
            peptide_abundance_std_error_study_variable: (int*float option)[]
            labeling                                  : Ontologies.Labeling
        }

    type PSMSection =
        {
            sequence           : string
            PSM_ID             : int
            accession          : string[]
            unique             : int
            database           : string
            database_version   : string
            search_engine      : string
            search_engine_score: float[]
            reliability        : int
            modifications      : string
            retention_time     : float
            charge             : int
            exp_mass_to_charge : float
            calc_mass_to_charge: float
            uri                : string
            spectra_ref        : string
            pre                : string
            post               : string
            start              : int
            ending             : int
        }

    type AlignedComplete =
        {
            ProtInf  : InferredProteinClassItemOut
            Qpsm     : PSMStatisticsResult
            Quant    : Dto.QuantificationResult
            TableSort: TableSort
        }

    let createAlignedComplete protInf qpsm quant tableSort =
        {
            ProtInf   = protInf
            Qpsm      = qpsm
            Quant     = quant
            TableSort = tableSort
        }

    let createSchema (nameAndType: (string*string)[]) =
        nameAndType
        |> Array.map (fun (name,fieldType) -> sprintf "%s=%s" name fieldType)
        |> String.concat ","

    let frameToTableSort (fieldNames: TableSortFieldNames) (frame: Frame<int,string>) =
        frame
        |> Frame.mapRows (fun r os -> 
            createTableSort
                (
                    os.GetAs<string>(fieldNames.Proteingroup)
                    |> String.split ';'
                    |> Seq.map (fun y ->
                        y
                        |> String.replace "\"" ""
                    )
                    |> Seq.sort
                    |> String.concat ";"
                )
                (
                    os.GetAs<string>(fieldNames.Experiment)
                    |> String.replace "\"" ""
                )
                (os.GetAs<float>(fieldNames.DistinctPeptideCount))
                (
                    if fieldNames.Quant_Heavy.IsSome then
                        Some (os.GetAs<float>(fieldNames.Quant_Heavy.Value))
                    else
                        None
                )
                (os.GetAs<float>(fieldNames.Quant_Light))
                (os.GetAs<float>(fieldNames.QValue))
                (os.GetAs<float>(fieldNames.StudySubject))
                (os.TryGetAs<float>(fieldNames.Subject_SEM) |> OptionalValue.asOption)
                (os.TryGetAs<float>(fieldNames.Subject_StDev) |> OptionalValue.asOption)
        )
        |> Series.values
        |> Seq.toArray

    let getFilePaths path identifier =
        IO.Directory.GetFiles (path, identifier)

    let readQpsm (path: string) =
        SeqIO.Seq.fromFileWithCsvSchema<PSMStatisticsResult>(path, '\t', true, schemaMode=SchemaReader.Csv.SchemaModes.Fill)
        |> Seq.toArray

    let readQuant (path: string) =
        SeqIO.Seq.fromFileWithCsvSchema<Dto.QuantificationResult>(path, '\t', false, schemaMode=SchemaReader.Csv.SchemaModes.Fill, skipLines=1)
        |> Seq.toArray

    let readProt (path: string) =
        SeqIO.Seq.fromFileWithCsvSchema<InferredProteinClassItemOut>(path, '\t', true)
        |> Seq.toArray
        |> Array.map (fun x ->
            {
                x with
                    ProteinGroup=
                        x.ProteinGroup
                        |> String.split ';'
                        |> Seq.map (fun y ->
                            y
                            |> String.replace "\"" ""
                        )
                        |> Seq.sort
                        |> String.concat ";"
            }
        )

    let readTab (path: string) (fieldNames: TableSortFieldNames) =
        let schemaTab =
            createSchema [|
                fieldNames.Proteingroup, "string";
                fieldNames.Experiment, "string";
                fieldNames.DistinctPeptideCount, "float";
                if fieldNames.Quant_Heavy.IsSome then
                    fieldNames.Quant_Heavy.Value, "float";
                fieldNames.Quant_Light, "float";
                fieldNames.QValue, "float";
                if fieldNames.Quant_Light <> fieldNames.StudySubject then
                    fieldNames.StudySubject, "float";
                fieldNames.Subject_SEM, "float";
                fieldNames.Subject_StDev, "float"
            |]
        Frame.ReadCsv(path=path ,separators="\t", schema=schemaTab)
        |> frameToTableSort fieldNames
        |> Array.groupBy (fun x -> x.Experiment)

    // TODO: The following three functions are rather slow and can be improved upon
    let findValueNumberedProt (expNames: (string*int)[]) (experimentValueMap: Map<string,float option> ) =
        expNames
        |> Array.map (fun (experiment,number) ->
            let value =
                experimentValueMap
                |> Map.tryFind experiment
                |> fun x -> 
                    match x with
                    | Some x -> x
                    | None -> None
            number, value
        )

    let findValueNumberedPep (expNames: (string*int)[]) (peptides: (Dto.QuantificationResult*string) []) (fieldName: string) =
        let fieldFunc (tableSort: Dto.QuantificationResult) =
            let res = ReflectionHelper.tryGetPropertyValue tableSort fieldName
            match res with
            |None -> failwith (sprintf "Field %s doesn't exist" fieldName)
            |Some x -> x
        expNames
        |> Array.map (fun (experiment,number) ->
            let value =
                peptides
                |> Array.tryFind (fun (pep,exp) -> exp = experiment)
                |> fun x ->
                    match x with
                    |None -> None
                    |Some y ->
                        Some ((fieldFunc (fst y)) :?> float)
            number, value
        )

    let fieldFuncPSM (psm: PSMStatisticsResult) (fieldName: string) =
        let res = ReflectionHelper.tryGetPropertyValue psm fieldName
        match res with
        |None -> failwith (sprintf "Field %s doesn't exist" fieldName)
        |Some x -> x :?> float

    let concatRuns (fillEmpty: string) (data: (int*float option)[]) =
        data
        |> Array.sortBy fst
        |> Array.map (fun (run,value) ->
            match value with
            | None -> fillEmpty
            | Some s ->
                if System.Double.IsNaN s then
                    "null"
                else
                    string s
        )
        |> String.concat "\t"

    let concatRunsForDifferent (fillEmpty: string) (data: (int*float option)[][]) =
        let length =
            data
            |> Array.head
            |> Array.length
        data
        |> Array.map (Array.sortBy fst)
        |> fun x ->
            [|0 .. (length-1)|]
            |> Array.collect (fun i ->
                x
                |> Array.map (fun var ->
                    let run,value = var.[i]
                    match value with
                    | None -> fillEmpty
                    | Some s ->
                        if System.Double.IsNaN s then
                            "null"
                        else
                            string s
                )
            )
            |> String.concat "\t"

    let sortAndPick (sortBy: 'T -> 'Key) (pick: 'T -> 'C) (input: 'T[]) =
        input
        |> Array.sortBy sortBy
        |> Array.map pick

    let formatOneMD (strF: int -> string -> string) (input: string[]) =
        [|
            for i=0 to (input.Length - 1) do
                yield strF (i+1) input.[i]
        |]
        |> String.concat "\n"

    let formatTwoMD (strF: int -> int -> string -> string) (input: string[][])=
        [|
            for i=0 to (input.Length - 1) do
                for j=0 to (input.[i].Length - 1) do
                    yield strF (i+1) (j+1) input.[i].[j]
        |]
        |> String.concat "\n"

    let formatOne (n1: int) (strF: int -> string) =
        [|
            for i=1 to n1 do
                yield strF i
        |]
        |> String.concat "\t"

    let formatOneForThree (n1: int) (strF1: int -> string) (strF2: int -> string) (strF3: int -> string) =
        [|
            for i=1 to n1 do
                yield ([strF1 i; strF2 i ; strF3 i] |> String.concat "\t")
        |]
        |> String.concat "\t"

    let formatTwo (n1: int) (n2: int) (strF: int -> int -> string) =
        [|
            for i=1 to n1 do
                for j=1 to n2 do
                    yield strF i j
        |]
        |> String.concat "\t"

    let matchOption (format: 'a -> string) (sb: Text.StringBuilder) (item: 'a option) =
        match item with
        |Some x -> 
            x
            |> format
            |> fun x -> sb.AppendLine(x)
        |None -> sb

    /// Retrieves the scan time based on the fitted parameter values (HULQ output).
    let getTargetScanTime (qp:Dto.QuantificationResult) = 
        try
        if qp.GlobalMod = 0 then
            qp.Params_Light.[1] 
        else
            qp.Params_Heavy.[1] 
        with
        | _ -> nan

    let groupsWithSameProteinAccession (protGroups: string[][]) =
        protGroups
        |> Array.distinct
        |> Array.groupBy (fun x -> x |> Array.head)
        |> Array.map snd
        |> Array.filter (fun groups ->
            groups.Length > 1
        )
        |> Array.map (fun groups ->
            groups
            |> Array.mapi (fun i group ->
                group, (sprintf "_%i" i)
            )
        )
        |> Array.concat
        |> Map.ofArray