namespace ProteomIQon

open System
open MzIO
open FSharp.Stats
open FSharpAux
open FSharpAux.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute
open FSharp.Plotly

module MzTABAux =

    type TableSort =
        {
            [<FieldAttribute("\"Key #1\"")>]
            Protein: string
            [<FieldAttribute("\"Key #2\"")>]
            Experiment: string
            DistinctPeptideCount: float
            [<FieldAttribute(3)>]
            Quant_Heavy: float
            [<FieldAttribute(4)>]
            Quant_Light: float
            QValue: float
            Ratio: float
            Ratio_CV: float
            Ratio_SEM: float
            Ratio_StDev: float
        }

    type InferredProteinClassItemOut =
        {
            GroupOfProteinIDs: string
            PeptideSequence  : string
            Class            : string
            TargetScore      : float
            DecoyScore       : float
            QValue           : float
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
                    GroupOfProteinIDs =
                        x.GroupOfProteinIDs
                        |> String.split ';'
                        |> Seq.map (fun y ->
                            y
                            |> String.replace "\"" ""
                        )
                        |> Seq.sort
                        |> String.concat ";"
            }
        )

    let readTab (path: string) =
        SeqIO.Seq.fromFileWithCsvSchema<TableSort>(path, '\t', true)
        |> Seq.toArray
        |> Array.map (fun x ->
            {
                x with
                    Protein =
                        x.Protein
                        |> String.split ';'
                        |> Seq.map (fun y ->
                            y
                            |> String.replace "\"" ""
                        )
                        |> Seq.sort
                        |> String.concat ";"
                    Experiment =
                        x.Experiment
                        |> String.replace "\"" ""
            }
        )
        |> Array.groupBy (fun x -> x.Experiment)

    let findValueNumberedProt (expNames: (string*int)[]) (proteins: TableSort []) (fieldName: string) =
        let fieldFunc (tableSort: TableSort) =
            let res = ReflectionHelper.tryGetPropertyValue tableSort fieldName
            match res with
            |None -> failwith (sprintf "Field %s doesn't exist" fieldName)
            |Some x -> x
        expNames
        |> Array.map (fun (experiment,number) ->
            let value =
                proteins
                |> Array.tryFind (fun prot -> prot.Experiment = experiment)
                |> fun x ->
                    match x with
                    |None -> None
                    |Some y ->
                        Some ((fieldFunc y) :?> float)
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