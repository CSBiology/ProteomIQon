namespace ProteomIQon

open BioFSharp
open BioFSharp.Mz.ProteinInference
open BioFSharp.IO
open FSharpAux
open System.Text.RegularExpressions
open PeptideClassification
open FSharpAux.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Csv
open FSharpAux.IO.SchemaReader.Attribute
open GFF3
open Domain

module ProteinInference =

    /// This type represents one element of the final output. It's the protein, with all the peptides that were measured and are pointing to it.
    type OutputProtein =
        {
            ProteinID        : string
            EvidenceClass    : PeptideEvidenceClass
            PeptideSequences : string
        }

    /// Packages info into one element of the final output. It's the protein, with all the peptides that were measured and are pointing to it.
        static member createOutputProtein protID evidenceClass sequences =
            {
                ProteinID        = protID
                EvidenceClass    = evidenceClass
                PeptideSequences = sequences
            }

    /// Represents one peptide-entry with its evidence class and the proteins it points to.
    type ClassInfo =
        {
            Sequence : string
            Class    : PeptideClassification.PeptideEvidenceClass
            Proteins : string []
        }

    ///checks if GFF line describes gene
    let isGene (item: GFFLine<seq<char>>) =
        match item with
        | GFFEntryLine x -> x.Feature = "gene"
        | _ -> false

    ///checks if GFF line describes rna
    let isRNA (item: GFFLine<seq<char>>) =
        match item with
        | GFFEntryLine x -> if x.Feature = "mRNA" then Some x else None
        | _ -> None

    /// Reads geographical information about protein from gff entry and builds the modelinfo of it
    /// This function takes an RNA gff3 entry and therefore will contain the splice variant id of the gene in its result.
    /// This splice variant id should be the same in the given FastA-file.
    let createProteinModelInfoFromEntry i locus (entry:GFFEntry) =
        let attributes = entry.Attributes
        /// Same as in FastA file
        let spliceVariantID =
            match Map.tryFind "Name" attributes with
            | Some res ->
                res.Head
            | None ->
                failwithf "could not find spliceVariantId for locus %s" locus
        let chromosomeID = entry.Seqid

        let direction =
            match entry.Strand with
            |'+' -> PeptideClassification.StrandDirection.Forward
            |'-' -> PeptideClassification.StrandDirection.Reverse

        PeptideClassification.createProteinModelInfo spliceVariantID chromosomeID direction locus i Seq.empty Seq.empty

    /// By reading GFF creates the protein models (relationships of proteins to each other) which basically means grouping the rnas over the gene loci
    /// TODO: Don't group over order but rather group over id
    let assignTranscriptsToGenes regexPattern (gffLines: seq<GFF3.GFFLine<seq<char>>>)  =
        gffLines
        // transcripts are grouped by the gene they originate from
        |> Seq.groupWhen isGene
        |> Seq.map (fun group ->
            match Seq.head group with
            | GFFEntryLine x ->
                let locus = x.Attributes.["ID"].Head // =genename, this value is used to assign mRNAs of the same gene together
                group
                |> Seq.choose isRNA //the transcripts of the gene are chosen
                |> Seq.mapi (fun i element ->
                    // every transcript of gene gets its own number i and other info is collected from element and used for info of protein
                    let modelInfo = createProteinModelInfoFromEntry i locus element
                    let r = System.Text.RegularExpressions.Regex.Match(modelInfo.Id,regexPattern)
                    // the gff3 id has to be matched with the sequence in the fasta file. therefore the regexpattern is used
                    if r.Success then
                        r.Value,
                        modelInfo
                    else
                        failwithf "could not match gff3 entry id %s with regexpattern %s. Either gff3 file is corrupt or regexpattern is not correct" modelInfo.Id regexPattern
                    )

            | _ -> Seq.empty
            )
        |> Seq.concat
        |> Map.ofSeq

    /// Given a ggf3 and a fasta file, creates a collection of all theoretically possible peptides and the proteins they might
    /// originate from
    ///
    /// ProteinClassItem: For a single peptide Sequence, contains information about all proteins it might originate from and
    /// its evidence class.
    ///
    /// No experimental data
    let createClassItemCollection gff3Path fastAPath regexPattern =

        let fastASequences =
            try
                //fileDir + "Chlamy_Cp.fastA"
                fastAPath
                |> FastA.fromFile BioArray.ofAminoAcidString
                |> Seq.map (fun item ->
                    item.Header
                    ,item.Sequence)
            with
            | err ->
                printfn "Could not read FastA file %s"fastAPath
                failwithf "%s" err.Message

        /// Create proteinModelInfos: Group all genevariants (RNA) for the gene loci (gene)
        ///
        /// The proteinModelInfos
        let proteinModelInfos =
           try
                GFF3.fromFileWithoutFasta gff3Path
                |> assignTranscriptsToGenes regexPattern
            with
            | err ->
                printfn "ERROR: Could not read gff3 file %s" gff3Path
                failwithf "%s" err.Message
        //reads from file to an array of FastaItems.

        /// Assigned fasta sequences to model Infos
        let proteinModels =
            try
                fastASequences
                |> Seq.mapi (fun i (header,sequence) ->
                    let regMatch = System.Text.RegularExpressions.Regex.Match(header,regexPattern)
                    if regMatch.Success then
                        match Map.tryFind regMatch.Value proteinModelInfos with
                        | Some modelInfo ->
                            Some (createProteinModel modelInfo sequence)
                        | None ->
                            failwithf "Could not find protein with id %s in gff3 File" regMatch.Value
                    else
                        printfn "Could not extract id of header \"%s\" with regexPattern %s" header regexPattern
                        createProteinModelInfo header "Placeholder" StrandDirection.Forward (sprintf "Placeholder%i" i) 0 Seq.empty Seq.empty
                        |> fun x -> createProteinModel x sequence
                        |> Some
                    )
            with
            | err ->
                printfn "Could not assign FastA sequences to RNAs"
                failwithf "%s" err.Message
        try
            let ppRelationModel =
                let digest sequence =
                    Digestion.BioArray.digest (Digestion.Table.getProteaseBy "Trypsin") 0 sequence
                    |> Digestion.BioArray.concernMissCleavages 0 3
                    |> Seq.map (fun p -> p.PepSequence |> List.toArray) // TODO not |> List.toArray

                PeptideClassification.createPeptideProteinRelation digest proteinModels

            let spliceVariantCount = PeptideClassification.createLocusSpliceVariantCount ppRelationModel

            let classified =
                ppRelationModel.GetArrayOfKeys
                |> Array.map (fun peptide ->
                    let proteinInfo = (ppRelationModel.TryGetByKey peptide).Value
                    let proteinIds = Seq.map (fun (x : PeptideClassification.ProteinModelInfo<string,string,string>) -> x.Id) proteinInfo
                    let (c,x) = PeptideClassification.classify spliceVariantCount (peptide, (ppRelationModel.TryGetByKey peptide).Value)
                    {
                        Sequence = BioArray.toString x;
                        Class = c;
                        Proteins = (proteinIds |> Seq.toArray)
                    }
                    )

            classified |> Array.map (fun ci -> ci.Sequence,(createProteinClassItem ci.Proteins ci.Class ci.Sequence)) |> Map.ofArray
        with
        | err ->
            printfn "\nERROR: Could not build classification map"
            failwithf "\t%s" err.Message

    type PSMInput =
        {
            //TO-DO: is this reversion still needed?
            [<FieldAttribute("StringSequence")>]//[<SeqConverter>]
            Seq:string
        }

    let removeModification pepSeq =
        String.filter (fun c -> System.Char.IsLower c |> not && c <> '[' && c <> ']') pepSeq

    let proteinGroupToString (proteinGroup:string[]) =
        Array.reduce (fun x y ->  x + ";" + y) proteinGroup

    let readAndInferFile classItemCollection protein peptide groupFiles outdirectory rawFolderPath=
        let rawFilePaths = System.IO.Directory.GetFiles (rawFolderPath, "*.qpsm")
                           |> Array.toList
        let classifiedProteins =
            rawFilePaths
            |> List.map (fun filePath ->
                try
                    let out =
                        let foldername = (rawFolderPath.Split ([|"\\"|], System.StringSplitOptions.None))
                        outdirectory + @"\" + foldername.[foldername.Length - 1] + "\\" + (System.IO.Path.GetFileNameWithoutExtension filePath) + ".prot"
                    let psmInputs =
                        Seq.fromFileWithCsvSchema<PSMInput>(filePath, '\t', true,schemaMode = SchemaModes.Fill)
                        |> Seq.toList
                    filePath,
                    psmInputs
                    |> List.map (fun s ->
                        let s =
                            s.Seq

                        match Map.tryFind (removeModification s) classItemCollection with
                        | Some (x:ProteinClassItem<'sequence>) -> createProteinClassItem x.GroupOfProteinIDs x.Class s
                        | None ->
                            failwithf "Could not find sequence %s in classItemCollection"  s
                        ),
                        out
                with
                | err ->
                    printfn "Could not map sequences of file %s to proteins:" filePath
                    failwithf "%s" err.Message
                )

        if groupFiles then
            let combinedClasses =
                List.collect (fun (_,pepSeq,_) -> pepSeq) classifiedProteins
                |> BioFSharp.Mz.ProteinInference.inferSequences protein peptide
            classifiedProteins
            |> List.iter (fun (inp,prots,out) ->
                let pepSeqSet = prots |> List.map (fun x -> x.PeptideSequence) |> Set.ofList
                combinedClasses
                |> Seq.choose (fun ic ->
                    let filteredPepSet =
                        ic.PeptideSequence
                        |> Array.filter (fun pep -> Set.contains pep pepSeqSet)
                    if filteredPepSet = [||] then
                        None
                    else
                        Some (ic.Class,ic.GroupOfProteinIDs,proteinGroupToString filteredPepSet)
                    )
                |> Seq.map (fun (c,prots,pep) -> OutputProtein.createOutputProtein prots c pep)
                |> FSharpAux.IO.SeqIO.Seq.toCSV "\t" true
                |> Seq.write out
            )
        else
            classifiedProteins
            |> List.iter (fun (inp,sequences,out) ->
                printfn "\t\t start inferring %s" inp
                BioFSharp.Mz.ProteinInference.inferSequences protein peptide sequences
                |> Seq.map (fun x -> OutputProtein.createOutputProtein x.GroupOfProteinIDs x.Class (proteinGroupToString x.PeptideSequence))
                |> FSharpAux.IO.SeqIO.Seq.toCSV "\t" true
                |> Seq.write out
            )

    let inferProteins gff3Location fastaLocation (proteinInferenceParams: ProteinInferenceParams) outDirectory rawFolderPath =
        printfn "Start building ClassItemCollection"
        let classItemCollection = createClassItemCollection gff3Location fastaLocation proteinInferenceParams.ProteinIdentifierRegex
        printfn "Classify and Infer Proteins"
        readAndInferFile classItemCollection proteinInferenceParams.Protein proteinInferenceParams.Peptide proteinInferenceParams.GroupFiles outDirectory rawFolderPath