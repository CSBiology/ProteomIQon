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

    type InferenceTask =
        {
            InputPath  : string
            OutputPath : string
        }

        static member  createInferenceTask inp out =
            {
                InputPath  = inp
                OutputPath = out
            }

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
        printfn "\tregexPattern: %s" regexPattern

        printfn "\tReading fastA Sequence"
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

        printfn "\tReading gff3 file"
        /// Create proteinModelInfos: Group all genevariants (RNA) for the gene loci (gene)
        ///
        /// The proteinModelInfos
        let proteinModelInfos =
           try
                GFF3.fromFileWithoutFasta gff3Path
                |> assignTranscriptsToGenes regexPattern
            with
            | err ->
                printfn "\nERROR: Could not read gff3 file %s" gff3Path
                failwithf "\t%s" err.Message
        printfn "\tAssign FastaSeqs to ProteinModelInfo"
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
        //printfn "\tcount fastaSeqs: %i, proteinModels: %i" (Seq.length fastASequences) (Seq.length proteinModels)
        printfn "\tBuild Classification Map"
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

    // TO-DO: is this reversion still needed?
    type SeqConverter() =
        inherit ConverterAttribute()
        override this.convertToObj =
            Converter.Single(fun str ->
                str
                |> String.rev
                |> fun x -> box x)

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

    let readAndInferFile classItemCollection protein peptide groupFiles (tasks:InferenceTask list) =
        let timer = System.Diagnostics.Stopwatch()
        printfn "\tMap peptide Sequences to Proteins"
        let classifiedProteins =
            tasks
            |> List.map (fun task ->
                try
                    let out =
                        task.OutputPath
                    let inp = task.InputPath
                    let psmInputs =
                        Seq.fromFileWithCsvSchema<PSMInput>(inp, '\t', true,schemaMode = SchemaModes.Fill)
                        |> Seq.toList
                    inp,
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
                    printfn "Could not map sequences of file %s to proteins:" task.InputPath
                    failwithf "%s" err.Message
                )

        if groupFiles then
            printfn "\tcreate Combined List"
            let combinedClasses =
                List.collect (fun (_,pepSeq,_) -> pepSeq) classifiedProteins
                |> BioFSharp.Mz.ProteinInference.inferSequences protein peptide
            printfn "\finished Creating combined list"
            classifiedProteins
            |> List.iter (fun (inp,prots,out) ->
                let pepSeqSet = prots |> List.map (fun x -> x.PeptideSequence) |> Set.ofList
                printfn "start with %s" inp
                timer.Start()
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
                printfn "wrote to %s" out
                printfn "Elapsed Time: %A" timer.Elapsed
                timer.Reset()
            )
        else
            classifiedProteins
            |> List.iter (fun (inp,sequences,out) ->
                printfn "\t\t start inferring %s" inp
                BioFSharp.Mz.ProteinInference.inferSequences protein peptide sequences
                |> Seq.map (fun x -> OutputProtein.createOutputProtein x.GroupOfProteinIDs x.Class (proteinGroupToString x.PeptideSequence))
                |> FSharpAux.IO.SeqIO.Seq.toCSV "\t" true
                |> Seq.write out
                printfn "\t\t finished inferring %s" inp)

    let initInferenceFile outdirectory rawFilePath =
        printfn "\t\tSearching for %s" rawFilePath
        if System.IO.File.Exists rawFilePath then
            if System.IO.FileInfo(rawFilePath).Extension = ".qpsm" then
                printfn "\t\t\tIdentified Input as qpsm file"
                let name = System.IO.Path.GetFileNameWithoutExtension rawFilePath
                let outFileName =
                    outdirectory + @"\" + name + ".prot"
                printfn "\t\t\tWill be written to %s" outFileName
                Some (InferenceTask.createInferenceTask rawFilePath outFileName)
            else None
        else
            None
            //failwith (sprintf "File %s could not be found" rawFilePath)

    //let inferProteins gff3Location fastaLocation regexPattern protein peptide groupFiles outDirectory (files:System.IO.FileInfo list) =
    let inferProteins gff3Location fastaLocation (proteinInferenceParams: ProteinInferenceParams) outDirectory (files:System.IO.FileInfo list) =
        printfn "Gathering Tasks"
        let inferenceTasks =
            files
            |> List.collect (fun fileInfo ->
                match fileInfo.Name with
                | "InputFile" ->
                    match initInferenceFile outDirectory fileInfo.FullName with
                    | Some task ->
                        Some task |> List.singleton
                    | None ->
                        failwithf "file %s not found or it does not have the proper format" fileInfo.FullName
                | "InputDirectory" ->
                    printfn "\tBrowse directory %s" fileInfo.FullName
                    let outDirectory =
                        outDirectory + @"\" + (System.IO.FileInfo(fileInfo.FullName).Name)
                    if System.IO.Directory.Exists outDirectory |> not then
                        System.IO.Directory.CreateDirectory outDirectory |> ignore
                    let files = (System.IO.Directory.GetFiles(fileInfo.FullName) |> Array.toList)
                    let tasks =
                        files
                        |> List.map (initInferenceFile outDirectory)
                    printfn "\t For directory %s, could find %i files with proper format" fileInfo.FullName (List.length tasks)
                    tasks
                | name -> failwithf "Input file name %s is not accepted. Either declare input as `InputFile` if it's a single file or `InputDirectory` if it's a directory" name)
            |> List.choose id
        printfn "Finished Gathering tasks"
        printfn ""
        inferenceTasks
        |> List.iter (fun task ->
            printfn "Tasks"
            printfn "\tin: %s" task.InputPath
            printfn "\tout: %s" task.OutputPath
            )
        printfn ""
        let timer = System.Diagnostics.Stopwatch()
        timer.Start()
        printfn "Start building ClassItemCollection"
        let classItemCollection = createClassItemCollection gff3Location fastaLocation proteinInferenceParams.ProteinIdentifierRegex
        printfn "Finished building ClassItemCollection"
        printfn "Time elapsed: %A" timer.Elapsed
        printfn ""
        timer.Stop();timer.Reset();timer.Start()
        printfn "Classify and Infer Proteins"
        readAndInferFile classItemCollection proteinInferenceParams.Protein proteinInferenceParams.Peptide proteinInferenceParams.GroupFiles inferenceTasks
        printfn "Finished Protein Inference"
        printfn "Time elapsed: %A" timer.Elapsed