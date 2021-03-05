namespace ProteomIQon

open System.Data.SQLite
open BioFSharp
open BioFSharp.Mz.SearchDB
open Domain
open Core 
open System.IO
open BioFSharp.Mz
open FSharpAux.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Csv
open FSharpAux.IO.SchemaReader.Attribute
open BioFSharp.Mz.PercolatorWrapper
open BioFSharp.Mz.PercolatorWrapper.Parameters
open System.Reflection


module PSMStatistics = 
    open System.IO

    type PercolatorIn = {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
        PSMId                        : string
        Label                        : int
        ScanNr                       : int
        Charge                       : int []
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
        Peptide                      : string
        Protein                      : string 
        }


    type PercolatorPSMOut = {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
        [<FieldAttribute(0)>]
        PSMId                        : string;
        [<FieldAttribute(1)>]
        PercolatorScore              : float
        [<FieldAttribute(2)>]
        QValue                       : float
        [<FieldAttribute(3)>]
        PosteriorErrorProbability    : float
        [<FieldAttribute(4)>]
        StringSequence               : string
        }
    
    /// Returns a list of proteins retrieved by PepsequenceID
    let initProteinAndClvIdxLookUp (cn:SQLiteConnection) tr =
        let selectCleavageIdxByPepSeqID   = Db.SQLiteQuery.prepareSelectCleavageIndexByPepSequenceID cn tr
        let selectProteinByProtID         = Db.SQLiteQuery.prepareSelectProteinByID cn tr
        (fun pepSequenceID ->
                selectCleavageIdxByPepSeqID pepSequenceID
                |> List.map (fun (clid,protID,peps,clvgS,clvgE,missClv) -> selectProteinByProtID protID,clvgS,clvgE,missClv )
        )
    
    /// Retrieves the flanking amino acids of a peptide Sequence originating of a given Protein Sequence
    let getFlankingAminoAcids missClvStart missClvEnd (protSequence:string) =
        let Nterminal = if missClvStart = 0 then "-" else protSequence.[missClvStart-1].ToString()
        let CTerminal = if missClvEnd = protSequence.Length-1 then "-" else  protSequence.[missClvEnd+1].ToString()
        Nterminal,CTerminal

    let restorePSMID psmID =
        FSharpAux.String.split '_' psmID  
        |> Array.item 0
        |> FSharpAux.String.split '-'
        |> String.concat " "

    ///
    let initToPercolatorIn maxCharge fastaHeaderToName proteinAndClvIdxLookUp (psm:Dto.PeptideSpectrumMatchingResult) =
        let charge =
            let ch = psm.Charge
            Array.init maxCharge (fun i -> if i = ch-1 then 1 else 0)
        let lookUp = proteinAndClvIdxLookUp psm.PepSequenceID
        let missCleavages,nTerminal,cTerminal =
            match lookUp with
            | [] -> 0,"",""
            | items ->
                items
                |> List.minBy (fun ((protID,name,seq),clvgS,clvgE,missClv) -> missClv)
                |> fun ((protID,name,seq),clvgS,clvgE,missClv) ->
                    let nTerm,cTerm = getFlankingAminoAcids clvgS clvgE seq
                    if psm.Label = 1 then
                        missClv,nTerm + ".", "." + cTerm
                    else
                        missClv, cTerm + ".", "." + nTerm
        let proteinNames =
            let parseLabel label = if label = 1 then "" else "decoy_"
            match lookUp with
            | [] -> sprintf "%snoProt" (parseLabel psm.Label)
            | items ->
                items
                |> List.map (fun ((protID,name,seq),clvgS,clvgE,missClv) -> sprintf "%s%s" (parseLabel psm.Label) (fastaHeaderToName name))
                |> String.concat ";"
        
        let flankedPepSequence =
            let sequence = if psm.Label = 1 then psm.StringSequence else FSharpAux.String.rev psm.StringSequence
            nTerminal + sequence + cTerminal
        {
            PSMId                        = psm.PSMId
            Label                        = psm.Label
            ScanNr                       = psm.ScanNr
            Charge                       = charge
            PrecursorMZ                  = psm.PrecursorMZ
            TheoMass                     = psm.TheoMass
            AbsDeltaMass                 = psm.AbsDeltaMass
            PeptideLength                = psm.PeptideLength
            MissCleavages                = missCleavages
            SequestScore                 = psm.SequestScore
            SequestNormDeltaBestToRest   = psm.SequestNormDeltaBestToRest
            SequestNormDeltaNext         = psm.SequestNormDeltaNext
            AndroScore                   = psm.AndroScore
            AndroNormDeltaBestToRest     = psm.AndroNormDeltaBestToRest
            AndroNormDeltaNext           = psm. AndroNormDeltaNext
            XtandemScore                 = psm.XtandemScore
            XtandemNormDeltaBestToRest   = psm.XtandemNormDeltaBestToRest
            XtandemNormDeltaNext         = psm.XtandemNormDeltaNext
            Peptide                      = flankedPepSequence
            Protein                      = proteinNames
        }
             
    ///
    let pepValueCalcAndProteinInference (processParams:PSMStatisticsParams) (outputDir:string) (d:string) (psms:string) =

        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension psms)

        logger.Trace (sprintf "Input file: %s" psms)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)

        logger.Trace (sprintf "Now performing peptide spectrum matching: %s Results will be written to: %s" psms outputDir)

        let percolatorInFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension psms).Replace(" ","") + ".pin"
            Path.Combine [|outputDir;fileName|]

        let percolatorOutFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension percolatorInFilePath) + ".spsm"
            Path.Combine [|outputDir;fileName|]

        let percolatorDecoyOutFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension percolatorInFilePath) + "_decoy.spsm"
            Path.Combine [|outputDir;fileName|]

        let outFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension psms) + ".qpsm"
            Path.Combine [|outputDir;fileName|]

        logger.Trace (sprintf "percolatorInFilePath:%s" percolatorInFilePath)
        logger.Trace (sprintf "percolatorOutFilePath:%s" percolatorOutFilePath)
        logger.Trace (sprintf "percolatorDecoyOutFilePath:%s" percolatorDecoyOutFilePath)
        logger.Trace (sprintf "outFilePath:%s" outFilePath)

        logger.Trace "Copy peptide DB into Memory"
        let cn =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        let memoryDB = SearchDB.copyDBIntoMemory cn
        cn.Dispose() 
        let pepDBTr = memoryDB.BeginTransaction()
        logger.Trace "Copy peptide DB into Memory: finished"
        
        logger.Trace "Read scored PSMs."
        let psms =
            FSharpAux.IO.SchemaReader.Csv.CsvReader<Dto.PeptideSpectrumMatchingResult>().ReadFile(psms,'\t',false,0)
            |> Array.ofSeq
        logger.Trace "Read scored PSMs: finished"

        logger.Trace "Prepare processing functions."
        let maxCharge = psms |> Array.map (fun x -> x.Charge) |> Array.max
        let proteinAndClvIdxLookUp = initProteinAndClvIdxLookUp memoryDB pepDBTr
        let toPercolatorIn = initToPercolatorIn maxCharge processParams.FastaHeaderToName proteinAndClvIdxLookUp
        logger.Trace "Finished preparing processing functions."
        
        logger.Trace "Converting psms to percolatorIn format."
        let percolatorIn =
            psms
            |> Array.map toPercolatorIn
            |> Array.filter (fun x -> nan.Equals(x.SequestNormDeltaNext) = false && nan.Equals(x.AndroNormDeltaNext) = false )

        let pepSequenceIDToMissCleavagesAndProt =
            percolatorIn
            |> Array.filter (fun candidatePSM -> candidatePSM.Label = 1)
            |> Array.map (fun candidatePSM -> candidatePSM.PSMId,(candidatePSM.Protein,candidatePSM.MissCleavages))
            |> Map.ofArray
        logger.Trace "Converting psms to percolatorIn format: finished"
        
        logger.Trace "Writing percolatorIn.tab. to disk"
        percolatorIn
        |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
        |> Seq.map (fun x -> FSharpAux.String.replace ";" "\t" x)
        |> FSharpAux.IO.FileIO.writeToFile false percolatorInFilePath
        logger.Trace "Writing percolatorIn.tab. to disk: finished"
 

        logger.Trace "Executing Percolator"
        let percolatorParams =
            [
            PercolatorParams.GeneralOptions  [(GeneralOptions.PostProcessing_TargetDecoyCompetition)]
            PercolatorParams.FileInputOptions [(FileInputOptions.PINTAB (System.IO.FileInfo(percolatorInFilePath)))]
            PercolatorParams.FileOutputOptions [(FileOutputOptions.POUTTAB_PSMs (System.IO.FileInfo(percolatorOutFilePath)));];
            PercolatorParams.FileOutputOptions [(FileOutputOptions.POUTTAB_DecoyPSMs (System.IO.FileInfo(percolatorDecoyOutFilePath)))]
            ]
        let executePercolator =
            let percPath = 
                let assembly = Assembly.GetExecutingAssembly()
                System.IO.FileInfo(assembly.Location).DirectoryName + @"\percolator-v3-01\binaries\percolator.exe"
            logger.Trace (sprintf "\tlooking for percolator at %s" percPath)
            let percolator = new PercolatorWrapper(OperatingSystem.Windows,percPath)
            percolator.Percolate percolatorParams
        logger.Trace "Executing Percolator:finished"

        try
            let scoredPSMs =
                FSharpAux.IO.SchemaReader.Csv.CsvReader<PercolatorPSMOut>(SchemaMode=SchemaModes.Fill,Verbose=false).ReadFile(percolatorOutFilePath,'\t',false,1)
                |> Seq.filter (fun scoredPSM -> scoredPSM.QValue < processParams.QValueThreshold && scoredPSM.PosteriorErrorProbability < processParams.PepValueThreshold)
                |> Seq.map (fun scoredPSM -> scoredPSM .PSMId,scoredPSM )
                |> Map.ofSeq

            let result: Dto.PSMStatisticsResult [] =
                psms
                |> Array.choose (fun candidatePSM ->
                                    match candidatePSM.Label with
                                    | x when x = 1 ->
                                        match Map.tryFind candidatePSM.PSMId scoredPSMs with
                                        | Some validPSM ->
                                            let psmID = restorePSMID validPSM.PSMId 
                                            let proteins,missCleavages = pepSequenceIDToMissCleavagesAndProt.[validPSM.PSMId]
                                            Some {
                                            PSMId                       = psmID
                                            GlobalMod                   = candidatePSM.GlobalMod
                                            PepSequenceID               = candidatePSM.PepSequenceID
                                            ModSequenceID               = candidatePSM.ModSequenceID
                                            Label                       = candidatePSM.Label
                                            ScanNr                      = candidatePSM.ScanNr
                                            ScanTime                    = candidatePSM.ScanTime
                                            Charge                      = candidatePSM.Charge
                                            PrecursorMZ                 = candidatePSM.PrecursorMZ
                                            TheoMass                    = candidatePSM.TheoMass
                                            AbsDeltaMass                = candidatePSM.AbsDeltaMass
                                            PeptideLength               = candidatePSM.PeptideLength
                                            MissCleavages               = missCleavages
                                            SequestScore                = candidatePSM.SequestScore
                                            SequestNormDeltaBestToRest  = candidatePSM.SequestNormDeltaBestToRest
                                            SequestNormDeltaNext        = candidatePSM.SequestNormDeltaNext
                                            AndroScore                  = candidatePSM.AndroScore
                                            AndroNormDeltaBestToRest    = candidatePSM.AndroNormDeltaBestToRest
                                            AndroNormDeltaNext          = candidatePSM.AndroNormDeltaNext
                                            XtandemScore                = candidatePSM.XtandemScore
                                            XtandemNormDeltaBestToRest  = candidatePSM.XtandemNormDeltaBestToRest
                                            XtandemNormDeltaNext        = candidatePSM.XtandemNormDeltaNext
                                            PercolatorScore             = validPSM.PercolatorScore
                                            QValue                      = validPSM.QValue
                                            PEPValue                    = validPSM.PosteriorErrorProbability
                                            StringSequence              = candidatePSM.StringSequence
                                            ProteinNames                = proteins 
                                            }
                                        | None  -> 
                                            None
                                    | _ -> 
                                        None
                                )
            logger.Trace (sprintf "Number of results: %i" result.Length)
            result
            |> FSharpAux.IO.SeqIO.Seq.CSV "\t" true true
            |> FSharpAux.IO.FileIO.writeToFile false outFilePath

        with
        | ex ->
            logger.Trace (sprintf "%A" ex.Message)
            printfn "%A" ex.Message
            ()
        logger.Trace "Done."
         