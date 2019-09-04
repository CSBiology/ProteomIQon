namespace ProteomIQon

open System.Data.SQLite
open BioFSharp
open BioFSharp.Mz.SearchDB
open Domain
open Core 
open Logary
open System.IO
open BioFSharp.Mz
open MzLite 

open Core.MzLite.Reader
open Core.MzLite.Peaks
open Core.MzLite
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
            Peptide                      = flankedPepSequence
            Protein                      = proteinNames
        }
             
    ///
    let pepValueCalcAndProteinInference (processParams:PSMStatisticsParams) (outputDir:string) (cn:SQLiteConnection) (psms:string) =
        printfn "Now performing peptide spectrum matching: %s Results will be written to: %s" psms outputDir

        let percolatorInFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension psms).Replace(" ","") + ".pin"
            Path.Combine [|outputDir;fileName|]

        let percolatorOutFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension percolatorInFilePath) + ".spsm"
            Path.Combine [|outputDir;fileName|]

        let outFilePath = 
            let fileName = (Path.GetFileNameWithoutExtension psms) + ".qpsm"
            Path.Combine [|outputDir;fileName|]

        printfn "outFilePath:%s" outFilePath

        printfn "Copy peptide DB into Memory"
        let memoryDB = SearchDB.copyDBIntoMemory cn 
        let pepDBTr = memoryDB.BeginTransaction()
        printfn "Copy peptide DB into Memory: finished"
        
        printfn "Read scored PSMs."
        let psms =
            FSharpAux.IO.SchemaReader.Csv.CsvReader<Dto.PeptideSpectrumMatchingResult>().ReadFile(psms,'\t',false,0)
            |> Array.ofSeq
        printfn "Read scored PSMs: finished"

        printfn "Prepare processing functions."
        let maxCharge = psms |> Array.map (fun x -> x.Charge) |> Array.max
        let proteinAndClvIdxLookUp = initProteinAndClvIdxLookUp memoryDB pepDBTr
        let toPercolatorIn = initToPercolatorIn maxCharge processParams.FastaHeaderToName proteinAndClvIdxLookUp
        printfn "Finished preparing processing functions."
        
        printfn "Converting psms to percolatorIn format."
        let percolatorIn =
            psms
            |> Array.map toPercolatorIn
            |> Array.filter (fun x -> nan.Equals(x.SequestNormDeltaNext) = false && nan.Equals(x.AndroNormDeltaNext) = false )

        let pepSequenceIDToMissCleavagesAndProt =
            percolatorIn
            |> Array.filter (fun candidatePSM -> candidatePSM.Label = 1)
            |> Array.map (fun candidatePSM -> candidatePSM.PSMId,(candidatePSM.Protein,candidatePSM.MissCleavages))
            |> Map.ofArray
        printfn "Converting psms to percolatorIn format: finished"
        
        printfn "Writing percolatorIn.tab. to disk"
        percolatorIn
        |> FSharpAux.IO.SeqIO.Seq.toCSV "\t" true
        |> Seq.map (fun x -> FSharpAux.String.replace ";" "\t" x)
        |> FSharpAux.IO.FileIO.writeToFile false percolatorInFilePath
        printfn "Writing percolatorIn.tab. to disk:finished"
 

        printfn "Executing Percolator"
        let percolatorParams =
            [
            PercolatorParams.GeneralOptions  [(GeneralOptions.PostProcessing_TargetDecoyCompetition)]
            PercolatorParams.FileInputOptions [(FileInputOptions.PINTAB (System.IO.FileInfo(percolatorInFilePath)))]
            PercolatorParams.FileOutputOptions [(FileOutputOptions.POUTTAB_PSMs (System.IO.FileInfo(percolatorOutFilePath)));];
            ]
        let executePercolator =
            let percPath = 
                let assembly = Assembly.GetExecutingAssembly()
                System.IO.FileInfo(assembly.Location).DirectoryName + @"\percolator-v3-01\bin\percolator.exe"
            printfn "\tlooking for percolator at %s" percPath
            let percolator = new PercolatorWrapper(OperatingSystem.Windows,percPath)
            percolator.Percolate percolatorParams
        printfn "Executing Percolator:finished"

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
                                                Some {PSMId=psmID;GlobalMod=candidatePSM.GlobalMod;PepSequenceID=candidatePSM.PepSequenceID;ModSequenceID=candidatePSM.ModSequenceID;Label=candidatePSM.Label;ScanNr=candidatePSM.ScanNr;Charge=candidatePSM.Charge;PrecursorMZ=candidatePSM.PrecursorMZ;TheoMass=candidatePSM.TheoMass;AbsDeltaMass=candidatePSM.AbsDeltaMass;PeptideLength=candidatePSM.PeptideLength;MissCleavages=missCleavages;SequestScore=candidatePSM.SequestScore;SequestNormDeltaBestToRest=candidatePSM.SequestNormDeltaBestToRest;SequestNormDeltaNext=candidatePSM.SequestNormDeltaNext;AndroScore=candidatePSM.AndroScore;AndroNormDeltaBestToRest=candidatePSM.AndroNormDeltaBestToRest;AndroNormDeltaNext=candidatePSM.AndroNormDeltaNext;PercolatorScore=validPSM.PercolatorScore;QValue=validPSM.QValue;PEPValue=validPSM.PosteriorErrorProbability;StringSequence=candidatePSM.StringSequence;ProteinNames= proteins }
                                        | None          -> None
                                    | _ -> None
                                )

            result
            |> FSharpAux.IO.SeqIO.Seq.toCSV "\t" true
            |> FSharpAux.IO.FileIO.writeToFile false outFilePath

        with
        | ex ->
            printfn "%A" ex.Message
            ()
        printfn "Done."
         