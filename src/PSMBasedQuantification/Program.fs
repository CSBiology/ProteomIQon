namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PSMBasedQuantification
open System.Reflection

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv =
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let results = parser.Parse argv
        let i' = results.GetResult InstrumentOutput
        let ii' = results.GetResult ScoredPSMs
        let p' = results.GetResult ParamFile
        let o' = results.GetResult OutputDirectory
        let d' = results.GetResult PeptideDataBase
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let i = Path.Combine(directory, i')
        let ii = Path.Combine(directory, ii')
        let p = Path.Combine(directory, p')
        let o = Path.Combine(directory, o')
        let d = Path.Combine(directory, d')
        Logging.generateConfig o
        let logger = Logging.createLogger "PSMBasedQuantification"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "ScoredPSMsPath -ii = %s" ii)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let p =
            Json.ReadAndDeserialize<Dto.QuantificationParams> p
            |> Dto.QuantificationParams.toDomain
        let dbConnection =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        logger.Trace "Set Index on data base if not present."
        SearchDB'.setIndexOnModSequenceAndGlobalMod dbConnection |> ignore
        logger.Trace "Set Index on data base if not present: finished"
        if File.Exists i then
            logger.Info "single file"
            quantifyPeptides p o dbConnection i ii
        elif Directory.Exists i && Directory.Exists ii then
            logger.Info "multiple files"
            let mzfiles =
                Directory.GetFiles(i,("*.mzlite"))
            let pepfiles =
                Directory.GetFiles(ii,("*.qpsm"))
            logger.Trace (sprintf "mz files : %A" mzfiles)
            logger.Trace (sprintf "PEP files : %A" pepfiles)
            let mzFilesAndPepFiles =
                mzfiles
                |> Array.choose (fun mzFilePath ->
                                match pepfiles  |> Array.tryFind (fun pepFilePath -> (Path.GetFileNameWithoutExtension pepFilePath) = (Path.GetFileNameWithoutExtension mzFilePath)) with
                                | Some pepFile -> Some(mzFilePath,pepFile)
                                | None ->
                                    logger.Trace (sprintf "no qpsmFileFor %s" mzFilePath)
                                    None

                              )
            if mzfiles.Length <> mzFilesAndPepFiles.Length then
                logger.Info (sprintf "There are %i mzFiles but %i files containing scored PSMs" mzfiles.Length pepfiles.Length)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            mzFilesAndPepFiles
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> FSharpAux.PSeq.iter (fun (i,ii) -> quantifyPeptides p o dbConnection i ii)
        else
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        dbConnection.Dispose()

        logger.Info "Done"
        0