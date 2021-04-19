namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open AlignmentBasedQuantification
open System.Reflection
open ProteomIQon.Core.InputPaths

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv =
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult InstrumentOutput    |> getPathRelativeToDir
        let ii = results.GetResult AlignedPeptides    |> getPathRelativeToDir
        let iii = results.GetResult Metrics           |> getPathRelativeToDir
        let iv = results.GetResult QuantifiedPeptides |> getPathRelativeToDir
        let o = results.GetResult OutputDirectory     |> getPathRelativeToDir
        let p = results.GetResult ParamFile           |> getPathRelativeToDir
        let d = results.GetResult PeptideDataBase     |> getPathRelativeToDir
        Logging.generateConfig o
        let logger = Logging.createLogger "AlignmentBasedQuantification"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "AlignedPeptides -ii = %s" ii)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let p =
            Json.ReadAndDeserialize<Dto.AlignmentBasedQuantificationParams> p
            |> Dto.AlignmentBasedQuantificationParams.toDomain
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
            quantifyPeptides p o d i iii iv ii 
        elif Directory.Exists i && Directory.Exists ii then
            logger.Info "multiple files"
            let mzfiles =
                Directory.GetFiles(i,("*.mzlite"))
            let alignmentFiles =
                Directory.GetFiles(ii,("*.align"))
            let metricFiles =
                Directory.GetFiles(iii,("*.alignmetric"))
            let quantFiles =
                Directory.GetFiles(iv,("*.quant"))
            logger.Trace (sprintf "mz files : %A" mzfiles)
            logger.Trace (sprintf "align files : %A" alignmentFiles)
            let mzFilesAndPepFiles =
                mzfiles
                |> Array.choose (fun mzFilePath ->
                    match alignmentFiles |> Array.tryFind (fun alignfilePath -> (Path.GetFileNameWithoutExtension alignfilePath) = (Path.GetFileNameWithoutExtension mzFilePath)) with
                    | Some alignfilePath -> 
                        match metricFiles |> Array.tryFind (fun metricFilePath -> (Path.GetFileNameWithoutExtension metricFilePath) = (Path.GetFileNameWithoutExtension alignfilePath)) with
                        | Some(metricFilePath) -> 
                            match quantFiles |> Array.tryFind (fun quantFilePath -> (Path.GetFileNameWithoutExtension quantFilePath) = (Path.GetFileNameWithoutExtension metricFilePath)) with
                            | Some(quantFilePath) ->
                                Some(mzFilePath,alignfilePath,metricFilePath,quantFilePath)
                            | None -> 
                                logger.Trace (sprintf "no quant file for %s" mzFilePath)
                                None
                        | None ->
                            logger.Trace (sprintf "no alignment metric file for %s" mzFilePath)
                            None
                    | None ->
                        logger.Trace (sprintf "no alignment file for %s" mzFilePath)
                        None
                    )
            if mzfiles.Length <> mzFilesAndPepFiles.Length then
                logger.Info (sprintf "There are %i mzFiles but only %i can be mapped to there corresponding .align, .alignmetric and .quant files" mzfiles.Length mzFilesAndPepFiles.Length)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            mzFilesAndPepFiles
            |> FSharpAux.PSeq.map (fun (mzFilePath,alignfilePath,metricFilePath,quantFilePath) -> quantifyPeptides p o d mzFilePath quantFilePath metricFilePath alignfilePath)
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        else
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."

        logger.Info "Done"
        0