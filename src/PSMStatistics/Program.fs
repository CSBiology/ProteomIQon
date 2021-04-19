namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PSMStatistics
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
        let i = results.GetResult PSMs            |> getPathRelativeToDir
        let o = results.GetResult OutputDirectory |> getPathRelativeToDir
        let p = results.GetResult ParamFile       |> getPathRelativeToDir
        let d = results.GetResult PeptideDataBase |> getPathRelativeToDir
        Logging.generateConfig o
        let logger = Logging.createLogger "PSMStatistics"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Info (sprintf "Peptide data base -d = %s" d)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let p = 
            Json.ReadAndDeserialize<Dto.PSMStatisticsParams> p
            |> Dto.PSMStatisticsParams.toDomain
        let dbConnection =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        if File.Exists i then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Processing %s" i)
            psmStats p o dbConnection i
        elif Directory.Exists i then 
            printfn "multiple files"
            let files = 
                Directory.GetFiles(i,("*.psm"))
            logger.Trace (sprintf "Processing multiple files: %A" files)
            let c = 
                match results.TryGetResult Parallelism_Level with 
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            files 
            |> FSharpAux.PSeq.map (psmStats p o dbConnection) 
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        else 
            failwith "The given path to the PSMs is neither a valid file path nor a valid directory path."
        dbConnection.Dispose()
        logger.Info "Done"
        0
