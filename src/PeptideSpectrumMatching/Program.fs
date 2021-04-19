namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PeptideSpectrumMatching
open System.Reflection
open System
open ProteomIQon.Core.InputPaths

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv = 
        printfn "%A" argv
        printfn "ServerGC: %A" System.Runtime.GCSettings.IsServerGC
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult InstrumentOutput |> getPathRelativeToDir
        let o = results.GetResult OutputDirectory  |> getPathRelativeToDir
        let p = results.GetResult ParamFile        |> getPathRelativeToDir
        let d = results.GetResult PeptideDataBase  |> getPathRelativeToDir
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "PeptideSpectrumMatching"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Info (sprintf "Peptide data base -d = %s" d)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let dbConnection =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        let p = 
            Json.ReadAndDeserialize<Dto.PeptideSpectrumMatchingParams> p
            |> Dto.PeptideSpectrumMatchingParams.toDomain
        if File.Exists i then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Scoring spectra for %s" i)
            scoreSpectra p o dbConnection i
        elif Directory.Exists i then 
            logger.Info (sprintf "multiple files")
            let files = 
                Directory.GetFiles(i,("*.mzlite"))
            logger.Trace (sprintf "Scoring multiple files: %A" files)
            let c = 
                match results.TryGetResult Parallelism_Level with 
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            files 
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> FSharpAux.PSeq.iter (scoreSpectra p o dbConnection)
        else 
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        dbConnection.Dispose()
        logger.Info "Done"
        0