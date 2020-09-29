namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open BioFSharp.Mz
open SpectralLibrary

module console1 =

    [<EntryPoint>]
    let main argv =
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let results = parser.Parse argv
        let i = results.GetResult InstrumentOutput
        let d = results.GetResult PeptideDataBase
        let q = results.GetResult PSMStatisticsResult
        let o = results.GetResult OutputDirectory
        let p = results.GetResult ParamFile
        Logging.generateConfig o
        let logger = Logging.createLogger "TableSort"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "Peptide data base -d = %s" d)
        logger.Info (sprintf "PSMStatisticsResult -g = %s" q)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "InputParameterPath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let dbConnection =
            if File.Exists d then
                logger.Trace (sprintf "Database found at given location (%s)" d)
                SearchDB.getDBConnection d
            else
                failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        let p =
            Json.ReadAndDeserialize<Dto.SpectralLibraryParams> p
            |> Dto.SpectralLibraryParams.toDomain
        if File.Exists i && File.Exists q then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Scoring spectra for %s and %s" i q)
            createSpectralLibrary o p dbConnection (i,q)
        elif Directory.Exists i && Directory.Exists q then
            logger.Info (sprintf "multiple files")
            let instrumentFiles =
                Directory.GetFiles(i,("*.mzlite"))
            let resultFiles =
                Directory.GetFiles(q,("*.qpsm"))
            let matchedFiles =
                instrumentFiles
                |> Array.collect (fun instrumentFile ->
                    resultFiles
                    |> Array.choose (fun resultFile ->
                        if Path.GetFileNameWithoutExtension instrumentFile = Path.GetFileNameWithoutExtension resultFile then
                            Some (instrumentFile, resultFile)
                        else
                            None
                    )
                )
            logger.Trace (sprintf "Scoring multiple files: %A" matchedFiles)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            matchedFiles
            |> FSharpAux.PSeq.map (createSpectralLibrary o p dbConnection)
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        elif (Directory.Exists i && File.Exists q) || (File.Exists i && Directory.Exists q) then
            failwith "Both, instrument output and PSMStatistics result must be either a path or a folder."
        else
            failwith "The given paths to the instrument output and PSMStatistics result are neither valid file paths nor valid directory paths."
        logger.Info "Done"
        0