namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open BioFSharp.Mz
open SpectralLibrary
open System.Reflection
open ProteomIQon.Core.InputPaths

module console1 =

    [<EntryPoint>]
    let main argv =
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult InstrumentOutput     |> getPathRelativeToDir
        let d = results.GetResult PeptideDataBase      |> getPathRelativeToDir
        let ii = results.GetResult PSMStatisticsResult |> getPathRelativeToDir
        let iii = results.GetResult QuantResult        |> getPathRelativeToDir
        let o = results.GetResult OutputDirectory      |> getPathRelativeToDir
        let p = results.GetResult ParamFile            |> getPathRelativeToDir
        Logging.generateConfig o
        let logger = Logging.createLogger "SpectralLibrary"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "Peptide data base -d = %s" d)
        logger.Info (sprintf "PSMStatisticsResult -ii = %s" ii)
        logger.Info (sprintf "QuantResult -iii = %s" iii)
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
        if File.Exists i && File.Exists ii then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Scoring spectra for %s and %s" i ii)
            createSpectralLibrary o p dbConnection  (i, ii, iii)
        elif Directory.Exists i && Directory.Exists ii then
            logger.Info (sprintf "multiple files")
            let instrumentFiles =
                Directory.GetFiles(i,("*.mzlite"))
            let resultFiles =
                Directory.GetFiles(ii,("*.qpsm"))
            let quantFiles =
                Directory.GetFiles(iii,("*.quant"))
            let matchedFiles =
                let instrumentPsm =
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
                instrumentPsm
                |> Array.collect (fun (instrument,psm) ->
                    quantFiles
                    |> Array.choose (fun quantFiles ->
                        if Path.GetFileNameWithoutExtension instrument = Path.GetFileNameWithoutExtension quantFiles then
                            Some (instrument, psm, quantFiles)
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
        elif (Directory.Exists i && File.Exists ii) || (File.Exists i && Directory.Exists ii) then
            failwith "Both, instrument output and PSMStatistics result must be either a path or a folder."
        else
            failwith "The given paths to the instrument output and PSMStatistics result are neither valid file paths nor valid directory paths."
        logger.Info "Done"
        0