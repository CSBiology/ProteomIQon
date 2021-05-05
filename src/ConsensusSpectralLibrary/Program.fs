namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
open ConsensusSpectralLibrary
open System.Reflection
open ProteomIQon.Core.InputPaths

module console1 =

    [<EntryPoint>]
    let main argv =
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i  = results.GetResult InstrumentOutput |> getPathRelativeToDir
        let ii = results.GetResult SpectralLibrary  |> getPathRelativeToDir
        let o = results.GetResult OutputDirectory   |> getPathRelativeToDir
        let p = results.GetResult ParamFile         |> getPathRelativeToDir
        Logging.generateConfig o
        let logger = Logging.createLogger "SpectralLibrary"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "SpectralLibrary -ii = %s" ii)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "InputParameterPath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let p =
            Json.ReadAndDeserialize<Dto.ConsensusSpectralLibraryParams> p
            |> Dto.ConsensusSpectralLibraryParams.toDomain
        Directory.CreateDirectory(o) |> ignore
        if File.Exists i then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Generating consensus library for %s.\nIt is recommended to use at least two libraries for a consensus library" i)
            createSpectralLibrary logger p o ii i 
        elif Directory.Exists i then
            logger.Info (sprintf "multiple files")
            let swathFiles =
                Directory.GetFiles(i,("*.mzlite"))
            logger.Trace (sprintf "Generating consensus library for: %A" swathFiles)
            swathFiles
            |> Array.iter (fun i -> createSpectralLibrary logger p o ii i) 
        else
            failwith "The given paths to the instrument output and PSMStatistics result are neither valid file paths nor valid directory paths."
        logger.Info "Done"
        0
