namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open ConsensusSpectralLibrary
open System.Reflection

module console1 =

    [<EntryPoint>]
    let main argv =
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let results = parser.Parse argv
        let i'  = results.GetResult InstrumentOutput
        let ii' = results.GetResult SpectralLibrary
        let o' = results.GetResult OutputDirectory
        let p' = results.GetResult ParamFile
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let i = Path.Combine(directory, i')
        let ii = Path.Combine(directory, ii')
        let o = Path.Combine(directory, o')
        let p = Path.Combine(directory, p')
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
