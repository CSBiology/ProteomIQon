namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PeptideSpectrumMatching
open System.Reflection

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv = 
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name)) 
        let results = parser.Parse argv
        let i' = results.GetResult InstrumentOutput
        let o' = results.GetResult OutputDirectory
        let p' = results.GetResult ParamFile
        let d' = results.GetResult PeptideDataBase
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let i = Path.Combine(directory, i')
        let o = Path.Combine(directory, o')
        let p = Path.Combine(directory, p')
        let d = Path.Combine(directory, d')
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "PeptideSpectrumMatching"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Info (sprintf "Peptide data base -d = %s" d)
        logger.Trace (sprintf "CLIArguments: %A" results)
        if File.Exists d then
           logger.Trace (sprintf "Database found at given location (%s)" d)
        else
           failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
        let p = 
            Json.ReadAndDeserialize<Dto.PeptideSpectrumMatchingParams> p
            |> Dto.PeptideSpectrumMatchingParams.toDomain
        if File.Exists i then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Scoring spectra for %s" i)
            scoreSpectra p o d i
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
            |> FSharpAux.PSeq.iter (scoreSpectra p o d)
        else 
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."

        logger.Info "Done"
        0