namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open PSMStatistics
open System.Reflection

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv = 
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name)) 
        let results = parser.Parse argv
        let i' = results.GetResult PSMs
        let o' = results.GetResult OutputDirectory
        let p' = results.GetResult ParamFile
        let d' = results.GetResult PeptideDataBase
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let i = Path.Combine(directory, i')
        let p = Path.Combine(directory, p')
        let o = Path.Combine(directory, o')
        let d = Path.Combine(directory, d')
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
        if File.Exists i then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Processing %s" i)
            pepValueCalcAndProteinInference p o d i
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
            |> FSharpAux.PSeq.map (pepValueCalcAndProteinInference p o d) 
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        else 
            failwith "The given path to the PSMs is neither a valid file path nor a valid directory path."

        logger.Info "Done"
        0
