namespace ProteomIQon

open System
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
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)       
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse(argv)
        //let results = parser.ParseCommandLine(inputs = argv, raiseOnUsage = false)
        let i = results.GetResult PSMs            |> List.map getPathRelativeToDir
        let o = results.GetResult OutputDirectory |> getPathRelativeToDir
        let p = results.GetResult ParamFile       |> getPathRelativeToDir
        let d = results.GetResult PeptideDataBase |> getPathRelativeToDir
        let dc = results.Contains DiagnosticCharts
        Logging.generateConfig o
        let logger = Logging.createLogger "PSMStatistics"
        logger.Info (sprintf "InputFilePath -i = %A" i)
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
        let files = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.psm"))) i
            |> Array.ofSeq
        let c = 
            match results.TryGetResult Parallelism_Level with 
            | Some c    -> c
            | None      -> 1
        if files.Length = 1 then
            logger.Info (sprintf "single file")
            logger.Trace (sprintf "Processing %s" files.[0])
            psmStats dc p o dbConnection files.[0]
        else
            logger.Info (sprintf "multiple files")
            logger.Trace (sprintf "Processing multiple files: %A" files)
            logger.Trace (sprintf "Program is running on %i cores" c)
            files 
            |> FSharpAux.PSeq.map (psmStats dc p o dbConnection) 
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        let labelled = p.StoreIsotopePatterns && not (SearchDB.getSDBParamsByCn dbConnection).IsotopicMod.IsEmpty
        dbConnection.Dispose()
        if p.StoreIsotopePatterns then
            // The ions of all files, each once; with isotopic modifications both global modifications are stored, so the global modification does not distinguish ions.
            let ions =
                files
                |> Array.collect (fun psms ->
                    let qpsms = Path.Combine [|o; (Path.GetFileNameWithoutExtension psms) + ".qpsm"|]
                    if File.Exists qpsms then readIons qpsms else [||])
                |> Array.distinctBy (fun (sequence, globalMod, charge) -> sequence, (if labelled then 0 else globalMod), charge)
            // Round-robin parts of the union, at most one per core, so that the work does not depend on which file has most of the ions.
            let cores = max 1 c
            let parts = ions |> Array.indexed |> Array.groupBy (fun (i, _) -> i % cores) |> Array.map (snd >> Array.map snd)
            logger.Info (sprintf "Storing the isotope patterns of %i ions in the peptide data base, in %i part(s)" ions.Length parts.Length)
            parts
            |> FSharpAux.PSeq.map (IsotopePatterns.store logger d)
            |> FSharpAux.PSeq.withDegreeOfParallelism cores
            |> Array.ofSeq
            |> ignore
        logger.Info "Done"
        0
