namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
open QuantBasedAlignment
open System.Reflection
open ProteomIQon.Core.InputPaths

module console1 =
    open BioFSharp.Mz

    [<EntryPoint>]
    let main argv =
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i   = results.GetResult TargetFiles |> List.map getPathRelativeToDir
        let ii  = results.GetResult SourceFiles |> List.map getPathRelativeToDir 
        let o   = results.GetResult OutputDirectory    |> getPathRelativeToDir
        let p   = results.GetResult ParamFile          |> getPathRelativeToDir
        let dc  = results.Contains DiagnosticCharts
        Logging.generateConfig o
        let logger = Logging.createLogger "QuantBasedAlignment"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        //let p =
        //    Json.ReadAndDeserialize<Dto.QuantificationParams> p
        //    |> Dto.QuantificationParams.toDomain
        let targetFiles = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.quant"))) i
            |> Array.ofSeq
        let sourceFiles = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.quant"))) i
            |> Array.ofSeq
        if targetFiles.Length = 1 then
            logger.Info "single file detected"
            alignFiles dc {Placeholder=true} o sourceFiles targetFiles.[0]
        else
            logger.Info "directory found"
            logger.Trace (sprintf ".quant files : %A" targetFiles)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            targetFiles
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> FSharpAux.PSeq.iter (fun (t) -> 
                let sourceFiles' = 
                    sourceFiles
                    |> Array.filter (fun x -> x <> t)
                alignFiles dc {Placeholder=true} o sourceFiles' t)
            |> ignore
        logger.Info "Done"
        0
