namespace ProteomIQon

open System
open System.IO
open System.Reflection
open Argu
open ProteomIQon.Core
open ProteomIQon.Core.InputPaths
open AlignmentBasedQuantStatistics
open CLIArgumentParsing

module console1 =

    [<EntryPoint>]
    let main argv = 
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i = results.GetResult Quant |> List.map getPathRelativeToDir
        let ii = results.GetResult Align |> List.map getPathRelativeToDir
        let iii = results.GetResult AlignedQuant |> List.map getPathRelativeToDir
        let o = results.GetResult OutputDirectory  |> getPathRelativeToDir
        let p = results.GetResult ParamFile        |> getPathRelativeToDir
        let dc = results.Contains DiagnosticCharts
        let mf = results.Contains MatchFiles
        let c =
            match results.TryGetResult Parallelism_Level with
            | Some c    -> c
            | None      -> 1
        Directory.CreateDirectory(o) |> ignore
        Logging.generateConfig o
        let logger = Logging.createLogger "AlignmentBasedQuantStatistics"
        logger.Info (sprintf "QuantFilePath -i = %A" i)
        logger.Info (sprintf "AlignFilePath -i = %A" ii)
        logger.Info (sprintf "AlignedQuantFilePath -i = %A" iii)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let p = 
            Json.ReadAndDeserialize<Dto.PeptideSpectrumMatchingParams> p
            |> Dto.PeptideSpectrumMatchingParams.toDomain
        if i.Length = 1 && File.Exists i.[0] then
            learnScore ([|i.[0], ii.[0], iii.[0]|]) logger 1 dc o
        elif (i.Length = 1 && Directory.Exists i.[0]) || i.Length > 1 then
            let quantFiles =
                if i.Length = 1 then
                    Directory.GetFiles(i.[0],("*.quant"))
                else
                    i
                    |> Array.ofList
            let alignFiles =
                if ii.Length = 1 then
                    Directory.GetFiles(ii.[0],("*.align"))
                else
                    ii
                    |> Array.ofList
            let alignQuantFiles =
                if iii.Length = 1 then
                    Directory.GetFiles(iii.[0],("*.quant"))
                else
                    iii
                    |> Array.ofList
            let matchedFiles =
                    if mf then 
                        quantFiles
                        |> Array.choose (fun quantFilePath ->
                            match alignFiles |> Array.tryFind (fun alignfilePath -> (Path.GetFileNameWithoutExtension alignfilePath) = (Path.GetFileNameWithoutExtension quantFilePath)) with
                            | Some alignfilePath -> 
                                match alignQuantFiles |> Array.tryFind (fun alignQuantFilePath -> (Path.GetFileNameWithoutExtension alignQuantFilePath) = (Path.GetFileNameWithoutExtension alignfilePath)) with
                                | Some(alignQuantFilePath) ->
                                    Some(quantFilePath,alignfilePath,alignQuantFilePath)
                                | None -> 
                                    logger.Trace (sprintf "no alignQuant file for %s" alignfilePath)
                                    None
                            | None ->
                                logger.Trace (sprintf "no alignment file for %s" quantFilePath)
                                None
                            )
                    else 
                        [|for i = 0 to i.Length-1 do yield quantFiles.[i], alignFiles.[i], alignQuantFiles.[i]|]
            learnScore matchedFiles logger c dc o
        else
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
          
        logger.Info "Done"
        0
        