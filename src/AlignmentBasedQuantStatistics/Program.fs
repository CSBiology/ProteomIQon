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
        let l = results.GetResult QuantLearn |> List.map getPathRelativeToDir
        let ll = results.GetResult AlignLearn |> List.map getPathRelativeToDir
        let lll = results.GetResult AlignedQuantLearn |> List.map getPathRelativeToDir
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
        logger.Info (sprintf "AlignFilePath -ii = %A" ii)
        logger.Info (sprintf "AlignedQuantFilePath -iii = %A" iii)
        logger.Info (sprintf "QuantFilePathForLearning -l = %A" l)
        logger.Info (sprintf "AlignFilePathForLearning -ll = %A" ll)
        logger.Info (sprintf "AlignedQuantFilePathForLearning -lll = %A" lll)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        let p = 
            Json.ReadAndDeserialize<Dto.AlignmentBasedQuantStatisticsParams> p
            |> Dto.AlignmentBasedQuantStatisticsParams.toDomain
        if i.Length = 1 && File.Exists i.[0] && l.Length = 1 && File.Exists l.[0] then
            logger.Trace "Single file"
            assignScoreAndQValue (i.[0], ii.[0], iii.[0]) [|l.[0], ll.[0], lll.[0]|] dc p o
            |> ignore
        elif ((i.Length = 1 && Directory.Exists i.[0]) || i.Length > 1) && ((l.Length = 1 && Directory.Exists i.[0]) || l.Length > 1) then
            logger.Trace "Multiple files"
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
            let quantFilesLearning =
                if l.Length = 1 then
                    Directory.GetFiles(l.[0],("*.quant"))
                else
                    l
                    |> Array.ofList
            let alignFilesLearning =
                if ll.Length = 1 then
                    Directory.GetFiles(ll.[0],("*.align"))
                else
                    ll
                    |> Array.ofList
            let alignQuantFilesLearning =
                if lll.Length = 1 then
                    Directory.GetFiles(lll.[0],("*.quant"))
                else
                    lll
                    |> Array.ofList
            logger.Trace $"Quant Files: {quantFiles}"
            logger.Trace $"Alignment Files: {alignFiles}"
            logger.Trace $"Aligned Quant Files: {alignQuantFiles}"
            logger.Trace $"Quant Files Learning: {quantFilesLearning}"
            logger.Trace $"Alignment Files Learning: {alignFilesLearning}"
            logger.Trace $"Aligned Quant Files Learning: {alignQuantFilesLearning}"
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
                    [|for i = 0 to quantFiles.Length-1 do yield quantFiles.[i], alignFiles.[i], alignQuantFiles.[i]|]
            let matchedFilesLearning =
                if mf then 
                    quantFilesLearning
                    |> Array.choose (fun quantFilePath ->
                        match alignFilesLearning |> Array.tryFind (fun alignfilePath -> (Path.GetFileNameWithoutExtension alignfilePath) = (Path.GetFileNameWithoutExtension quantFilePath)) with
                        | Some alignfilePath -> 
                            match alignQuantFilesLearning |> Array.tryFind (fun alignQuantFilePath -> (Path.GetFileNameWithoutExtension alignQuantFilePath) = (Path.GetFileNameWithoutExtension alignfilePath)) with
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
                    [|for i = 0 to quantFiles.Length-1 do yield quantFilesLearning.[i], alignFilesLearning.[i], alignQuantFilesLearning.[i]|]
            matchedFiles
            |> FSharpAux.PSeq.map (fun (matchedFile) -> assignScoreAndQValue matchedFile matchedFilesLearning dc p o)
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        elif ((i.Length = 1 && Directory.Exists i.[0]) || i.Length > 1) && (l.Length = 1 && File.Exists l.[0]) then
            logger.Trace "Multiple files with single file for learning"
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
            logger.Trace $"Quant Files: {quantFiles}"
            logger.Trace $"Alignment Files: {alignFiles}"
            logger.Trace $"Aligned Quant Files: {alignQuantFiles}"
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
                    [|for i = 0 to quantFiles.Length-1 do yield quantFiles.[i], alignFiles.[i], alignQuantFiles.[i]|]
            matchedFiles
            |> FSharpAux.PSeq.map (fun (matchedFile) -> assignScoreAndQValue matchedFile [|l.[0], ll.[0], lll.[0]|] dc p o)
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> Array.ofSeq
            |> ignore
        elif ((i.Length = 1 && File.Exists i.[0]) || i.Length > 1) && ((l.Length = 1 && Directory.Exists i.[0]) || l.Length > 1) then
            logger.Trace "Single file with multiple files for learning"
            let quantFilesLearning =
                if l.Length = 1 then
                    Directory.GetFiles(l.[0],("*.quant"))
                else
                    l
                    |> Array.ofList
            let alignFilesLearning =
                if ll.Length = 1 then
                    Directory.GetFiles(ll.[0],("*.align"))
                else
                    ll
                    |> Array.ofList
            let alignQuantFilesLearning =
                if lll.Length = 1 then
                    Directory.GetFiles(lll.[0],("*.quant"))
                else
                    lll
                    |> Array.ofList
            logger.Trace $"Quant Files Learning: {quantFilesLearning}"
            logger.Trace $"Alignment Files Learning: {alignFilesLearning}"
            logger.Trace $"Aligned Quant Files Learning: {alignQuantFilesLearning}"
            let matchedFilesLearning =
                if mf then 
                    quantFilesLearning
                    |> Array.choose (fun quantFilePath ->
                        match alignFilesLearning |> Array.tryFind (fun alignfilePath -> (Path.GetFileNameWithoutExtension alignfilePath) = (Path.GetFileNameWithoutExtension quantFilePath)) with
                        | Some alignfilePath -> 
                            match alignQuantFilesLearning |> Array.tryFind (fun alignQuantFilePath -> (Path.GetFileNameWithoutExtension alignQuantFilePath) = (Path.GetFileNameWithoutExtension alignfilePath)) with
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
                    [|for i = 0 to quantFilesLearning.Length-1 do yield quantFilesLearning.[i], alignFilesLearning.[i], alignQuantFilesLearning.[i]|]
            assignScoreAndQValue (i.[0], ii.[0], iii.[0]) matchedFilesLearning dc p o
            |> ignore
        else
            failwith "The given path to the instrument output is neither a valid file path nor a valid directory path."
          
        logger.Info "Done"
        0
        