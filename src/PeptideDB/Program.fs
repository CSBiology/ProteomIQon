namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open System.Reflection

module console1 =

    [<EntryPoint>]
    let main argv =
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name))
        let results = parser.Parse argv
        let inputFasta = results.TryGetResult FastaPath
        let outputDir = results.TryGetResult OutputDirectory
        let paramF = results.TryGetResult ParamFile
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        match inputFasta, outputDir, paramF with
        | Some i', Some o' , Some p' ->
            let i = Path.Combine(directory, i')
            let o = Path.Combine(directory, o')
            let p = Path.Combine(directory, p')
            Logging.generateConfig outputDir.Value
            let logger = Logging.createLogger "PeptideDB"
            logger.Trace (sprintf "CLIArguments %A" results)
            logger.Info (sprintf "OutputFilePath %s" o)
            logger.Info (sprintf "ParamFilePath %s" p)
            Directory.CreateDirectory(o) |> ignore
            let processParams =
                Json.ReadAndDeserialize<Dto.PeptideDBParams> p
                |> Dto.PeptideDBParams.toDomain
            
            PeptideDB.createPeptideDB processParams o i
        | _ -> failwith "parameterfile or output directory have no valid path or are no valid file"
        0