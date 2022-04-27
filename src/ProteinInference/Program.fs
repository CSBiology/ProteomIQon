namespace ProteomIQon

open System
open System.IO
open CLIArgumentParsing
open Argu
open BioFSharp.Mz
open System.Reflection
open FSharpAux
open ProteomIQon.Core.InputPaths

module console1 =

    [<EntryPoint>]
    let main argv = 
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i    = results.GetResult InputFolder |> List.map getPathRelativeToDir
        let d    = results.GetResult PeptideDataBase |> getPathRelativeToDir
        let o    = results.GetResult OutputDirectory |> getPathRelativeToDir
        let p    = results.GetResult ParamFile |> getPathRelativeToDir
        let dc   = results.Contains DiagnosticCharts
        let gff3 =
            let path = results.TryGetResult GFF3
            match path with
            | Some path -> Some (getPathRelativeToDir path)
            | None      -> None
        Logging.generateConfig o
        let logger = Logging.createLogger "ProteinInference"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "Peptide data base -d = %s" d)
        //logger.Info (sprintf "InputGFF3Path -g = %s" gff3)
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
        let proteinInferenceParams = 
                Json.ReadAndDeserialize<Dto.ProteinInferenceParams> p
                |> Dto.ProteinInferenceParams.toDomain
        let files = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.qpsm"))) i
            |> Array.ofSeq
        logger.Trace (sprintf "Inputfiles: %A" files)
        ProteinInference.inferProteins dc gff3 dbConnection proteinInferenceParams o files
        logger.Info "Done"
        0