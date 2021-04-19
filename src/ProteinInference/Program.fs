namespace ProteomIQon

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
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name)) 
        let results = parser.Parse argv
        let i'    = results.GetResult InputFolder
        let d'    = results.GetResult PeptideDataBase
        let gff3' = results.TryGetResult GFF3
        let o'    = results.GetResult OutputDirectory
        let p'    = results.GetResult ParamFile
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let d = Path.Combine(directory, d')
        let gff3 =
            match gff3' with
            | Some path ->
                Some (Path.Combine(directory, path))
            | None -> None
        let o = Path.Combine(directory, o')
        let p = Path.Combine(directory, p')
        Logging.generateConfig o
        let logger = Logging.createLogger "ProteinInference"
        logger.Info (sprintf "InputFilePath -i = %A" i')
        logger.Info (sprintf "Peptide data base -d = %s" d)
        //logger.Info (sprintf "InputGFF3Path -g = %s" gff3)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "InputParameterPath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        //let subFolder = (i.Split ([|"\\"|], System.StringSplitOptions.None))
        //Directory.CreateDirectory(o + "\\" + subFolder.[subFolder.Length - 1]) |> ignore
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
            let parsePath (extension: string) fileOrDirectoryPath  = 
                if File.Exists fileOrDirectoryPath then
                    [|fileOrDirectoryPath |]
                elif Directory.Exists fileOrDirectoryPath  then 
                    let fps = Directory.GetFiles(fileOrDirectoryPath ,(extension))
                    fps
                else 
                    [||]
            let parsePaths (extension: string) (fileOrDirectoryPaths:seq<string>) =
                fileOrDirectoryPaths
                |> Seq.collect (parsePath extension)
            parsePaths "qpsm" i'
            |> Array.ofSeq
        ProteinInference.inferProteins gff3 dbConnection proteinInferenceParams o files
        logger.Info "Done"
        0