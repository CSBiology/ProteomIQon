namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
open TableSort
open System.Reflection

module console1 =

    [<EntryPoint>]
    let main argv = 
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let results = parser.Parse argv
        let quant' = results.GetResult QuantFile
        let prot'  = results.GetResult ProtFile
        let o'     = results.GetResult OutputDirectory
        let p'     = results.GetResult ParamFile
        let directory = System.IO.Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)
        let quant = Path.Combine(directory, quant')
        let prot = Path.Combine(directory, prot')
        let o = Path.Combine(directory, o')
        let p = Path.Combine(directory, p')
        Logging.generateConfig o
        let logger = Logging.createLogger "TableSort"
        logger.Info (sprintf "QuantFilePath -i = %s" quant)
        logger.Info (sprintf "ProtFilePath = %s" prot)
        logger.Info (sprintf "OutputFolderPath -o = %s" o)
        logger.Info (sprintf "InputParameterPath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let p =
            Json.ReadAndDeserialize<Dto.TableSortParams> p
            |> Dto.TableSortParams.toDomain
        if Directory.Exists quant && Directory.Exists prot then
            logger.Info (sprintf "multiple files")
            let quantFiles =
                Directory.GetFiles(quant,("*.quant"))
            let protFiles =
                Directory.GetFiles(prot,("*.prot"))
            let matchedFiles =
                quantFiles
                |> Array.collect (fun quantFile ->
                    protFiles
                    |> Array.choose (fun protFile ->
                        if Path.GetFileNameWithoutExtension quantFile = Path.GetFileNameWithoutExtension protFile then
                            Some (quantFile, protFile)
                        else
                            None
                    )
                )
            matchedFiles
            |> Array.unzip
            |> fun (quant, prot) -> sortTables quant prot o p
            |> ignore
        elif (File.Exists quant) || (File.Exists prot) then
            failwith "Both, quant file path and prot file path, must be a folder."
        else
            failwith "The given paths are no valid folder paths"
        logger.Trace "Done"
        0
