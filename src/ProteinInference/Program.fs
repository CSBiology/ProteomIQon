namespace ProteomIQon

open System.IO
open CLIArgumentParsing
open Argu
module console1 =

    [<EntryPoint>]
    let main argv = 
        printfn "%A" argv

        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name)) 
        let results = parser.Parse argv
        let i     = results.GetResult InputFolder
        let fastA = results.GetResult FastA
        let gff3  = results.GetResult GFF3
        let o     = results.GetResult OutputDirectory
        let p     = results.GetResult ParamFile
        Logging.generateConfig o
        let logger = Logging.createLogger "ProteinInference"
        logger.Info (sprintf "InputFilePath -i = %s" i)
        logger.Info (sprintf "InputFastAPath -f = %s" fastA)
        logger.Info (sprintf "InputGFF3Path -g = %s" gff3)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        logger.Info (sprintf "InputParameterPath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        let subFolder = (i.Split ([|"\\"|], System.StringSplitOptions.None))
        Directory.CreateDirectory(o + "\\" + subFolder.[subFolder.Length - 1]) |> ignore
        let proteinInferenceParams = 
                Json.ReadAndDeserialize<Dto.ProteinInferenceParams> p
                |> Dto.ProteinInferenceParams.toDomain
        ProteinInference.inferProteins gff3 fastA proteinInferenceParams o i
        logger.Info "Hit any key to exit."
        System.Console.ReadKey() |> ignore
        0