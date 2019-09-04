namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
    
    type CLIArguments =
        | [<Mandatory>] [<Unique>] [<AltCommandLine("-i")>] InstrumentOutput of filePath:string
        | [<Mandatory>] [<Unique>] [<AltCommandLine("-o")>] OutputDirectory  of directoryPath:string 
        | [<Mandatory>] [<Unique>] [<AltCommandLine("-p")>] ParamFile of filePath:string
        | [<Unique>]  [<AltCommandLine("-c")>] Parallelism_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InstrumentOutput _ -> "Specify the mass spectrometry output, either specify a file directory containing the files to be analyzed or specify the file path of a single mzlite file."
                | OutputDirectory  _ -> "Specify output directory"
                | ParamFile _        -> "Specify param file For centroidization"
                | Parallelism_Level _   -> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."

