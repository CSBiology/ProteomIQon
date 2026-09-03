namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>] MSFraggerOutput of path:string list
        | [<Mandatory>] [<AltCommandLine("-ii")>] MzLiteOutput of path:string list
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<Mandatory>] [<AltCommandLine("-d")>] PeptideDataBase of path:string 
        | [<Unique>] [<AltCommandLine("-c")>] Parallelism_Level of level:int
        | [<AltCommandLine("-l")>] Log_Level of level:int
        | [<AltCommandLine("-v")>] Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | MSFraggerOutput _  -> "specify validated psm Output"
                | MzLiteOutput _     -> "specify mass spectrometry Output"
                | OutputDirectory  _ -> "specify output directory"
                | PeptideDataBase  _ -> "specify the path to the database"
                | Log_Level _        -> "set the log level."
                | Verbosity_Level _  -> "set the verbosity level."
                | Parallelism_Level _-> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."