namespace ProteomIQon

open Argu

module CLIArgumentParsing = 

    type CLIArguments =
        | [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<AltCommandLine("-p")>] ParamFile of path:string
        | [<AltCommandLine("-l")>] Log_Level of level:int
        | [<AltCommandLine("-v")>] Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | OutputDirectory  _ -> "specify peptide data base output directory"
                | ParamFile _        -> "specify param file For centroidization"
                | Log_Level _        -> "set the log level."
                | Verbosity_Level _  -> "set the verbosity level."
