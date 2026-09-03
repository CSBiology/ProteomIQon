namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<AltCommandLine("-i")>] InputFile of path:string
        | [<AltCommandLine("-o")>] OutputFile  of path:string 
        | [<AltCommandLine("-pc")>] ProteinColumn of path:string
        | [<AltCommandLine("-rac")>] RatioColumnEnding of string
        | [<AltCommandLine("-rfc")>] ReferenceColumnEnding of string
        | [<AltCommandLine("-l")>] Log_Level of level:int
        | [<AltCommandLine("-v")>] Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InputFile _               -> "specify the input file path"
                | OutputFile  _             -> "specify the output file path"
                | ProteinColumn _           -> "specify the name of the protein name column"
                | RatioColumnEnding _       -> "specify the ending of the ratio name column"
                | ReferenceColumnEnding _   -> "specify the ending of the reference quant column (e.g. heavy or light)"
                | Log_Level _               -> "set the log level."
                | Verbosity_Level _         -> "set the verbosity level."