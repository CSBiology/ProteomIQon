namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>] TableSortFile of path:string
        | [<Mandatory>] [<AltCommandLine("-ii")>] ProtFile of path:string
        | [<Mandatory>] [<AltCommandLine("-iii")>] QuantFile of path:string
        | [<Mandatory>] [<AltCommandLine("-iv")>] QpsmFile of path:string
        //| [<Mandatory>] [<AltCommandLine("-d")>] PeptideDataBase of path:string 
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputFile  of path:string 
        | [<Mandatory>] [<AltCommandLine("-p")>] ParamFile of path:string
        | [<AltCommandLine("-l")>] Log_Level of level:int
        | [<AltCommandLine("-v")>] Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | TableSortFile _    -> "specify TableSort file location"
                | ProtFile _         -> "specify Prot files location"
                | QuantFile _        -> "specify Quant files location"
                | QpsmFile _         -> "specify Qpsm files location"
                | OutputFile  _      -> "specify output file path"
                | ParamFile _        -> "specify param file For centroidization"
                | Log_Level _        -> "set the log level."
                | Verbosity_Level _  -> "set the verbosity level."
