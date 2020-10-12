namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-itab")>] TableSortFile of path:string
        | [<Mandatory>] [<AltCommandLine("-iprot")>] ProtFile of path:string
        | [<Mandatory>] [<AltCommandLine("-ipep")>] QuantFile of path:string
        | [<Mandatory>] [<AltCommandLine("-ipsm")>] QpsmFile of path:string
        | [<Mandatory>] [<AltCommandLine("-d")>] PeptideDataBase of path:string 
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputFile  of path:string 
        | [<Mandatory>] [<AltCommandLine("-p")>] ParamFile of path:string
        | [<AltCommandLine("-l")>] Log_Level of level:int
        | [<AltCommandLine("-v")>] Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | TableSortFile _    -> "specify TableSort file location"
                | ProtFile _         -> "specify Prot file location"
                | QuantFile _        -> "specify Quant file location"
                | QpsmFile _         -> "specify Qpsm file location"
                | OutputFile  _      -> "specify output file path"
                | ParamFile _        -> "specify param file For centroidization"
                | Log_Level _        -> "set the log level."
                | PeptideDataBase  _ -> "Specify the file path of the peptide data base."
                | Verbosity_Level _  -> "set the verbosity level."
