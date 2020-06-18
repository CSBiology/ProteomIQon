namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<AltCommandLine("-iq")>] QuantFile of path:string
        | [<AltCommandLine("-ip")>] ProtFile of path:string
        | [<AltCommandLine("-o")>]  OutputDirectory of path:string 
        | [<AltCommandLine("-p")>]  ParamFile of path:string
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | QuantFile _        -> "specify Quant file directory"
                | ProtFile _         -> "specify Prot file directory"
                | OutputDirectory  _ -> "specify output directory"
                | ParamFile _        -> "specify param file for TableSort"

