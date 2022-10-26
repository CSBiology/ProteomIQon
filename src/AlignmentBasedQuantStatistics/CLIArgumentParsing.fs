namespace ProteomIQon

open System.IO
open Argu

module CLIArgumentParsing = 

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>]   Quant of path:string list
        | [<Mandatory>] [<AltCommandLine("-ii")>]  Align of path:string list
        | [<Mandatory>] [<AltCommandLine("-iii")>] AlignedQuant of path:string list
        | [<Mandatory>] [<AltCommandLine("-o")>]   OutputDirectory  of path:string 
        | [<Mandatory>] [<AltCommandLine("-p")>]   ParamFile of path:string
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | Quant _            -> "Specify the quant files or directory"
                | Align _            -> "Specify the align files or directory"
                | AlignedQuant _     -> "Specify the aligned quant files or directory"
                | OutputDirectory  _ -> "Specify the output directory."
                | ParamFile _        -> "Specify parameter file for peptide spectrum matching."









