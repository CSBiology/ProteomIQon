namespace ProteomIQon

open System.IO
open Argu

module CLIArgumentParsing = 

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>] QuantFiles of path:string list
        | [<Mandatory>] [<AltCommandLine("-ii")>] ProtFiles of path:string list
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputDirectory  of path:string 
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | QuantFiles _    -> "Specify the alignment based quantification output, either specify a file directory containing the files to be analyzed or specify the file path of a single quant file."
                | ProtFiles _           -> "Specify the protein inference output, either specify a file directory containing the files to be analyzed or specify the file path of a single quant file."
                | OutputDirectory  _    -> "Specify the output directory."









