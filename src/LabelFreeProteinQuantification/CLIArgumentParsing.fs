namespace ProteomIQon

open System.IO
open Argu

module CLIArgumentParsing = 

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>] ProteinAssignedQuantPepIons of path:string list
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<Mandatory>] [<AltCommandLine("-p")>] ParamFile of path:string
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | ProteinAssignedQuantPepIons _ -> "Specify the a single, list or directory containing .quantAndProt files."
                | OutputDirectory  _            -> "Specify the output directory."
                | ParamFile _                   -> "Specify parameter file."









