namespace ProteomIQon

open Argu

module CLIArgumentParsing = 

    type CLIArguments =
        | [<AltCommandLine("-i")>] FastaPath  of path:string 
        | [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<AltCommandLine("-p")>] ParamFile of path:string
        //| [<AltCommandLine("-l")>] Log_Level of level:int
        //| [<AltCommandLine("-v")>] Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | FastaPath _        -> "Specify fasta file path"
                | OutputDirectory _  -> "Specify peptide data base output directory"
                | ParamFile _        -> "Specify param file for the creation of the "
                //| Log_Level _        -> "set the log level."
                //| Verbosity_Level _  -> "set the verbosity level."
