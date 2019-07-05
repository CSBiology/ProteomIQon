namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
    
    type IOType =
        | Directory of string
        | Path      of string

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>] InstrumentOutput of path:string
        | [<Mandatory>] [<AltCommandLine("-d")>] PeptideDataBase of path:string 
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<Mandatory>] [<AltCommandLine("-p")>] ParamFile of path:string
        | [<Unique>]    [<AltCommandLine("-c")>] Parallelism_Level of level:int
        | [<Unique>]    [<AltCommandLine("-l")>] Log_Level of level:int
        | [<Unique>]    [<AltCommandLine("-v")>] Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InstrumentOutput _    -> "Specify the mass spectrometry output, either specify a file directory containing the files to be analyzed or specify the file path of a single mzlite file."
                | PeptideDataBase  _    -> "Specify the file path of the peptide data base."
                | OutputDirectory  _    -> "Specify the output directory."
                | ParamFile _           -> "Specify parameter file for peptide spectrum matching."
                | Parallelism_Level _   -> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."
                | Log_Level _           -> "Set the log level."
                | Verbosity_Level _     -> "Set the verbosity level."
