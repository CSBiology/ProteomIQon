namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>]  TargetFiles of path:string list
        | [<Mandatory>] [<AltCommandLine("-ii")>] SourceFiles of path:string list
        | [<Mandatory>] [<AltCommandLine("-o")>]  OutputDirectory  of path:string 
        | [<Unique>]    [<AltCommandLine("-dc")>] DiagnosticCharts 
        | [<Unique>]    [<AltCommandLine("-c")>]  Parallelism_Level of level:int
        | [<Unique>]    [<AltCommandLine("-l")>]  Log_Level of level:int
        | [<Unique>]    [<AltCommandLine("-v")>]  Verbosity_Level of level:int
        | [<Unique>]    [<AltCommandLine("-s")>]  Seed of seed:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | TargetFiles _         -> "Specify quantified peptides, specify a directory containing .quant files."
                | SourceFiles _         -> "Specify quantified peptides, specify a directory containing .quant files."
                | OutputDirectory  _    -> "Specify the output directory."
                | DiagnosticCharts      -> "Set to save diagnostic charts to the output directory."
                | Parallelism_Level _   -> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."
                | Log_Level _           -> "Set the log level."
                | Verbosity_Level _     -> "Set the verbosity level."
                | Seed _                -> "Seed of the cross-validation folds (default 1)."
