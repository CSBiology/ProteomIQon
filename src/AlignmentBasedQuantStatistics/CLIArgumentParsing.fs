namespace ProteomIQon

open System.IO
open Argu

module CLIArgumentParsing = 

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>]   Quant of path:string list
        | [<Mandatory>] [<AltCommandLine("-ii")>]  Align of path:string list
        | [<Mandatory>] [<AltCommandLine("-iii")>] AlignedQuant of path:string list
        | [<Mandatory>] [<AltCommandLine("-l")>]   QuantLearn of path:string list
        | [<Mandatory>] [<AltCommandLine("-ll")>]  AlignLearn of path:string list
        | [<Mandatory>] [<AltCommandLine("-lll")>] AlignedQuantLearn of path:string list
        | [<Mandatory>] [<AltCommandLine("-o")>]   OutputDirectory  of path:string 
        | [<Mandatory>] [<AltCommandLine("-p")>]   ParamFile of path:string
        | [<Unique>]    [<AltCommandLine("-c")>]   Parallelism_Level of level:int
        | [<Unique>]    [<AltCommandLine("-dc")>]  DiagnosticCharts 
        | [<Unique>]    [<AltCommandLine("-mf")>]  MatchFiles
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | Quant _             -> "Specify the quant files or directory"
                | Align _             -> "Specify the align files or directory"
                | AlignedQuant _      -> "Specify the aligned quant files or directory"
                | QuantLearn _        -> "Specify the quant files or directory for learning"
                | AlignLearn _        -> "Specify the align files or directory for learning"
                | AlignedQuantLearn _ -> "Specify the aligned quant files or directory for learning"
                | OutputDirectory  _  -> "Specify the output directory."
                | ParamFile _         -> "Specify parameter file for peptide spectrum matching."
                | Parallelism_Level _ -> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."
                | DiagnosticCharts    -> "Set to save diagnostic charts to the output directory."
                | MatchFiles          -> "If this flag is set the files specified by Quant, Align AlignedQuant are matched according to their file name, otherwise they are matched by their position in the input lists."
                







