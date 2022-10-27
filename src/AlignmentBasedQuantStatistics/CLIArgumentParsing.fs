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
        | [<Unique>]    [<AltCommandLine("-c")>]    Parallelism_Level of level:int
        | [<Unique>]    [<AltCommandLine("-dc")>]   DiagnosticCharts 
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | Quant _             -> "Specify the quant files or directory"
                | Align _             -> "Specify the align files or directory"
                | AlignedQuant _      -> "Specify the aligned quant files or directory"
                | OutputDirectory  _  -> "Specify the output directory."
                | ParamFile _         -> "Specify parameter file for peptide spectrum matching."
                | Parallelism_Level _ -> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."
                | DiagnosticCharts _    -> "Set to save diagnostic charts to the output directory."








