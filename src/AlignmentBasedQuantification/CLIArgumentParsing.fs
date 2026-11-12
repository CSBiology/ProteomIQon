namespace ProteomIQon

open Argu

module CLIArgumentParsing =
    open System.IO

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>]    InstrumentOutput of path:string
        | [<Mandatory>] [<AltCommandLine("-ii")>]   AlignedPeptides of path:string
        | [<Mandatory>] [<AltCommandLine("-iv")>]   QuantifiedPeptides of path:string
        | [<Mandatory>] [<AltCommandLine("-d")>]    PeptideDataBase of path:string
        | [<Mandatory>] [<AltCommandLine("-o")>]    OutputDirectory  of path:string
        | [<Mandatory>] [<AltCommandLine("-p")>]    ParamFile of path:string
        | [<Unique>]    [<AltCommandLine("-mf")>]   MatchFiles
        | [<Unique>]    [<AltCommandLine("-dc")>]   DiagnosticCharts
        | [<Unique>]    [<AltCommandLine("-c")>]    Parallelism_Level of level:int
        | [<Unique>]    [<AltCommandLine("-l")>]    Log_Level of level:int
        | [<Unique>]    [<AltCommandLine("-v")>]    Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InstrumentOutput _    -> "Specify the mass spectrometry output, either specify a file directory containing the .mzlite files to be analyzed or specify the file path of a single mzlite file."
                | AlignedPeptides  _    -> "Specify the .align files of QuantBasedAlignment, either a directory or the file path of a single .align file."
                | QuantifiedPeptides  _ -> "Specify the quantified peptides (.quant) of the PSM based quantification, either a directory or the file path of a single .quant file."
                | PeptideDataBase  _    -> "Specify the file path of the peptide data base."
                | OutputDirectory  _    -> "Specify the output directory."
                | ParamFile _           -> "Specify the AlignmentBasedQuantification parameter file."
                | MatchFiles            -> "If this flag is set the files specified by InstrumentOutput, AlignedPeptides and QuantifiedPeptides are matched according to their file name, otherwise they are matched by their position in the input lists."
                | DiagnosticCharts      -> "Set to save diagnostic charts to the output directory."
                | Parallelism_Level _   -> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."
                | Log_Level _           -> "Set the log level."
                | Verbosity_Level _     -> "Set the verbosity level."
