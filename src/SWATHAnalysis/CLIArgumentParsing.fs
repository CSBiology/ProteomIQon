namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
    

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>]  InstrumentOutput of path:string
        | [<Mandatory>] [<AltCommandLine("-ii")>] ConsensusLibary of path:string
        | [<Mandatory>] [<AltCommandLine("-d")>]  PeptideDataBase of path:string 
        | [<Mandatory>] [<AltCommandLine("-o")>]  OutputDirectory  of path:string 
        | [<Mandatory>] [<AltCommandLine("-p")>]  ParamFile of path:string
        | [<Unique>]    [<AltCommandLine("-mf")>] MatchFiles
        | [<Unique>]    [<AltCommandLine("-dc")>] DiagnosticCharts 
        | [<Unique>]    [<AltCommandLine("-c")>]  Parallelism_Level of level:int
        | [<Unique>]    [<AltCommandLine("-l")>]  Log_Level of level:int
        | [<Unique>]    [<AltCommandLine("-v")>]  Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InstrumentOutput _    -> "Specify the mass spectrometry output, either specify a file directory containing the files to be analyzed or specify the file path of a single mzlite file."
                | ConsensusLibary  _    -> "Specify the consensuslibrary, either specify a file directory containing the files to be analyzed or specify the file path of a single qpsm file. If InstrumentOutPut(i) and ScoredPSMs (ii) both reference a directory the files are automatically aligned by their file name."
                | PeptideDataBase  _    -> "Specify the file path of the peptide data base."
                | OutputDirectory  _    -> "Specify the output directory."
                | ParamFile _           -> "Specify parameter file for peptide spectrum matching."
                | MatchFiles            -> "If this flag is set the files specified using the QuantifiedPeptides and InferredProteins are matched according to their file name, otherwise they are matched by their position in the input lists."
                | DiagnosticCharts _    -> "Set to save diagnostic charts to the output directory."
                | Parallelism_Level _   -> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."
                | Log_Level _           -> "Set the log level."
                | Verbosity_Level _     -> "Set the verbosity level."
