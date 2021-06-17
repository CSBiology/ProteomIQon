namespace ProteomIQon

open System.IO
open Argu

module CLIArgumentParsing = 

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>] InstrumentOutput of path:string
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<Mandatory>] [<AltCommandLine("-p")>] ParamFile of path:string
        | [<Unique>]  [<AltCommandLine("-c")>] Parallelism_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InstrumentOutput _    -> "Specify the mass spectrometry output, either specify a file directory containing the files to be analyzed or specify the file path of a single mzlite file."
                | OutputDirectory  _    -> "Specify the output directory."
                | ParamFile _           -> "Specify parameter file for peptide spectrum matching."
                | Parallelism_Level _   -> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."








