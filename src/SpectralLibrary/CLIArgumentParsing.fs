namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>] InstrumentOutput    of path:string
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputDirectory     of path:string
        | [<Mandatory>] [<AltCommandLine("-p")>] ParamFile           of path:string
        | [<Mandatory>] [<AltCommandLine("-d")>] PeptideDataBase     of path:string
        | [<Mandatory>] [<AltCommandLine("-q")>] PSMStatisticsResult of path:string
        | [<Unique>]    [<AltCommandLine("-c")>] Parallelism_Level of level:int

    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InstrumentOutput _    -> "Specify mass spectrometry Output"
                | OutputDirectory  _    -> "Specify output directory"
                | ParamFile _           -> "Specify param file For centroidization"
                | PeptideDataBase _     -> "Specify the file path of the peptide data base."
                | PSMStatisticsResult _ -> "Specify the qpsm file/folder"
                | Parallelism_Level _   -> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."
