namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<AltCommandLine("-i")>] InstrumentOutput of path:string list
        | [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<AltCommandLine("-p")>] ParamFile of path:string
        | [<Unique>] [<AltCommandLine("-c")>] Parallelism_Level of level:int
        | [<Unique>] [<AltCommandLine("-f")>] FixFiles 
        | [<AltCommandLine("-l")>] Log_Level of level:int
        | [<AltCommandLine("-v")>] Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InstrumentOutput _ -> "specify mass spectrometry Output"
                | OutputDirectory  _ -> "specify output directory"
                | ParamFile _        -> "specify param file For centroidization"
                | Log_Level _        -> "set the log level."
                | Verbosity_Level _  -> "set the verbosity level."
                | Parallelism_Level _-> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."
                | FixFiles           -> "If this flag is set the files specified by InstrumentOutput have the &quot removed."