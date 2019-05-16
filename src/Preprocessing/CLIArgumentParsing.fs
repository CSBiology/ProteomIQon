namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
    
    type CLIArguments =
        | [<Unique>] [<AltCommandLine("-i")>] InstrumentOutput of filePath:string
        | [<Unique>] [<AltCommandLine("-o")>] OutputDirectory  of directoryPath:string 
        | [<Unique>] [<AltCommandLine("-p")>] ParamFile of filePath:string
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InstrumentOutput _ -> "specify mass spectrometry Output"
                | OutputDirectory  _ -> "specify output directory"
                | ParamFile _        -> "specify param file For centroidization"

