namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<AltCommandLine("-i")>] SpectralLibrary of path:string
        | [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<AltCommandLine("-p")>] ParamFile of path:string

    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | SpectralLibrary _  -> "specify spectral libraries"
                | OutputDirectory  _ -> "specify output directory"
                | ParamFile _        -> "specify param file for consensus library generation"
