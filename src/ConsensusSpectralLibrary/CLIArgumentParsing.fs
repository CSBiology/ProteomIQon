namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<AltCommandLine("-i")>] InstrumentOutput of path:string
        | [<AltCommandLine("-ii")>] SpectralLibrary of path:string
        | [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<AltCommandLine("-p")>] ParamFile of path:string
        | [<Unique>]    [<AltCommandLine("-dc")>] DiagnosticCharts 

    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InstrumentOutput _  -> "specify spectral libraries"
                | SpectralLibrary _  -> "specify spectral libraries"
                | OutputDirectory  _ -> "specify output directory"
                | ParamFile _        -> "specify param file for consensus library generation"
                | DiagnosticCharts _    -> "Set to save diagnostic charts to the output directory."