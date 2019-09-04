namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<AltCommandLine("-i")>] InputFolder      of path:string
        | [<AltCommandLine("-f")>] FastA            of path:string
        | [<AltCommandLine("-g")>] GFF3             of path:string
        | [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<AltCommandLine("-p")>] ParamFile        of path:string

    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InputFolder _      -> "specify folder with input files"
                | FastA _            -> "specify FastA file"
                | GFF3 _             -> "specify GFF3 file"
                | OutputDirectory  _ -> "specify output directory"
                | ParamFile _        -> "specify param file for protein inference"

