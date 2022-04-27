namespace ProteomIQon

open Argu

module CLIArgumentParsing = 
    open System.IO
  
    type CLIArguments =
        | [<AltCommandLine("-i")>] InputFolder      of path:string list
        | [<AltCommandLine("-d")>] PeptideDataBase  of path:string 
        | [<AltCommandLine("-g")>] GFF3             of path:string 
        | [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<AltCommandLine("-p")>] ParamFile        of path:string
        | [<Unique>] [<AltCommandLine("-dc")>] DiagnosticCharts 

    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InputFolder _      -> "specify folder with input files"
                | PeptideDataBase _  -> "Specify the file path of the peptide data base."
                | GFF3 _             -> "specify GFF3 file"
                | OutputDirectory  _ -> "specify output directory"
                | ParamFile _        -> "specify param file for protein inference"
                | DiagnosticCharts _ -> "Set to save diagnostic charts to the output directory."

