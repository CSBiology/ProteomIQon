namespace ProteomIQon

open System.IO
open Argu

module CLIArgumentParsing = 

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>] QuantifiedPeptides of path:string list
        | [<Mandatory>] [<AltCommandLine("-ii")>] InferredProteins of path:string list
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<Unique>]    [<AltCommandLine("-mf")>] MatchFiles 
        | [<Unique>]    [<AltCommandLine("-c")>] Parallelism_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | QuantifiedPeptides _  -> "Specify quantified peptides, specify either a directory containing or a space separated list of .quant files."
                | InferredProteins  _   -> "Specify inferred proteins, specify either a directory containing or a space separated list of .prot files."
                | OutputDirectory  _    -> "Specify the output directory."
                | MatchFiles            -> "If this flag is set the files specified using the QuantifiedPeptides and InferredProteins are matched according to their file name, otherwise they are matched by their position in the input lists."
                | Parallelism_Level _   -> "Set the number of cores the programm can use. Parallelization occurs on file level. This flag is only of effect if a input directory (-i) is specified."





