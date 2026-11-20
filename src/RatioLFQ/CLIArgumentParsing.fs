namespace ProteomIQon

open Argu

module CLIArgumentParsing =

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>] IonTable of path:string
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputDirectory of path:string
        | [<Mandatory>] [<AltCommandLine("-p")>] ParamFile of path:string
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | IonTable _        -> "Specify the peptide ion table, GlobModAggregation.txt of LabeledProteinQuantification written without a transform."
                | OutputDirectory _ -> "Specify the output directory. RatioLFQ.txt and the log file are written into it."
                | ParamFile _       -> "Specify parameter file."
