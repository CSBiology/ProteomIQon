namespace ProteomIQon

open System.IO
open Argu

module CLIArgumentParsing = 

    type CLIArguments =
        | [<Mandatory>] [<AltCommandLine("-i")>] InstrumentOutput of path:string list
        | [<Mandatory>] [<AltCommandLine("-d")>] PeptideDataBase of path:string 
        | [<Mandatory>] [<AltCommandLine("-o")>] OutputDirectory  of path:string 
        | [<Mandatory>] [<AltCommandLine("-p")>] ParamFile of path:string
        | [<Unique>]    [<AltCommandLine("-c")>] Parallelism_Level of level:int
        | [<Unique>]    [<AltCommandLine("-f")>] Max_Fragment_Mass of mass:float
        | [<Unique>]    [<AltCommandLine("-s")>] Slices of count:int
        | [<Unique>]    [<AltCommandLine("-k")>] Keep_Spectra
        | [<Unique>] Mobility_Merge_Ppm of ppm:float
        | [<Unique>]    [<AltCommandLine("-l")>] Log_Level of level:int
        | [<Unique>]    [<AltCommandLine("-v")>] Verbosity_Level of level:int
    with
        interface IArgParserTemplate with
            member s.Usage =
                match s with
                | InstrumentOutput _    -> "Specify the mass spectrometry output, either a directory that contains mzlite files or the paths of single mzlite files."
                | PeptideDataBase  _    -> "Specify the file path of the peptide data base."
                | OutputDirectory  _    -> "Specify the output directory."
                | ParamFile _           -> "Specify the parameter file for peptide spectrum matching."
                | Keep_Spectra          -> "Keep the preprocessed spectra in memory from the first slice on, so the later slices read nothing from the file. About 2.4 kB per spectrum and charge."
                | Mobility_Merge_Ppm _  -> "m/z cluster width for mandatory within-spectrum mobility summation. Defaults to FragmentTolerancePPM from the search parameters (normally 20 ppm). Override for the instrument's accuracy. Uses summed intensity and intensity-weighted m/z. 0 sums exact coordinates and does not disable summation."
                | Slices _              -> "Number of slices the fragment index is built and searched in, one after the other. More slices need less memory and one more pass over the spectra each. Default 1."
                | Max_Fragment_Mass _   -> "Highest neutral fragment mass that enters the index, in Dalton. Without it the ceiling is the highest precursor window including isotope errors and tolerance, plus 1 Da. A lower value saves memory and is warned about when it is below that bound."
                | Parallelism_Level _   -> "Number of files searched together, one search thread per file. Default 1."
                | Log_Level _           -> "Set the log level."
                | Verbosity_Level _     -> "Set the verbosity level."
