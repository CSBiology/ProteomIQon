(**
---
title: MzMLToMzLite
category: Tools
categoryindex: 1
index: 1
---
*)
(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# mzMl to mzlite Conversion
**Disclaimer** this tool converts [mzML](https://www.psidev.info/mzML) to [mzLite](https://github.com/CSBiology/MzIO/blob/developer/src/MzIO.SQL/MzIOSQL.fs) 
if you want to convert the other way around please visit the [MzLiteToMzML]({{root}}tools/MzLiteToMzML.html) documentation. 
We recommend the use of [msconvert](https://www.nature.com/articles/nbt.2377) to convert your raw data into .mzML. 
A user friendly way to execute msconvert is available through [Galaxy Europe](https://galaxyproject.eu/posts/2019/03/24/msconvert/).

The success of modern proteomics was made possible by constant progression in the field of mass spectrometry. Over the course of the past years quite a few
manufacturers of mass spectrometers have managed to establish themselfes in the field of biological research. Since aquisition and accession of mass spectra are performance critical processes, 
various performance optimized, but vendor specific and closed source formats have been developed to store raw MS data. 
This comes to the disadvantage for toolchain developers which want to provide tools for every scientist regardless of the format of their raw data.    

In a effort to provide an open format for the storage of MS data the format [mzML](https://www.psidev.info/mzML) was developed. While this XML based format is straight forward to implement
it falls behind in performance critical scenarios. To be competitive in performance and to comply to the [FAIR principles](https://www.go-fair.org/fair-principles/) we chose
to use [mzLite](https://github.com/CSBiology/MzIO/blob/developer/src/MzIO.SQL/MzIOSQL.fs), an open and SQLite based implementation of mzML, within our toolchain.

The tool mzMLToMzLite allows to convert mzML files to mzLite files. Additionally, it allows the user to perform peak picking or filtering of mass spectra.

## Parameters
The following table gives an overview of the parameter set:

| **Parameter**                  | **Default Value**                                                                                                                         | **Description**                                                    |
|--------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------|
| Compress                       | MzIO.Binary.BinaryDataCompressionType.NoCompression                                                                                       | Indicates if peak data should be compressed                        |
| StartRetentionTime             | None                                                                                                                                      | Minimum scan time for spectra to be copied                         |
| EndRetentionTime               | None                                                                                                                                      | Maximum scan time for spectra to be copied                         |
| MS1PeakPicking                 | PeakPicking.ProfilePeaks                                                                                                                  | Parameter to configure peak picking                                |
| MS2PeakPicking                 | PeakPicking.ProfilePeaks                                                                                                                  | Parameter to configure peak picking                                |

## Parameter Generation

Parameters are handed to the cli tool as a .json file. you can download the default file [here](https://github.com/CSBiology/ProteomIQon/blob/master/src/ProteomIQon/defaultParams/mzMLToMzLiteParams.json), 
or use an F# script, which can be downloaded or run in Binder at the top of the page, to write your own parameter file:
*)

#r "nuget: ProteomIQon, 0.0.1"

open ProteomIQon
open ProteomIQon.Domain

let defaultMzMLToMzLiteParams :Dto.PreprocessingParams =   
    {
        Compress                    = MzIO.Binary.BinaryDataCompressionType.NoCompression
        StartRetentionTime          = None
        EndRetentionTime            = None
        MS1PeakPicking              = PeakPicking.ProfilePeaks 
        MS2PeakPicking              = PeakPicking.ProfilePeaks 
    }


let serialized = 
    defaultMzMLToMzLiteParams
    |> Json.serialize

(***condition:ipynb***)
#if IPYNB
(**
If you are running this tool in Binder, you can copy the output of the following codeblock and save it in a JSON file.
*)
serialized
#endif // IPYNB

(**
## Executing the Tool
**Disclaimer** this tool converts [mzML](https://www.psidev.info/mzML) to [mzlite](https://github.com/CSBiology/MzIO/blob/developer/src/MzIO.SQL/MzIOSQL.fs) 
if you want to convert the other way around please visit the [MzLiteToMzML]({{root}}tools/MzLiteToMzML.html) documentation. 
We recommend the use of [msconvert](https://www.nature.com/articles/nbt.2377) to convert your raw data into .mzML. 
A user friendly way to execute msconvert is available through [Galaxy Europe](https://galaxyproject.eu/posts/2019/03/24/msconvert/).
To rescore all MS/MS to identify 'true' psms call: 

*)

(**
	proteomiqon-mzmltomzlite -i "path/to/your/run.mzML" -o "path/to/your/outDirectory" -p "path/to/your/params.json"
*)

(**
It is also possible to call the tool on a list of .mzML files. If you have a mulitcore cpu it is possible to convert multiple runs in parallel using the -c flag:
*)

(**
	proteomiqon-mzmltomzlite -i "path/to/your/run1.mzML" "path/to/your/run2.mzML" "path/to/your/run3.mzML" -o "path/to/your/outDirectory" -p "path/to/your/params.json" -c 3
*)

(**
A detailed description of the CLI arguments the tool expects can be obtained by calling the tool:
*)

(**
	proteomiqon-mzmltomzlite --help
*)

