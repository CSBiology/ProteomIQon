(**
---
title: PSMBasedQuantification
category: Tools
categoryindex: 1
index: 4
---
*)
(**
[![Binder]({{root}}img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath={{fsdocs-source-basename}}.ipynb)&emsp;
[![Script]({{root}}img/badge-script.svg)]({{root}}{{fsdocs-source-basename}}.fsx)&emsp;
[![Notebook]({{root}}img/badge-notebook.svg)]({{root}}{{fsdocs-source-basename}}.ipynb)

# Peptide Spectrum Matching based Quantification
**Disclaimer** this tool needs a [peptide database]({{root}}tools/peptideDb.html) and [peptide spectrum matches]({{root}}tools/PeptideSpectrumMatching.html)
which [passed fdr thresholds]({{root}}tools/PSMStatistics.html).

Once a MS/MS spectrum is mapped to a peptide sequence the quantity of the fragmented peptide ion comes into view. 

Given an MS run in the mzLite or mzml format and a list of [fdr controlled peptide spectrum matches]({{root}}tools/PSMStatistics.html), 
this tool iterates accross all identified MS/MS scans and groups them by the assigned peptide ion. The scan times of each MS/MS spectrum 
are then weighted according to the quality of each match to build an reliable estimator for the scan time of the peptide ion in question.
This scan time estimator, combined with the monoisotopic m/z, is then used to extract an ion chromatogram. Using wavelet based peak detection techniques we identify
all peaks present in the XIC and select the most probable peak our target for quantification. Using parameter estimation techniques we subsequently use peak fitting
to fit a set of two gaussian models to the detected peak, from whom the one with the better fit is selected. This allows us not only to report how well the
signal fitted to the theoretical expected peak shape but also to obtain accurate estimates for the peak area, our estimator for peptide ion abundance.

<img src="{{root}}img/LabeledQuant.png" width="1000" height="750" />

The quantification tool was designed to allow label-free quantification as well as quantification of full metabolic labeled samples. For this we 
use the known identity of one of the the peptide ions and calculate the m/z of the unobserved differentially labeled counterpart to extract and 
quantify the corresponding XIC. 

## Parameters
The following table gives an overview of the parameter set:

| **Parameter**                  | **Default Value**                                                                                                                         | **Description**                                                    |
|--------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------|
| PerformLabeledQuantification   | true                                                                                                                                      | Indicates if a labeled quantification should be performed          |
| XicExtraction                  | {TopKPSMs = None; ScanTimeWindow = 2.; MzWindow_Da = 0.07; XicProcessing = Wavelet waveletParams}                                         | Parameters to tune the Xic extraction                              |
| BaseLineCorrection             | Some { MaxIterations = 10; Lambda = 6; P = 0.05 }                                                                                         | optional parameter, defines if and how baseline correction should be performed |

        
## Parameter Generation

Parameters are handed to the cli tool as a .json file. you can download the default file [here](https://github.com/CSBiology/ProteomIQon/blob/master/src/ProteomIQon/defaultParams/peptideSpectrumMatchingParams.json), 
or use an F# script, which can be downloaded or run in Binder at the top of the page, to write your own parameter file:
*)

#r "nuget: BioFSharp.Mz, 0.1.5-beta"
#r "nuget: ProteomIQon, 0.0.1"

open ProteomIQon
open ProteomIQon.Domain
open FSharp.Stats.Signal

open BioFSharp.Mz

let defaultQuantificationParams :Dto.QuantificationParams = 
    ///
    let waveletParams :WaveletParameters = 
        {
            Borderpadding           = None    
            BorderPadMethod         = Padding.BorderPaddingMethod.Random 
            InternalPaddingMethod   = Padding.InternalPaddingMethod.LinearInterpolation 
            HugeGapPaddingMethod    = Padding.HugeGapPaddingMethod.Zero
            HugeGapPaddingDistance  = 100.
            MinPeakDistance         = None
            MinPeakLength           = Some 0.1
            MaxPeakLength           = 1.5 
            NoiseQuantile           = 0.01 
            MinSNR                  = 0.01  
        }

    let XicExtraction = 
        {
            TopKPSMs                     = None
            ScanTimeWindow               = 2.
            MzWindow_Da                  = 0.07 
            XicProcessing                = Wavelet waveletParams
        }

    let BaseLineCorrection = 
        {
            MaxIterations                = 10 
            Lambda                       = 6 
            P                            = 0.05
        }
    {
        PerformLabeledQuantification = true
        XicExtraction                = XicExtraction
        BaseLineCorrection           = Some BaseLineCorrection
    }


let serialized = 
    defaultQuantificationParams
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
**Disclaimer** this tool needs a [peptide database]({{root}}tools/peptideDb.html) and [peptide spectrum matches]({{root}}tools/PeptideSpectrumMatching.html)
which [passed fdr thresholds]({{root}}tools/PSMStatistics.html).

*)

(**
	proteomiqon-psmbasedquantification -i "path/to/your/run.mzml" -ii "path/to/your/run.spsm" -d "path/to/your/database.sqlite" -o "path/to/your/outDirectory" -p "path/to/your/params.json"
*)

(**
It is also possible to call the tool on a lists of MS and scored psm files. If you have a mulitcore cpu it is possible to score multiple runs in parallel using the -c flag:
*)

(**
	proteomiqon-psmbasedquantification -i "path/to/your/run1.mzml" "path/to/your/run2.mzml" "path/to/your/run3.mzml" -ii "path/to/your/run1.spsm" "path/to/your/run2.spsm" "path/to/your/run3.spsm" -d "path/to/your/database.sqlite" -o "path/to/your/outDirectory" -p "path/to/your/params.json" -c 3
*)

(**
A detailed description of the CLI arguments the tool expects can be obtained by calling the tool:
*)

(**
	proteomiqon-psmbasedquantification --help
*)
