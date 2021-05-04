(**
// can't yet format YamlFrontmatter (["title: PSMBasedQuantification"; "category: Tools"; "categoryindex: 1"; "index: 6"], Some { StartLine = 2 StartColumn = 0 EndLine = 6 EndColumn = 8 }) to pynb markdown

[![Binder](https://csbiology.github.io/ProteomIQon/img/badge-binder.svg)](https://mybinder.org/v2/gh/csbiology/ProteomIQon/gh-pages?filepath=tools/PSMBasedQuantification.ipynb)&emsp;
[![Script](https://csbiology.github.io/ProteomIQon/img/badge-script.svg)](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.fsx)&emsp;
[![Notebook](https://csbiology.github.io/ProteomIQon/img/badge-notebook.svg)](https://csbiology.github.io/ProteomIQon/tools/PSMBasedQuantification.ipynb)

# Peptide Spectrum Matching based Quantification
**Disclaimer** this tool needs a [peptide database](https://csbiology.github.io/ProteomIQon/tools/peptideDb.html) and [peptide spectrum matches](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.html)
which [passed fdr thresholds](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html).

Once a MS/MS spectrum is mapped to a peptide sequence the quantity of the fragmented peptide ion comes into view. 

Given an MS run in the mzLite or mzml format and a list of [fdr controlled peptide spectrum matches](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html), 
this tool iterates accross all identified MS/MS scans and groups them by the assigned peptide ion. The scan times of each MS/MS spectrum 
are then weighted according to the quality of each match to build an reliable estimator for the scan time of the peptide ion in question.
This scan time estimator, combined with the monoisotopic m/z, is then used to extract an ion chromatogram. Using wavelet based peak detection techniques we identify
all peaks present in the XIC and select the most probable peak our target for quantification. Using parameter estimation techniques we subsequently use peak fitting
to fit a set of two gaussian models to the detected peak, from whom the one with the better fit is selected. This allows us not only to report how well the
signal fitted to the theoretical expected peak shape but also to obtain accurate estimates for the peak area, our estimator for peptide ion abundance.

<img src="https://csbiology.github.io/ProteomIQon/img/LabeledQuant.png" width="1000" height="750" />
<img src="https://csbiology.github.io/ProteomIQon/img/LabeledQuant.png" width="1000" height="750" />

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
(**
## Executing the Tool
**Disclaimer** this tool needs a [peptide database](https://csbiology.github.io/ProteomIQon/tools/peptideDb.html) and [peptide spectrum matches](https://csbiology.github.io/ProteomIQon/tools/PeptideSpectrumMatching.html)
which [passed fdr thresholds](https://csbiology.github.io/ProteomIQon/tools/PSMStatistics.html).


	proteomiqon-psmbasedquantification -i "path/to/your/run.mzlite" -ii "path/to/your/run.qpsm" -d "path/to/your/database.sqlite" -o "path/to/your/outDirectory" -p "path/to/your/params.json"

It is also possible to call the tool on a lists of MS and scored psm files. If you have a mulitcore cpu it is possible to score multiple runs in parallel using the -c flag:

	proteomiqon-psmbasedquantification -i "path/to/your/run1.mzlite" "path/to/your/run2.mzlite" "path/to/your/run3.mzlite" -ii "path/to/your/run1.qpsm" "path/to/your/run2.qpsm" "path/to/your/run3.qpsm" -d "path/to/your/database.sqlite" -o "path/to/your/outDirectory" -p "path/to/your/params.json" -c 3

A detailed description of the CLI arguments the tool expects can be obtained by calling the tool:

	proteomiqon-psmbasedquantification --help

*)