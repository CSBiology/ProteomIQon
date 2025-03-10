#r "netstandard"
#r @"../../../bin\MzMLToMzLite\net5.0\Newtonsoft.Json.dll"
#r @"../../../bin\alignmentbasedquantification\net5.0\FSharp.Stats.dll"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
#r @"../../../bin\alignmentbasedquantification\net5.0\BioFSharp.Mz.dll"
#r @"../../../bin\alignmentbasedquantification\net5.0\FSharpAux.IO.dll"

open FSharp.Stats
open FSharp.Stats.Signal
open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain
open ProteomIQon.FSharpStats'.Wavelet


let defaultPreprocessingParams :Dto.AlignmentBasedQuantificationParams = 
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
            MzWindow_Da                  = Domain.Window.Estimate
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
        PerformLocalWarp             = true
        XicExtraction                = XicExtraction
        BaseLineCorrection           = Some BaseLineCorrection
    }


let serialized = 
    defaultPreprocessingParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\AlignmentBasedQuantificationParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\AlignmentBasedQuantificationParams.json")
    |> Json.deserialize<Dto.QuantificationParams>
    |> QuantificationParams.toDomain




