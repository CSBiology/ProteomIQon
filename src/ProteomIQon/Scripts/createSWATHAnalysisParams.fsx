// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../packages\FSharp.Stats\lib\netstandard2.0\FSharp.Stats.dll"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
#r @"../../../packages\BioFSharp.Mz\lib\netstandard2.0\BioFSharp.Mz.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"

open FSharp.Stats
open FSharp.Stats.Signal
open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain
open ProteomIQon.FSharpStats'.Wavelet

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

let defaultSWATHAnalysisParams: Dto.SWATHAnalysisParams =
    {
        PeptideList          = None
        MatchingTolerancePPM = 100.
        QueryOffsetRange     = 10.
        XicProcessing        = Wavelet waveletParams
    }

let serialized = 
    defaultSWATHAnalysisParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\SWATHAnalysisParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\SWATHAnalysisParams.json")
    |> Json.deserialize<Dto.SWATHAnalysisParams>
    |> SWATHAnalysisParams.toDomain