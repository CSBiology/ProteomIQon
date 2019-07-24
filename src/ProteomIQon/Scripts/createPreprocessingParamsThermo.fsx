// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"

#r @"../../../packages\BioFSharp\lib\netstandard2.0\BioFSharp.dll"
#r @"../../../packages\BioFSharp.Mz\lib\netstandard2.0\BioFSharp.Mz.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"

open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain
//[0..19] 
//|> List.splitInto 16
//|> List.length
BioFSharp.Mass.deltaMassByPpm 1000. 700.

let defaultPreprocessingParams :Dto.PreprocessingParams = 

    let ms1PeakPickingParams  = 
        {
            
            NumberOfScales          = 3
            YThreshold              = YThreshold.Fixed 1.
            Centroid_MzTolerance    = 0.1
            SNRS_Percentile         = 95.
            MinSNR                  = 1.
            PaddingParams           = None
        }

         
    {
        Compress                    = true
        StartRetentionTime          = None
        EndRetentionTime            = None 
        MS1PeakPicking              = PeakPicking.Centroid (CentroidizationMode.Wavelet ms1PeakPickingParams)
        MS2PeakPicking              = PeakPicking.ProfilePeaks //PeakPicking.Centroid (CentroidizationMode.Wavelet ms2PeakPickingParams)
    }


let serialized = 
    defaultPreprocessingParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\preprocessingParamsThermo.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\preprocessingParamsThermo.json")
    |> Json.deserialize<Dto.PreprocessingParams>
    |> PreprocessingParams.toDomain
