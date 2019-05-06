// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"

open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain


let defaultCentroidizationParams = 

    let ms1CentroidizationParams  = 
        {
            NumberOfScales        = 3
            YThreshold            = 1.
            MzTolerance           = 0.1
            SNRS_Percentile       = 95.
            MinSNR                = 1.        
        }
    
    let ms2CentroidizationParams = 
        {
            /// Centroidization
            NumberOfScales          = 10
            Centroid_MzTolerance    = 0.1
            SNRS_Percentile         = 95.
            MinSNR                  = 1.
            /// Padding
            MaximumPaddingPoints    = Some 7
            Padding_MzTolerance     = 0.05
            WindowSize              = 150
            SpacingPerc             = 95.
        } 
         
    {
        Centroid                    = true
        UseManufacturerCentroids    = false
        StartRetentionTime          = nan
        EndRetentionTime            = nan
        MS1Centroidization          = ms1CentroidizationParams
        MS2CentroidizationParams    = ms2CentroidizationParams
    }


let serialized = 
    defaultCentroidizationParams
    |> Json.serialize

System.IO.File.WriteAllText( @"C:\Users\david\Source\Repos\netCoreRepos\ProteomIQon\src\ProteomIQon\defaultParams\centroidizationParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(@"C:\Users\david\Source\Repos\netCoreRepos\ProteomIQon\src\ProteomIQon\defaultParams\centroidizationParams.json")
    |> Json.deserialize<Dto.CentroidizationParams>
    |> CentroidizationParams.toDomain
