// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
#r @"../../../packages\BioFSharp.Mz\lib\netstandard2.0\BioFSharp.Mz.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"

open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain

let defaultPreprocessingParams :Dto.QuantificationParams = 

    let XicExtraction = 
        {
            ScanTimeWindow               = 2.
            MzWindow_Da                  = 0.07 
            MinSNR                       = 0.1  
            PolynomOrder                 = 2
            WindowSize                   = Domain.WindowSize.Estimate
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
    defaultPreprocessingParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\QuantificationParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\QuantificationParams.json")
    |> Json.deserialize<Dto.QuantificationParams>
    |> QuantificationParams.toDomain
