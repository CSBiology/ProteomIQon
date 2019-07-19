// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../packages\BioFSharp\lib\netstandard2.0\BioFSharp.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"
#r @"../../../packages\BioFSharp.Mz\lib\netstandard2.0\BioFSharp.Mz.dll"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"

open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain
open BioFSharp.Mz

BioFSharp.Mass.deltaMassByPpm 1000. 800.

let defaultPreprocessingParams :Dto.PeptideSpectrumMatchingParams = 

    let chargeDetermParams :ChargeState.ChargeDetermParams = 
        {
        ExpectedMinimalCharge   = 2
        ExpectedMaximumCharge   = 5
        Width                   = 1.1
        MinIntensity            = 0.15
        DeltaMinIntensity       = 0.3
        NrOfRndSpectra          = 10000
        }

    let andromedaParams = 
        {
        PMinPMax                = 4,10
        MatchingIonTolerancePPM = 1000.       
        }
    {
        ChargeStateDeterminationParams  = chargeDetermParams 
        LookUpPPM                       = 30.
        MS2ScanRange                    = 100.,2000.
        nTerminalSeries                 = NTerminalSeries.B
        cTerminalSeries                 = CTerminalSeries.Y
        Andromeda                       = andromedaParams
    }


let serialized = 
    defaultPreprocessingParams
    |> Json.serialize


System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\peptideSpectrumMatchingParamsThermo.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\peptideSpectrumMatchingParams.json")
    |> Json.deserialize<Dto.PeptideSpectrumMatchingParams>
    |> PeptideSpectrumMatchingParams.toDomain
