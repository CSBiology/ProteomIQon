// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\peptidespectrummatching\net5.0\BioFSharp.dll"
#r @"../../../bin\peptidespectrummatching\net5.0\FSharpAux.IO.dll"
#r @"../../../bin\peptidespectrummatching\net5.0\BioFSharp.Mz.dll"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"
#r @"../../../bin\peptidespectrummatching\net5.0\FSharp.Stats.dll"
#r @"../../../bin\peptidespectrummatching\net5.0\Plotly.NET.dll"


open ProteomIQon
open ProteomIQon.Dto
open ProteomIQon.Domain
open BioFSharp.Mz
open FSharp.Stats

    
let defaultPSMParams :Dto.PeptideSpectrumMatchingParams = 

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
        MatchingIonTolerancePPM = 100.       
        }
    {
        ChargeStateDeterminationParams  = chargeDetermParams 
        LookUpPPM                       = 30.
        nTerminalSeries                 = NTerminalSeries.B
        cTerminalSeries                 = CTerminalSeries.Y
        Andromeda                       = andromedaParams
    }


let serialized = 
    defaultPSMParams
    |> Json.serialize


System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\peptideSpectrumMatchingParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\peptideSpectrumMatchingParams.json")
    |> Json.deserialize<Dto.PeptideSpectrumMatchingParams>
    |> PeptideSpectrumMatchingParams.toDomain
