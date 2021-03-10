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

BioFSharp.Mass.toMZ 1000. 1.
BioFSharp.Mass.deltaMassByPpm 1000. 800.
//// 
let optimizeWindowWidth polOrder (windowWidthToTest:int[]) noiseAutoCorr (signalOfInterest:float[]) =
    let signalOfInterest' = signalOfInterest |> vector
    //let noiseAutoCorr = Correlation.Vector.autoCorrelation 1 (blankSignal |> vector)
    let filterF w yData = FSharp.Stats.Signal.Filtering.savitzky_golay w polOrder 0 0 yData
    let windowWidthToTest' = windowWidthToTest |> Array.filter (fun x -> x%2 <> 0)
    let optimizedWindowWidth =
        windowWidthToTest'
        |> Array.map (fun w ->
                      let smoothedY = filterF w signalOfInterest
                      let noise = smoothedY - (signalOfInterest')
                      w, Correlation.Vector.autoCorrelation 1 noise
                     )
        |> Array.minBy (fun (w,ac) -> (ac - noiseAutoCorr) |> abs )
        |> fst
    optimizedWindowWidth

optimizeWindowWidth 2 [|5 .. 2 .. 60|] 0.1 [|0. .. 100.|]
let divideAndOdd polOrder x =
    let tmp = 
        let x = x / 2
        if x % 2 = 0 then (x - 1) else x
    if tmp <= polOrder + 1 && polOrder % 2 = 0 then polOrder + 3
    elif tmp <= polOrder + 1 && polOrder % 2 <> 0 then polOrder + 4
    else tmp 

divideAndOdd 2 40

//let xData= [|23.508557; 23.529638; 23.54829; 23.587451; 23.608678; 23.621234; 23.64388;
//23.653768; 23.662876; 23.676972; 23.684484; 23.692049; 23.701103; 23.708632;
//23.726701; 23.732395; 23.750797; 23.759155; 23.770524; 23.777755; 23.792884;
//23.799155; 23.80668; 23.819814; 23.825997; 23.839103; 23.842822; 23.848453;
//23.852064; 23.859514; 23.863299|]
//let yData= [|0.0; 0.0; 46836.5431; 0.0; 0.0; 960794.4469; 0.0; 0.0; 1021751.182; 0.0; 0.0;
//0.0; 1274490.831; 1644802.612; 0.0; 1805123.74; 0.0; 1217825.244; 1313588.049;
//1094465.574; 0.0; 0.0; 934788.4569; 0.0; 1445439.839; 0.0; 0.0; 0.0; 0.0; 0.0;
//0.0|]
let xData= [|17.014015; 17.034184; 17.053365; 17.073361; 17.092534; 17.11459; 17.133842;
17.153949; 17.174103; 17.21274; 17.232397|]
//upperBound: vector [|581232.0881; 16.52613333; 0.04556124508; nan|]
let yData= [|83292.01773; 41548.27553; 79729.07472; 104574.355; 203244.4582; 198756.6547;
225618.1578; 127809.2722; 122217.9089; 0.0; 0.0|]

//let interpolate =
//FSharp.Stats.Interpolation.LinearSpline.initInterpolate
open Plotly.NET
let InitialParamGuess= [|1805123.74; 23.73800891; 0.05807036296; nan|]
let lowerBound= vector [|1263586.618; 23.56063791; 0.02903518148; nan|]
let upperBound= vector [|2346660.862; 23.91537991; 0.08710554444; nan|]
[|5; 7; 9; 11; 13; 15; 17; 19; 21; 23; 25; 27; 29; 31; 33; 35; 37; 39; 41;
43; 45; 47; 49; 51; 53; 55; 57; 59|]
|> Array.filter (fun x -> x%2 <> 0)
Chart.Point (xData,yData)
|> Chart.Show
///
let v = Quantification.ParameterEstimation.varianceOf xData yData
let s = Quantification.ParameterEstimation.skewOf xData yData

let a = [|0. .. 10.|]
let b = [|0.;0.;2.;4.;5.;4.5;4.;3.;2.;1.;0.|]   

Chart.Point (a,b)
|> Chart.Show

let v = Quantification.ParameterEstimation.varianceOf a b
let s = Quantification.ParameterEstimation.skewOf a b

///    
let estTau stdev skew  =
    stdev * (skew/2.) (*** 0.333*)

Quantification.ParameterEstimation.estTau (sqrt(v)) s
estTau (sqrt(v)) s

    
    
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
        MatchingIonTolerancePPM = 100.       
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


System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\peptideSpectrumMatchingParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\peptideSpectrumMatchingParams.json")
    |> Json.deserialize<Dto.PeptideSpectrumMatchingParams>
    |> PeptideSpectrumMatchingParams.toDomain
