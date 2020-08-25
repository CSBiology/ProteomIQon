// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\net47\System.Data.Sqlite.dll"
#r @"../../../bin\ProteomIQon\net47\Newtonsoft.Json.dll"
#r @"../../../bin\ProteomIQon\net47\FSharp.Stats.dll"
                              

#r @"../../../bin\ProteomIQon\net47\FSharp.Stats.dll"
#r @"../../../bin\ProteomIQon\net47\ProteomIQon.dll"
#r @"../../../bin\ProteomIQon\net47\FSharp.Plotly.dll"

#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"
#r @"../../../bin\ProteomIQon\net47\BioFSharp.dll"
#r @"../../../bin\ProteomIQon\net47\BioFSharp.Mz.dll"


open ProteomIQon
open ProteomIQon.Domain
open ProteomIQon.Dto
open FSharp.Plotly

let res = DTW'.align DTW'.s1 DTW'.s2
let res2 = DTW'.align' DTW'.s1 DTW'.s2 46.

[
Chart.Point(DTW'.s1)
Chart.Point(DTW'.s2)
Chart.Point(res)
Chart.Point([snd res2,2.])
]
|> Chart.Combine
|> Chart.Show
open FSharp.Stats

let xV = [|1. .. 10000.|]                            
let yV = [|1. .. 10000.|]
let xy = Array.zip xV yV
10000 / 200




let weightedMean (weights:seq<'T>) (items:seq<'T>) =
    let sum,n = Seq.fold2 (fun (sum,n) w i -> w*i+sum,n + w ) (0.,0.) weights items 
    sum / n 

let x = [1.;2.;3.;3.;3.]
weightedMean [1.;1.;1.;1.;0.001] x
open FSharp.Stats
FSharp.Stats.Correlation.Vector.autoCorrelation 1 ([|0.;0.|] |> vector)
([1.;2.;3.;3.;3.] |> vector)

mzliteReader.ReadMassSpectra("sample=0") 

type dummy = 
    {
    Name : string
    nr : int
    nrf : float
    nrf1 : float
    nrf2 : float
    nrf3 : float
    nrf4 : float
    nrf5 : float
    nrf6 : float
    testArr: int []
    lul:string 
    }


let dat = 
    [|
        for i = 0 to 720000 do 
             
            yield {
                Name    = (i.ToString())
                nr      = i
                nrf     = FSharp.Stats.Distributions.Continuous.Normal.Sample 10. 2.
                nrf1    = FSharp.Stats.Distributions.Continuous.Normal.Sample 10. 2.
                nrf2    = FSharp.Stats.Distributions.Continuous.Normal.Sample 10. 2.
                nrf3    = FSharp.Stats.Distributions.Continuous.Normal.Sample 10. 2.
                nrf4    = FSharp.Stats.Distributions.Continuous.Normal.Sample 10. 2.
                nrf5    = FSharp.Stats.Distributions.Continuous.Normal.Sample 10. 2.
                nrf6    = FSharp.Stats.Distributions.Continuous.Normal.Sample 10. 2.
                testArr = [|1..10|]
                lul     = "hallloooooo;hallloooooo;hallloooooo;hallloooooo;hallloooooo;hallloooooo;hallloooooo;hallloooooo;hallloooooo" 
            }
    |]
                    
#time

let y = 
    dat 
    |> FSharpAux.IO.SeqIO.Seq.toCSV "\t" true
    |> Array.ofSeq

y
|> Seq.map (fun x -> FSharpAux.String.replace ";" "\t" x)
|> Array.ofSeq

|> FSharpAux.IO.FileIO.writeToFile false @"C:\Users\david\source\repos\test.txt"

dat 
|> Array.take 10
|> FSharpAux.IO.SeqIO.Seq.toCSV "\t" true
|> Seq.map (fun x -> FSharpAux.String.replace ";" "\t" x)
|> FSharpAux.IO.FileIO.writeToFile false @"C:\Users\david\source\repos\test4.txt"

let restorePSMID psmID =
    FSharpAux.String.split '_' psmID  
    |> Array.item 0
    |> FSharpAux.String.split '-'
    |> String.concat " "

restorePSMID "scan=104_87_2_4"

restorePSMID "sample=0-experiment=1-scan=2780_1579_2_0"

"Raid Area 51 plx"
|> Seq.filter (System.Char.IsWhiteSpace >> not)
|> Seq.sortBy (fun char -> System.Char.ToLower char)
|> System.String.Concat

//System.IO.FileInfo(assembly.Location).DirectoryName + @"\percolator-v3-01\bin\percolator.exe"




open BioFSharp.ModificationInfo
open BioFSharp.AminoAcids

let pyro_Glu'GluNterm' =
    createSearchModification "Pyro_Glu'Gln'" "28" "	Pyro-glu from Q" true "H3N"
        [Specific(Gln,ModLocation.Nterm);] SearchModType.Plus "pq"

let modi = BioFSharp.Mz.SearchDB.ModCombinator.convertSearchModification []  pyro_Glu'GluNterm' 

open BioFSharp.Formula
open BioFSharp
open BioFSharp.Mz
let massF:(IBioItem->float) = BioItem.initMonoisoMassWithMem


let BioSeq = "EGLDVHFVDEYEK" |> BioList.ofAminoAcidString
let BioSeqH = "EGLDVHFVDEYEK" |> BioList.ofAminoAcidString |> List.map (fun x ->AminoAcids.setModifications [ModificationInfo.Table.N15] x)

let res = Fragmentation.Series.yOfBioList massF BioSeq
let res2 = Fragmentation.Series.yOfBioList massF BioSeqH

#time
for i = 1 to 100000000 do "EGLDVHFVDEYEK" |> BioList.ofAminoAcidString |> List.map (fun x ->AminoAcids.setModifications [ModificationInfo.Table.N15] x)

let light = res.[0].MainPeak.Mass

let heavy = res2.[0].MainPeak.Mass

//let minus = (Formula.add (Formula.parseFormulaString "H3N")) Formula.emptyFormula  

//let plus  = (Formula.substract (Formula.parseFormulaString "H3N")) Formula.emptyFormula   

//let p     = (Formula.substract (Formula.parseFormulaString "H3N"))  (Formula.parseFormulaString "H1N") 
