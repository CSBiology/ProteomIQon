// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\net47\System.Data.Sqlite.dll"
#r @"../../../bin\ProteomIQon\net47\Newtonsoft.Json.dll"
#r @"../../../bin\ProteomIQon\net47\FSharp.Stats.dll"
                              
#r @"../../../bin\ProteomIQon\net47\MzLite.dll"
#r @"../../../bin\ProteomIQon\net47\MzLite.SQL.dll"
#r @"../../../bin\ProteomIQon\net47\MzLite.Processing.dll"
#r @"../../../bin\ProteomIQon\net47\ProteomIQon.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.dll"


open ProteomIQon.Core

let mzliteReader = MzLite.Reader.getReader @"C:\Users\david\Source\Repos\netCoreRepos\ProteomiconTest\20171129 FW LW tot001.mzlite"
let tr = mzliteReader.BeginTransaction()
open System
open System.IO
Directory.Exists @"Z:\5_6BackUpHyperCSB2\SFB-core\Acclimation_3hm\181005_cold1_3h_GD1_01_8993.d"
Directory.GetFiles( @"Z:\5_6BackUpHyperCSB2\SFB-core\Acclimation_3hm\181005_cold1_3h_GD1_01_8993.d",("*.d")) 
Directory.GetDirectories(@"Z:\5_6BackUpHyperCSB2\SFB-core\Acclimation_3hm\",("*.d"))

let spec = 
    mzliteReader.ReadMassSpectra("sample=0")
    |> Array.ofSeq
    |> Array.filter (fun x -> ProteomIQon.Core.MzLite.MassSpectrum.getMsLevel x = 2 )
true || true
let length = 
    spec 
    |> Array.map (fun spec -> mzliteReader.ReadSpectrumPeaks(spec.ID).Peaks.Length)

length .Length
length 
|> Array.filter (fun x -> x <> 0)
|> Array.length

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