// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\net47\System.Data.Sqlite.dll"
#r @"../../../bin\ProteomIQon\net47\Newtonsoft.Json.dll"
                              
#r @"../../../bin\ProteomIQon\net47\MzLite.dll"
#r @"../../../bin\ProteomIQon\net47\MzLite.SQL.dll"
#r @"../../../bin\ProteomIQon\net47\MzLite.Processing.dll"
#r @"../../../bin\ProteomIQon\net47\ProteomIQon.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.IO.dll"
#r @"../../../packages\FSharpAux.IO\lib\netstandard2.0\FSharpAux.dll"


open ProteomIQon.Core

let mzliteReader = MzLite.Reader.getReader @"C:\Users\david\Source\Repos\netCoreRepos\ProteomiconTest\20171129 FW LW tot001.mzlite"
let tr = mzliteReader.BeginTransaction()

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


mzliteReader.ReadMassSpectra("sample=0") 



let x = [|1. ..10000.|]
let y = [|1. .. 0.1 .. 10000.|]
y .Length
y
|> Array.map (fun yy -> x |> Array.findBack (fun x -> x <= yy),yy)
|> Array.groupBy fst
|> Array.map (fun (ms1,ms1And2) -> ms1,ms1And2 |> Array.map snd)


open ProteomIQon
open ProteomIQon.Dto

open ProteomIQon.Json
open ProteomIQon.Dto


let deserialized = 
    System.IO.File.ReadAllText(@"C:\Users\david\Source\Repos\netCoreRepos\ProteomiconTest\preprocessingParams.json")
    |> Json.deserialize<Dto.PreprocessingParams>
    |> PreprocessingParams.toDomain





