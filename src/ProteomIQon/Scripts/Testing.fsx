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

let mzliteReader = MzLite.Reader.getReader @"C:\Users\david\Source\Repos\netCoreRepos\ProteomiconTest\20171129 FW LW tot001.wiff"
let tr = mzliteReader.BeginTransaction()

mzliteReader.ReadMassSpectrum("sample=0 experiment=0 scan=3")

mzliteReader.ReadMassSpectra("sample=0") 












