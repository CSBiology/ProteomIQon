namespace ProteomIQon

open Domain
open Dto

module Json = 
  
  open Newtonsoft.Json
    
  let serialize obj =
    JsonConvert.SerializeObject obj

  let serializeAndWrite path obj =
    System.IO.File.WriteAllText(path,JsonConvert.SerializeObject obj)

  let deserialize<'a> str =
    JsonConvert.DeserializeObject<'a> str

  let ReadAndDeserialize<'a> path =
    System.IO.File.ReadAllText path
    |> JsonConvert.DeserializeObject<'a> 

