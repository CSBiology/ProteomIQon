namespace ProteomIQon

open Domain
open Core 
open System.IO
open BioFSharp.Mz

module PeptideDB = 

    let createPeptideDB (processParams:PeptideDBParams) (outputDir:string) =
        
        let searchDBParams :SearchDB.SearchDbParams = 
            {        
            Name                = processParams.Name                
            DbFolder            = outputDir           
            FastaPath           = processParams.FastaPath           
            FastaHeaderToName   = processParams.FastaHeaderToName   
            Protease            = processParams.Protease            
            MinMissedCleavages  = processParams.MinMissedCleavages  
            MaxMissedCleavages  = processParams.MaxMissedCleavages  
            MaxMass             = processParams.MaxMass             
            MinPepLength        = processParams.MinPepLength        
            MaxPepLength        = processParams.MaxPepLength        
            IsotopicMod         = processParams.IsotopicMod          
            MassMode            = processParams.MassMode            
            MassFunction        = processParams.MassFunction        
            FixedMods           = processParams.FixedMods                       
            VariableMods        = processParams.VariableMods        
            VarModThreshold     = processParams.VarModThreshold     
            } 
        printfn "Initializing PeptideDB at: %s" outputDir
        let db = SearchDB.connectOrCreateDB searchDBParams
        printfn "Successfully created PeptideDB"
        printfn "Done"
