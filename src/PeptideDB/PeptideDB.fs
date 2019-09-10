namespace ProteomIQon

open Domain
open Core
open System.IO
open BioFSharp.Mz

module PeptideDB =

    let createPeptideDB (processParams:PeptideDBParams) (outputDir:string) =

        let logger = Logging.createLogger (sprintf @"%s\PeptideDB_log.txt" outputDir) "PeptideDB_createPeptideDB"

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
        logger.Trace (sprintf "Initializing PeptideDB at: %s" outputDir)
        logger.Trace (sprintf "searchDBParams: %A" searchDBParams)
        let db = SearchDB.connectOrCreateDB searchDBParams
        logger.Trace "Successfully created PeptideDB"
        logger.Trace "Done"