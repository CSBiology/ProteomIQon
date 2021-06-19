namespace ProteomIQon

open Domain
open Core
open System.IO
open BioFSharp.Mz

module PeptideDB =

    let createPeptideDB (processParams:PeptideDBParams) (outputDir:string) (fastaPath:string)=

        let logger = Logging.createLogger (sprintf "PeptideDB_%s" processParams.Name)

        let searchDBParams :SearchDB.SearchDbParams =
            {
            Name                = processParams.Name
            DbFolder            = outputDir
            FastaPath           = fastaPath
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
        let dbConnection = SearchDB.connectOrCreateDB searchDBParams
        logger.Trace "Successfully created PeptideDB"
        logger.Trace "Set Index on data base if not present."
        SearchDB'.setIndexOnModSequenceAndGlobalMod dbConnection |> ignore
        logger.Trace "Set Index on data base if not present: finished"
        dbConnection.Dispose()
        logger.Trace "Done"