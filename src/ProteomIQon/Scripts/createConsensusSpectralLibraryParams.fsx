// Learn more about F# at http://fsharp.org. See the 'F# Tutorial' project
// for more guidance on F# programming.
#r "netstandard"
#r @"../../../bin\ProteomIQon\netstandard2.0\ProteomIQon.dll"

open ProteomIQon
open ProteomIQon.Dto

let defaultConsensusSpectralLibraryParams: Dto.ConsensusSpectralLibraryParams =
    {
        // InitialPeptideSelection
        BinningWindowWidth                          = 10.
        FractionOfMostAbundandIonsPerBin            = 0.1
        MinFragmentCount                            = 3
        MinFragmentLadderIdx                        = 3
        MinPeptideLength                            = 8
        // XicExtraction                            = // XicExtraction
        RtWindowWidth                               = 10.                               
        // Matching                                 = // Matching
        FragMatchingBinWidth                        = 0.01                        
        FragMatchingBinOffset                       = 0.0
        MS2ScanRange                                = 100.,2000.
        // Filtering                                = // Filtering
        MaxRatioMostAbundandVsSecondAbundandPeak    = 0.2
    }

let serialized = 
    defaultConsensusSpectralLibraryParams
    |> Json.serialize

System.IO.File.WriteAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\ConsensusSpectralLibraryParams.json",serialized)

let deserialized = 
    System.IO.File.ReadAllText(__SOURCE_DIRECTORY__ + @"/../defaultParams\ConsensusSpectralLibraryParams.json")
    |> Json.deserialize<Dto.ConsensusSpectralLibraryParams>
    |> ConsensusSpectralLibraryParams.toDomain