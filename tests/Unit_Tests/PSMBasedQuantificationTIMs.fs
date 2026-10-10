module ProteomIQon.UnitTests

open Expecto 
open ProteomIQon
open ProteomIQon.PSMBasedQuantificationTIMs
open MzIO.Processing
open MzIO.Binary
open MzIO.Processing.MzIOLinq
open BioFSharp.Mz

let private peak intensity rt mobility = 
    Peak1D(intensity, 500.0, mobility)
    |> fun p -> RtIndexEntry.AsPeak2D(p, rt)

[<Tests>]
let gaborGridTest = 
    testList "tryCreateGaborGrid" [
        testCase "aggregates peaks in the same bin and ignores invalid peaks" <| fun _ -> 
            let input = 
                [|
                    [|
                        peak 2.0 10.0000 1.0000
                        peak 3.0 10.0000 1.0004
                        peak System.Double.NaN 10.0080 1.0000
                        peak 0.0 10.0120 1.0000
                    |]
                |]
            let grid = 
                Expect.wantSome
                    (tryCreateGaborGrid input)
                    "Grid should have been created"
            
            Expect.equal grid.RtBinCount 1 "Unexpected RT bin count"
            Expect.equal grid.MobilityBinCount 1 "Unexpected Ionmobility bin count"
            Expect.equal grid.Intensities.[0].[0] 5.0 "Peaks are not aggregated correctly"
    ]

[<Tests>]
let qpsmHeaderTest =
    testList "QpsmInput.checkHeader" [
        testCase "accepts the PSMStatistics layout and rejects the former FragPipe layout" <| fun _ ->
            let withHeader (header: string) f =
                let path = System.IO.Path.GetTempFileName()
                try
                    System.IO.File.WriteAllText(path, header + "\n")
                    f path
                finally
                    System.IO.File.Delete path
            withHeader (String.concat "\t" QpsmInput.expectedColumns) QpsmInput.checkHeader
            let fragPipeHeader =
                "PSMId\tPepSequenceID\tModSequenceID\tScanTime\tCharge\tPrecursorMZ\tTheoMass\tAbsDeltaMass\tIonMobility\tPeptideLength\tMissCleavages\tExpectscore\tHyperscore\tProbability\tStringSequence\tProteinNames\tGlobalMod"
            withHeader fragPipeHeader (fun path ->
                Expect.throws (fun () -> QpsmInput.checkHeader path) "the FragPipe layout is rejected")
    ]
[<Tests>]
let qpsmEmptyTest =
    testList "QpsmInput.read" [
        testCase "an empty file and a file with the header only have no PSMs" <| fun _ ->
            for content in [ ""; String.concat "\t" QpsmInput.expectedColumns + "\n" ] do
                let path = System.IO.Path.GetTempFileName()
                try
                    System.IO.File.WriteAllText(path, content)
                    Expect.isEmpty (QpsmInput.read path) "no PSMs"
                finally
                    System.IO.File.Delete path
    ]
