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