module ProteomIQon.UnitTestsPeptideSpectrumMatchingTIMsRegressions

open System
open Expecto
open ProteomIQon
open PeptideIndex
open SpectrumSearch
open PeptideSpectrumMatchingTIMs
open BioFSharp
open BioFSharp.Mz
open ProteomIQon.Core

let settings =
    { PrecursorTolerancePpm=20.; FragmentTolerancePpm=20.; IsotopeErrors=[|0;1;2|]
      MaxFragmentCharge=2; TopNPeaks=150; MinimumRatio=0.01; RemovePrecursorRange=1.5
      Deisotope=true; MinimumPeaks=15; MinMatchedFragments=4; MinFragmentsModelling=1; ReportedHitsPerLabel=10 }

// Twenty-one indistinguishable binned candidates, but only the last matches exactly.
let table =
    let residue rank j =
        let first = if rank=20 then 100.009 else 100.0001 + float rank*0.00005
        match j with 0 -> first | 1 -> 50. | 2 -> 60. | 3 -> 70. | _ -> 201.-first
    { Count=21; Mass=Array.create 21 (381.+waterMass); Order=Array.init 21 id
      ModSequenceIdRow=Array.init 21 id; PepSequenceIdRow=Array.init 21 id; GlobalModRow=Array.zeroCreate 21
      SeqStartRow=Array.init 22 (fun i -> 5*i); Residues=Array.init 105 byte
      Codes=Array.init 105 (fun i -> {Token=string i; Masses=[|residue (i/5) (i%5)|]})
      MaxLength=5; RawSequenceRow=dict [] }

let peaks =
    { Mz=Array.append ([|100.009;150.009;210.009;280.009|] |> Array.map ((+) protonMass)) (Array.init 11 (fun i -> 1000.+float i*10.))
      Intensity=Array.create 15 100. }

let pass lo hi =
    let index=buildFragmentIndex table 0.01 500. 1 lo hi ignore
    match indexPass settings table index (Scratch()) 10 (table.Mass.[0]/2.+protonMass) peaks 2 lo hi with
    | SliceFound (scores,t,d) -> scores,t,d
    | x -> failwithf "Unexpected %A" x

[<Tests>]
let all = testList "PeptideSpectrumMatchingTIMs regressions" [
    testCase "splitting one mz intensity across mobility scans cannot change the integrated search spectrum" <| fun _ ->
        let split = global.MzIO.Commons.Arrays.ArrayWrapper [|
            global.MzIO.Binary.Peak1D(55., 457., 0.73)
            global.MzIO.Binary.Peak1D(57., 457., 0.74)
            global.MzIO.Binary.Peak1D(102., 457., 0.75)
            global.MzIO.Binary.Peak1D(150., 600., 0.74) |]
        let combined = global.MzIO.Commons.Arrays.ArrayWrapper [|
            global.MzIO.Binary.Peak1D(214., 457., 0.74)
            global.MzIO.Binary.Peak1D(150., 600., 0.74) |]
        let select (mz,intensity) = preprocess {settings with TopNPeaks=1;Deisotope=false} 900. 2 mz intensity
        let splitResult = searchPeakArraysWithTolerance 0. split |> select
        let combinedResult = searchPeakArraysWithTolerance 0. combined |> select
        Expect.equal splitResult combinedResult "the same integrated signal must produce the same normalized top peak"
        // The former coordinate-dropping reader fails this physical invariance: it picks
        // 600 instead of 457 when the larger total intensity is split over mobility scans.
        let oldSplit = MzIO.Peaks.unzipIMzliteArray split |> select
        let oldCombined = MzIO.Peaks.unzipIMzliteArray combined |> select
        Expect.equal oldSplit.Mz [|600.|] "counterexample to the old unsummed reader"
        Expect.equal oldCombined.Mz [|457.|] "old result depends on mobility partitioning"
    testCase "default reader integrates the actual mobility mzlite fixture and preserves its metadata" <| fun _ ->
        use reader = MzIO.Reader.getReader (IO.Path.Combine(AppContext.BaseDirectory,"minimalTIMs.mzlite"))
        MzIO.Reader.openConnection reader
        let headers = readMs2Headers reader (MzIO.Reader.getDefaultRunID reader) ignore
        let raw = reader.ReadSpectrumPeaks(headers.[0].Id).Peaks |> Seq.toArray
        Expect.equal raw.Length 1009 "fixture contains observations from multiple mobility scans"
        let spectrum =
            match readMs2Spectrum reader headers 0 with
            | Ok spectrum -> spectrum
            | Error ex -> raise ex
        Expect.equal spectrum.Mz.Length 492 "default reader integrates with the standard 20 ppm cluster width"
        Expect.equal (Array.distinct spectrum.Mz).Length spectrum.Mz.Length "no duplicated mz coordinates survive"
        Expect.floatClose Accuracy.high (Array.sum spectrum.Intensity) (raw |> Array.sumBy (fun p -> p.Intensity)) "conserve measured intensity"
        Expect.equal spectrum.Id headers.[0].Id "spectrum boundary and identifier preserved"
        Expect.equal spectrum.IonMobility headers.[0].IonMobility "scalar precursor mobility preserved"
    testCase "ppm mobility merging conserves intensity and does not chain across distant peaks" <| fun _ ->
        let input = global.MzIO.Commons.Arrays.ArrayWrapper [|
            global.MzIO.Binary.Peak1D(20., 500.008, 0.73)
            global.MzIO.Binary.Peak1D(10., 500., 0.74)
            global.MzIO.Binary.Peak1D(30., 500.004, 0.75) |]
        let mz,intensity = searchPeakArraysWithTolerance 10. input
        Expect.equal intensity [|40.;20.|] "16 ppm total span cannot form one cluster"
        Expect.floatClose Accuracy.high mz.[0] 500.003 "intensity-weighted centroid"
        Expect.equal mz.[1] 500.008 "next cluster remains separate"
        Expect.equal (Array.sum intensity) 60. "conserve intensity"
    testCase "ppm mobility merging has no fixed bin boundary and validates its tolerance" <| fun _ ->
        let mz,intensity = mergeNearbyPeaks 10. [|500.009;500.011|] [|10.;10.|]
        Expect.equal intensity [|20.|] "close peaks across a 0.01 Da boundary combine"
        Expect.floatClose Accuracy.high mz.[0] 500.01 "centroid"
        for ppm in [-1.;nan;infinity] do
            Expect.throws (fun () -> validateMobilityMergePpm ppm) "invalid tolerance"
    testCase "mobility integration conserves intensity at exact mz before top peak selection" <| fun _ ->
        let input = global.MzIO.Commons.Arrays.ArrayWrapper [|
            global.MzIO.Binary.Peak1D(150., 600., 0.74)
            global.MzIO.Binary.Peak1D(55., 457.260817, 0.743380)
            global.MzIO.Binary.Peak1D(57., 457.260817, 0.736677)
            global.MzIO.Binary.Peak1D(102., 457.260817, 0.729971)
            global.MzIO.Binary.Peak1D(10., 457.260818, 0.729971) |]
        let mz,intensity = searchPeakArraysWithTolerance 0. input
        Expect.equal mz [|457.260817;457.260818;600.|] "retain nearby distinct mz coordinates"
        Expect.equal intensity [|214.;10.;150.|] "sum across the stored mobility window"
        Expect.equal (Array.sum intensity) 374. "conserve total intensity"
        let processed = preprocess {settings with TopNPeaks=1;Deisotope=false} 900. 2 mz intensity
        Expect.equal processed.Mz [|457.260817|] "summation must precede top-N selection"
        Expect.equal processed.Intensity [|100.|] "normalize the integrated base peak"
    testCase "mobility integration ignores invalid observations without poisoning a valid sum" <| fun _ ->
        let input = global.MzIO.Commons.Arrays.ArrayWrapper [|
            global.MzIO.Binary.Peak1D(20., 500., 0.7)
            global.MzIO.Binary.Peak1D(nan, 500., 0.8)
            global.MzIO.Binary.Peak1D(10., nan, 0.8)
            global.MzIO.Binary.Peak1D(-5., 500., 0.8) |]
        Expect.equal (searchPeakArrays input) ([|500.|],[|20.|]) "valid intensities survive for every scoring pass"
    testCase "spectra without per-peak mobility retain the original peak arrays" <| fun _ ->
        let input = global.MzIO.Commons.Arrays.ArrayWrapper [|
            global.MzIO.Binary.Peak1D(10., 600.)
            global.MzIO.Binary.Peak1D(20., 500.)
            global.MzIO.Binary.Peak1D(30., 500.) |]
        Expect.equal (searchPeakArrays input) ([|600.;500.;500.|],[|10.;20.;30.|]) "ordinary spectrum path unchanged"
    testCase "a failed later slice cannot emit partial results or disappear from the failure count" <| fun _ ->
        let attempt={TooFewPeaks=false;Failed=false;Searched=false;Histogram=Array.zeroCreate 64;Targets=[||];Decoys=[||]}
        let header={Id="spectrum";Position=0;ScanTime=1.;PrecursorMz=table.Mass.[0]/2.+protonMass;Charge=2;IonMobility=nan}
        let spectrum : Ms2Spectrum={Id=header.Id;ScanNr=0;ScanTime=1.;PrecursorMz=header.PrecursorMz;Charge=2;IonMobility=nan;Mz=peaks.Mz;Intensity=peaks.Intensity}
        let logs=ResizeArray<string>()
        use stream=new IO.MemoryStream()
        use writer=new IO.StreamWriter(stream)
        let run={Path="test";Log=logs.Add;Reader=Unchecked.defaultof<_>;Headers=[|header|];AttemptStart=[|0;1|];Attempts=[|attempt|];Kept=null;Scratch=Scratch();Writer=writer;Stopwatch=Diagnostics.Stopwatch.StartNew();PhaseTicks=Profile.create()}
        let index1=buildFragmentIndex table 0.01 500. 1 0 10 ignore
        searchSliceWith (fun _ _ _ -> Ok spectrum) settings table [|2|] 10 index1 0 10 0 run
        Expect.isTrue attempt.Searched "first slice searched"
        let index2=buildFragmentIndex table 0.01 500. 1 10 21 ignore
        searchSliceWith (fun _ _ _ -> Error(IO.IOException("injected second-slice read error"))) settings table [|2|] 10 index2 10 21 1 run
        Expect.isTrue attempt.Failed "incomplete histogram must be unusable"
        let classic={CalcIonSeries=(fun _ -> failwith "must not score");ParseSequence=(fun _ _ -> []);ParseRank=(fun _ -> []);AndromedaPMinPMax=(4,10);AndromedaTolerancePpm=20.}
        finishRunWith (fun _ _ _ -> failwith "failed attempt must not read peaks again") settings classic table [|2|] run
        Expect.equal stream.Length 0L "no partial rows"
        Expect.isTrue (logs |> Seq.exists (fun s -> s.Contains("1 failed"))) "count failed attempt"
    testCase "failed input acquisition preserves earlier outputs and releases the opened SQLite reader" <| fun _ ->
        let path name=IO.Path.Combine(AppContext.BaseDirectory,name)
        let p=Json.ReadAndDeserialize<Dto.PeptideSpectrumMatchingTIMsParams> (path "defaultParams.json") |> Dto.PeptideSpectrumMatchingTIMsParams.toDomain
        use cn=SearchDB.getDBConnection (path "MinimalTIMs.db")
        let sdb,t=prepareTable p cn ignore
        let folder=IO.Path.Combine(IO.Path.GetTempPath(),"tims-psm-"+Guid.NewGuid().ToString("N"))
        IO.Directory.CreateDirectory(folder) |> ignore
        let good=IO.Path.Combine(folder,"good.mzlite")
        let bad=IO.Path.Combine(folder,"bad.mzlite")
        IO.File.Copy(path "minimalTIMs.mzlite",good)
        IO.File.Copy(path "minimalTIMs.mzlite",bad)
        try
            do
                use db=new System.Data.SQLite.SQLiteConnection(sprintf "Data Source=%s;Pooling=False" bad)
                db.Open()
                use cmd=new System.Data.SQLite.SQLiteCommand("ALTER TABLE Spectrum RENAME COLUMN Description TO BrokenDescription",db)
                cmd.ExecuteNonQuery() |> ignore
            let output=IO.Path.Combine(folder,"good.psm")
            IO.File.WriteAllText(output,"preserve this previous result")
            Expect.throws (fun () -> scoreRuns p folder ignore 0. 2 false sdb t [|good;bad|]) "header query should fail after reader acquisition"
            Expect.equal (IO.File.ReadAllText output) "preserve this previous result" "all inputs must open before outputs"
            use released=new IO.FileStream(bad,IO.FileMode.Open,IO.FileAccess.ReadWrite,IO.FileShare.None)
            Expect.isTrue released.CanWrite "failed acquisition released the SQLite handle"
        finally
            let resolved=IO.Path.GetFullPath(folder)
            let temp=IO.Path.GetFullPath(IO.Path.GetTempPath()).TrimEnd(IO.Path.DirectorySeparatorChar)
            if IO.Path.GetDirectoryName(resolved) <> temp || not (IO.Path.GetFileName(resolved).StartsWith("tims-psm-")) then
                failwith "Refusing cleanup outside the test temporary directory."
            IO.Directory.Delete(folder,true)
    testCase "index restriction keeps every precursor candidate including isotope windows" <| fun _ ->
        let t={table with Mass=Array.init 21 (fun i -> 300.+float i*0.5)}
        let masses=[|302.;306.|]
        let eligible=reachableRanks settings t masses
        for rank=0 to t.Count-1 do
            let expected=masses |> Array.exists (fun m -> precursorWindows settings t m |> Array.exists (fun w -> rank>=w.Lo && rank<w.Hi))
            Expect.equal eligible.[rank] expected "precursor-window union"
        let restricted=buildFragmentIndexFor eligible t 0.01 500. 1 0 t.Count ignore
        let full=buildFragmentIndex t 0.01 500. 1 0 t.Count ignore
        let entries index =
            [| for b=0 to index.BinCount-1 do
                for i=0 to index.BinLength.[b]-1 do
                    yield b,index.Chunks.[index.BinChunk.[b]].[index.BinLocal.[b]+i] |]
        Expect.equal (entries restricted) (entries full |> Array.filter (fun (_,entry) -> eligible.[entry &&& RankMask])) "same bin, rank, ion type and multiplicity"
        let bounds=sliceBoundariesFor eligible t 4 0.01 500.
        let sliced=[| for i=0 to 3 do yield! entries (buildFragmentIndexFor eligible t 0.01 500. 1 bounds.[i] bounds.[i+1] ignore) |]
        Expect.equal (Array.sort sliced) (entries restricted |> Array.sort) "slicing uses the same restricted population"
    testCase "reused interval buffer equals the original union for random peaks and charges" <| fun _ ->
        let rng=Random 913
        for _ in 1..100 do
            let peaks={Mz=Array.init 150 (fun _ -> 50.+rng.NextDouble()*1500.) |> Array.sort;Intensity=Array.init 150 (fun _ -> rng.NextDouble()*100.)}
            for z in [1;2;3] do
                let expected=massIntervals settings 0.01 z peaks
                let buffer=Array.zeroCreate (150*z)
                let n=massIntervalsInto settings 0.01 z peaks buffer
                Expect.equal (Array.take n buffer) expected "same intervals and maximum intensities"
    testCase "bounded selection and histogram equal exhaustive exact ordering" <| fun _ ->
        let rng=Random 1913
        for limit in [1;3;10;30] do
            let scratch=Scratch()
            scratch.Reset 200 16
            let expected=ResizeArray<Candidate>()
            let refine c={c with Exact=float (c.Rank%7);ExactMatched=c.MatchedB+c.MatchedY;Total=10}
            for rank=0 to 99 do
                for label=0 to 1 do
                    let slot=rank*2+label
                    let nb,ny=rng.Next(6),rng.Next(6)
                    scratch.Slots.[slot].CountB <- nb
                    scratch.Slots.[slot].CountY <- ny
                    scratch.Slots.[slot].SumB <- float nb*10.
                    scratch.Slots.[slot].SumY <- float ny*10.
                    if nb+ny>=1 then
                        expected.Add {Rank=rank;IsDecoy=label=1;MatchedB=nb;MatchedY=ny;Hyperscore=hyperscore nb ny (float nb*10.) (float ny*10.);Exact=nan;ExactMatched=0;Total=0}
            let hist,t,d=collect {settings with ReportedHitsPerLabel=limit} scratch [|{Lo=0;Hi=100;Offset=0}|] refine
            for isDecoy,actual in [false,t;true,d] do
                let reference=expected |> Seq.filter (fun c -> c.IsDecoy=isDecoy && c.MatchedB+c.MatchedY>=4) |> Seq.map refine |> Array.ofSeq
                Array.Sort(reference,candidateOrder)
                Expect.equal actual (Array.truncate limit reference) "full exact order including all ties"
            let model=expectationModel (expected |> Seq.map (fun c -> c.Hyperscore) |> Array.ofSeq)
            let direct=expectationModelOfHistogram hist
            for score in [0.;5.;10.;20.] do Expect.equal (direct score) (model score) "full modelling population"
    testCase "encoded residue parsing and cached masses preserve all fragment families" <| fun _ ->
        let path name=IO.Path.Combine(AppContext.BaseDirectory,name)
        let p=Json.ReadAndDeserialize<Dto.PeptideSpectrumMatchingTIMsParams> (path "defaultParams.json") |> Dto.PeptideSpectrumMatchingTIMsParams.toDomain
        use cn=SearchDB.getDBConnection (path "MinimalTIMs.db")
        let sdb,t=prepareTable p cn ignore
        let classic=toClassicScoring p sdb t
        for rank in 0 .. max 1 (t.Count/500) .. t.Count-1 do
            let parsed=classic.ParseSequence (int(t.GlobalMod rank)) (t.SequenceString rank)
            Expect.equal (classic.ParseRank rank) parsed "labels, modifications and dropped terminators"
            let expected=Fragmentation.Series.fragmentMasses p.nTerminalSeries p.cTerminalSeries sdb.MassFunction parsed
            let actual=classic.CalcIonSeries (classic.ParseRank rank)
            Expect.equal actual expected "identical target and decoy ion families"
    testCase "parameter validation rejects invalid tolerances and ceilings" <| fun _ ->
        let p=Json.ReadAndDeserialize<Dto.PeptideSpectrumMatchingTIMsParams> (IO.Path.Combine(AppContext.BaseDirectory,"defaultParams.json")) |> Dto.PeptideSpectrumMatchingTIMsParams.toDomain
        validateParameters {p with PrecursorTolerancePPM=0.} 0. 1
        for bad in [nan;infinity;-1.] do
            Expect.throws (fun () -> validateParameters {p with FragmentTolerancePPM=bad} 0. 1) "fragment tolerance"
            Expect.throws (fun () -> validateParameters p bad 1) "ceiling"
        Expect.throws (fun () -> validateParameters p 0. 0) "slice count"
        Expect.throws (fun () -> validateParameters {p with FallbackChargeStates=[]} 0. 1) "no fallback charge"
    testCase "full collect rescore merge keeps the late exact match for every slice count" <| fun _ ->
        let scores, target, decoy = pass 0 21
        Expect.isTrue (target |> Array.exists (fun c -> c.Rank=20 && c.ExactMatched=4)) "late exact target must survive"
        for count in [2;3;4;7;21] do
            let boundaries=Array.init (count+1) (fun i -> 21*i/count)
            let parts=Array.init count (fun i -> pass boundaries.[i] boundaries.[i+1])
            let t=parts |> Array.fold (fun acc (_,t,_) -> mergeTop acc t 10) [||]
            let d=parts |> Array.fold (fun acc (_,_,d) -> mergeTop acc d 10) [||]
            Expect.equal t target "target identities and exact scores independent of slices"
            Expect.equal d decoy "decoys independent of slices"
            let combined=Array.zeroCreate scores.Length
            for hist,_,_ in parts do
                for i=0 to hist.Length-1 do combined.[i] <- combined.[i]+hist.[i]
            Expect.equal combined scores "modelling population independent of slices"
            Expect.isLessThanOrEqual t.Length 10 "retained state remains bounded"
    testCase "zero tolerance includes every exact mass and both boundaries are inclusive" <| fun _ ->
        let exact=precursorWindows {settings with PrecursorTolerancePpm=0.; IsotopeErrors=[|0|]} table table.Mass.[0]
        Expect.equal (exact |> Array.sumBy (fun w -> w.Hi-w.Lo)) 21 "exact mass is not an empty interval"
        let small={table with Count=3; Mass=[|1000.-0.02; 1000.;1000.+0.02|]}
        let windows=precursorWindows {settings with IsotopeErrors=[|0|]} small 1000.
        Expect.equal (windows |> Array.sumBy (fun w -> w.Hi-w.Lo)) 3 "include upper and lower boundary"
    testCase "truncated root at a token boundary is rejected" <| fun _ ->
        let json="""{"ID":"x","properties":{"MS:1000511":{"Values":["2"]}}"""
        Expect.isNotNull (SpectrumDescription.tryParse json).Error "missing final brace is an error"
    testCase "empty, scalar, array and multiple-root descriptions are rejected" <| fun _ ->
        for json in ["";"null";"[]";"{}{}";"{\"ID\":\"x\""] do
            Expect.isNotNull (SpectrumDescription.tryParse json).Error json
    testCase "every truncation of a real description is rejected" <| fun _ ->
        let json=IO.File.ReadAllText(IO.Path.Combine(AppContext.BaseDirectory,"sample_description.json")).Trim()
        for length in 0 .. json.Length-1 do
            Expect.isNotNull (SpectrumDescription.tryParse (json.Substring(0,length))).Error (sprintf "prefix length %i" length)
    testCase "a read failure marks all charges of only the affected spectrum" <| fun _ ->
        let attempt () = {TooFewPeaks=false; Failed=false; Searched=true; Histogram=[|0;1;0|];Targets=[||];Decoys=[||]}
        let header id = {Id=id;Position=0;ScanTime=1.;PrecursorMz=500.;Charge=0;IonMobility=nan}
        let run={Path="test";Log=ignore;Reader=Unchecked.defaultof<_>;Headers=[|header "a";header "b"|]
                 AttemptStart=[|0;2;3|];Attempts=Array.init 3 (fun _ -> attempt());Kept=null;Scratch=Scratch();Writer=null
                 Stopwatch=Diagnostics.Stopwatch.StartNew();PhaseTicks=Profile.create()}
        match readMs2Spectrum run.Reader run.Headers 0 with
        | Ok _ -> failtest "read should fail"
        | Error ex -> markReadFailure run 0 ex
        Expect.equal (run.Attempts |> Array.map (fun a -> a.Failed)) [|true;true;false|] "partial search must be unusable"
    testCase "duplicate output stems are rejected before writing" <| fun _ ->
        Expect.throws (fun () -> validateOutputPaths "out" [|"a/sample.mzlite";"b/sample.mzlite"|]) "collision"
        validateOutputPaths "out" [|"a/first.mzlite";"b/second.mzlite"|]
    testCase "merge bounds its result even with an empty side" <| fun _ ->
        let _,t,_=pass 0 21
        Expect.equal (mergeTop [||] t 3).Length 3 "fresh only"
        Expect.equal (mergeTop t [||] 3).Length 3 "kept only"
    testCase "nonfinite peaks cannot poison preprocessing" <| fun _ ->
        let s=preprocess settings 500. 2 [|100.;nan;200.;300.|] [|100.;200.;infinity;50.|]
        Expect.equal s.Mz [|100.;300.|] "only finite valid peaks"
]
