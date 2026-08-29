namespace ProteomIQon

open BioFSharp
open BioFSharp.IO
open System.Data
open System.Data.SQLite
open BioFSharp.Mz
open BioFSharp.Mz.SearchDB
open FSharp.Stats
open FSharp.Stats.Fitting.NonLinearRegression
open FSharpAux
open ModificationInfo
open System
open FSharp.Stats.Signal

module SparsePeakArray'= 

    ///
    type SparsePeakArray = { 
        Data            : System.Collections.Generic.IDictionary<int,float>
        MzToBinIdx      : float -> int
        BinIdxToMz      : int -> float 
        }
   
    ///
    let dot (x:SparsePeakArray) (y:SparsePeakArray) =
        x.Data
        |> Seq.fold (fun (acc:float) xi -> 
            let present,yi = y.Data.TryGetValue xi.Key
            if present then acc + (yi * xi.Value) else acc
            ) 0.
           
    ///
    let initMzToBinIdx width offset x = int ((x / width) + offset)

    let initBinIdxToMz width offset x = 
        ((float x) - offset) * width 


    ///
    let peaksToNearestBinVector binWidth offset (minMassBoarder:float) (maxMassBoarder:float) (pkarr:BioFSharp.Mz.PeakArray<_>) =    
        let mzToBinIdx = initMzToBinIdx binWidth offset
        let binIdxToMz = initBinIdxToMz binWidth offset       
        let keyValues = 
            pkarr 
            |> Array.choose (fun p ->  
                if p.Mz < maxMassBoarder && p.Mz > minMassBoarder then 
                    let index = (mzToBinIdx p.Mz) 
                    Some (index, p.Intensity) 
                else 
                    None
                )
            |> Array.groupBy fst 
            |> Array.map (fun (idx,data) -> idx, data |> Array.sumBy snd)
        { 
            Data        = keyValues |> dict
            MzToBinIdx  = mzToBinIdx
            BinIdxToMz  = binIdxToMz
        } 

module DTW' =
                       
    type warpingPaths = {
        Distance: float
        Paths   : Matrix<float>
        }
    
    //window: see   : Only allow for shifts up to this amount away from the two diagonals.
    //maxDist: see  : Stop if the returned distance measure will be larger than this value.
    //maxStep: see  : Do not allow steps larger than this value.
    //maxLengthDiff : Return infinity if difference in length of two series is larger.
    //penalty       : Penalty to add if compression or expansion is applied (on top of the distance).
    //psi           : Relaxation to ignore begin and/or end of sequences (for cylical sequencies)
    ///Dynamic Time Warping.
    ///The full matrix of all warping paths is build.
    ///x: First sequence
    ///y: Second sequence
    ///returns: (DTW distance, DTW matrix)
    let warp window maxDist maxStep maxLengthDiff penalty psi (x:array<float>) (y:array<float>) = 
        let r,c = x.Length, y.Length
        match maxLengthDiff with
        | Some x when abs(r-c) > x -> None 
        | _  -> 
        let window = 
            match window with
            | Some x -> x
            | None   -> max r c             
        let maxStep = 
            match maxStep with
            | Some x -> Some (x**2.)
            | None   -> None 
        let maxDist = 
            match maxDist with
            | Some x -> Some (x**2.)
            | None   -> None
        let penalty = 
            match penalty with
            | Some x -> (x**2.)
            | None   -> 0. 
        let psi = 
            match psi with
            | Some x -> x
            | None   -> 0    
        let dtw = 
            let tmp = Matrix.create (r+1) (c+1) infinity 
            for i = 0 to psi do 
                tmp.[0,i] <- 0.
                tmp.[i,0] <- 0.
            tmp
        ///
        let rec innerLoop upper prevLastUnderMaxDist lastUnderMaxDis i0 i1 i j = 
            if j = upper then 
                i+1, lastUnderMaxDis
            else 
                let d = (x.[i] - y.[j])**2.
                match maxStep with 
                | Some x when d > x -> 
                    innerLoop upper prevLastUnderMaxDist lastUnderMaxDis i0 i1 i (j+1)
                | _ ->                   
                    dtw.[i1,j+1] <- d + (List.min [ dtw.[i0, j];dtw.[i0, j + 1] + penalty;dtw.[i1, j] + penalty]) 
                    match maxDist with 
                    | Some x when dtw.[i1, j + 1] <= x ->
                        innerLoop upper prevLastUnderMaxDist (j) i0 i1 i (j+1)
                    | Some x -> 
                        dtw.[i1, j + 1] <- infinity 
                        if prevLastUnderMaxDist < (j+1) then 
                            (i+1, lastUnderMaxDis)
                        else 
                        innerLoop upper prevLastUnderMaxDist lastUnderMaxDis i0 i1 i (j+1)
                    | None  -> 
                        innerLoop upper prevLastUnderMaxDist lastUnderMaxDis i0 i1 i (j+1)
        ///
        let rec outerLoop lastUnderMaxDist i = 
            if i = r then Some (dtw,i)
            else 
                let prevLastUnderMaxDist = 
                    if lastUnderMaxDist = -1 then
                        infinity |> int |> (*) -1
                    else
                       lastUnderMaxDist
                let last_under_max_dist, i0, i1 = -1, i, i+1
                let lower = 
                    let tmp = i - (max 0 (r - c)) - window + 1
                    max 0 tmp
                let upper = 
                    let tmp = i + (max 0 (c - r)) + window
                    min c tmp
                let nextIdx,prevLastUnderMaxDist = innerLoop upper prevLastUnderMaxDist lastUnderMaxDist i0 i1 i lower
                match maxDist with 
                | Some x when lastUnderMaxDist = -1 ->
                    None 
                | Some x -> 
                    outerLoop lastUnderMaxDist (i+1)
                | None ->
                    outerLoop lastUnderMaxDist (i+1)
        let (dtw') = outerLoop 0 0
        match dtw' with 
        | None -> 
            Some ({Distance=infinity;Paths=dtw}) 
        | Some (dtw',i1) ->
            Matrix.inplace_mapi (fun m n x -> sqrt x) dtw' |> ignore
            let d = 
                //if psi = 0 then 
                     dtw.[i1, (min c (c + window - 1))]  
            //    else 
            //        let ir = i1
            //        let ic = min (c + window - 1)
            //        vr = dtw.get[ir-psi:ir+1, ic]
            //        vc = dtw[ir, ic-psi:ic+1]
            //        mir = np.argmin(vr)
            //        mic = np.argmin(vc)
            //        if vr[mir] < vc[mic]:
            //            dtw[ir-psi+mir+1:ir+1, ic] = -1
            //            d = vr[mir]
            //        else:
            //            dtw[ir, ic - psi + mic + 1:ic+1] = -1
            //            d = vc[mic]
            Some ({Distance=d;Paths=dtw}) 
    
    open FSharpAux
    
    
    /// Returns the best warping path.
    let bestPath (paths:Matrix<float>) =
        let m,n = 
            let m',n' = paths.Dimensions
            m'-1,n'-1
        let pInit = if paths.[m,n] <> -1. then [m-1,n-1] else []
        let rec loop m n path =
            if m <= 0 || n <= 0 then
                match path with 
                | h::t -> t
                | []   -> []
            else 
            let pos = [|paths.[m - 1, n - 1]; paths.[m - 1, n]; paths.[m, n - 1]|]    
            let c = Array.foldi (fun i acc x -> if x < pos.[acc] then i else acc ) 0 pos
            let m',n' = 
                if c = 0 then
                    (m-1),(n-1)
                elif c = 1 then 
                    (m-1),(n)
                else
                    (m),(n-1)
            if paths.[m',n'] <> -1. then 
                loop m' n' ((m'-1,n'-1)::path)
            else loop m' n' (path)
        loop m n pInit        
    
    ///
    let warping_Path window maxDist maxStep maxLengthDiff penalty psi (x:array<float>) (y:array<float>) = 
        let res = warp window maxDist maxStep maxLengthDiff penalty psi x y
        match res with 
        | Some x -> bestPath x.Paths
        | None -> failwith "Warping not possible."
        
    ///
    let distance window maxDist maxStep maxLengthDiff penalty psi (x:array<float>) (y:array<float>) = 
        let res = warp window maxDist maxStep maxLengthDiff penalty psi x y
        match res with 
        | Some x -> x.Distance
        | None -> failwith "Warping not possible."
       
    // Aligns source to target. Given a sourceX it returns the mapping of the xVal closest to
    // sourceX to its targetX counterPart after alignment.
    let align (target:(float*float)[]) (source:(float*float)[]) =
        let xt,yt = target |> Array.unzip
        let xs,ys = source |> Array.unzip
        let optimalPath = warping_Path None None None None None None yt ys
        List.map (fun (ytIdx,ysIdx) -> xt.[ytIdx],ys.[ysIdx]) optimalPath 

    // Aligns source to target. Given a sourceX it returns the mapping of the xVal closest to
    // sourceX to its targetX counterPart after alignment.
    let align' (target:(float*float)[]) (source:(float*float)[]) sourceX =
        let xt,yt = target |> Array.unzip
        let xs,ys = source |> Array.unzip
        let optimalPath = warping_Path None None None None None None yt ys
        List.map (fun (ytIdx,ysIdx) -> xs.[ysIdx],xt.[ytIdx]) optimalPath 
        |> List.minBy (fun (xs,xt) -> abs(xs-sourceX))
    
    let s1 = 
        [|0.;1.;2.;3.;2.;1.;0.;0.;0.;0.;0.;0.;0.;0.;0.;0.;0.;0.;1.;2.;1.;0.|] 
        |> Array.mapi (fun i x -> float i,x)
    
    let s2 = 
        [|0.;0.;0.;0.;0.;1.;2.;2.5;2.;1.;0.;0.;0.;0.;0.;0.;0.;0.;0.;0.;0.1;0.2;0.1;0.;0.;0.5;1.;0.5;0.|] 
        |> Array.mapi (fun i x -> float i + 20.,x)
    
    let s2' = 
        align s1 s2

    
    let zNorm (y:float []) = 
        let m = Seq.mean y
        let s = Seq.stDev y 
        y
        |> Array.map (fun x -> (x-m)/s)
    
    

module SeqIO' = 
    open FSharpAux.IO.SeqIO

    ///Returns a function to format a given value as string
    let inline stringFunction (separator: string) (flatten: bool) (input: 'a) =
        let o = box input
        match o with
        //match string first so that it doesn't get treated as a char array
        | :? string ->
            fun (x: obj) ->
                let sb = new System.Text.StringBuilder()
                sb.Append x |> ignore
                let res = sb.ToString()
                sb.Clear() |> ignore
                res
        | :? System.Collections.IEnumerable ->
            if flatten then
                fun x ->
                    let sb = new System.Text.StringBuilder()
                    let a = x :?> System.Collections.IEnumerable
                    //iterates over Collections.IEnumerable to get entries as objects for the string builder
                    let b = [for i in a do yield box i]
                    b
                    |> Seq.iteri (fun i x ->
                        if i = 0 then
                            sb.AppendFormat("{0}", x) |> ignore
                        else
                            sb.AppendFormat(sprintf "%s{0}" separator, x) |> ignore
                        )
                    let res = sb.ToString()
                    sb.Clear() |> ignore
                    res
            else
                fun x -> 
                    let standardSeparators = [";";"\t";","]
                    let separator = 
                        standardSeparators
                        |> List.find (fun sep -> sep <> separator)
                    let sb = new System.Text.StringBuilder()
                    let a = x :?> System.Collections.IEnumerable
                    //iterates over Collections.IEnumerable to get entries as objects for the string builder
                    let b = [for i in a do yield box i]
                    b
                    |> Seq.iteri (fun i x ->
                        if i = 0 then
                            sb.AppendFormat("{0}", x) |> ignore
                        else
                            sb.AppendFormat(sprintf "%s{0}" separator, x) |> ignore
                        )
                    let res = sb.ToString()
                    sb.Clear() |> ignore
                    res            
        | _ ->
            fun (x: obj) ->
                let sb = new System.Text.StringBuilder()
                sb.Append x |> ignore
                let res = sb.ToString()
                sb.Clear() |> ignore
                res

    let csv (separator: string) header flatten (data: seq<'a>) =
        // Seq.CSVwith derives its schema from the first element and throws on an empty sequence;
        // derive the header from the static record type instead so an empty result still yields a valid file.
        if Seq.isEmpty data then
            if header && Microsoft.FSharp.Reflection.FSharpType.IsRecord typeof<'a> then
                Microsoft.FSharp.Reflection.FSharpType.GetRecordFields typeof<'a>
                |> Array.map (fun p -> p.Name)
                |> String.concat separator
                |> Seq.singleton
            else Seq.empty
        else Seq.CSVwith Seq.valueFunction stringFunction separator header flatten data
    
module FSharpStats' = 

    module Wavelet =

        type WaveletPeak = {
            XAxisIndex  : int 
            ScaleIndex  : int
            Scale       : float
            Correlation : float
            }

        type Gaussian = {
            Amplitude   : float
            XLoc        : float
            Stdev       : float
            Trace       : (float*float) []
            Function    : float -> float
            }

        type PeakGroup = {   
            Start       : float
            End         : float
            Data        : (float*float)[]
            Fits        : Gaussian list
            SumTrace    : float []
            Convolved   : bool
            }

        type Parameters = {
            Borderpadding           : int option
            BorderPadMethod         : Padding.BorderPaddingMethod 
            InternalPaddingMethod   : Padding.InternalPaddingMethod 
            HugeGapPaddingMethod    : Padding.HugeGapPaddingMethod 
            HugeGapPaddingDistance  : float 
            MinPeakDistance         : float option
            MinPeakLength           : float option
            MaxPeakLength           : float 
            NoiseQuantile           : float 
            MinSNR                  : float
            }

        ///
        let weightedMean (weights:seq<'T>) (items:seq<'T>) =
            let sum,n = Seq.fold2 (fun (sum,n) w i -> w*i+sum,n + w ) (0.,0.) weights items
            sum / n

        ///
        let private getColLowBound  (xData: float[]) centerIdx mzTol = 
            let rec loop  (xData: float []) center i mzTol =
                if i <= 0 then 0
                else                                              
                    match (xData.[center] - xData.[i]) with 
                    | x when x >= mzTol -> i+1   
                    | x when x <= mzTol -> loop xData center (i-1) mzTol
            loop xData centerIdx (centerIdx-1) mzTol

        /// 
        let private getColHighBound (xData: float []) centerIdx mzTol = 
            let rec loop (xData: float []) center i mzTol = 
                if i >= xData.Length-1 then xData.Length-1
                else                                             
                    match xData.[i] - xData.[center] with 
                    | x when x >= mzTol -> i-1
                    | x when x <= mzTol -> loop xData center (i+1) mzTol
            loop xData centerIdx (centerIdx+1) mzTol

        ///
        let private getScales (minScale:float) maxScale averageDistance =
            let start : float = Math.Max(minScale,averageDistance*2.)
            let tmp =
                [|start ..  Math.Max(0.05, round 2 (averageDistance/2.)) .. maxScale|]
            tmp

        ///
        let private aic n nParams likelihood = 
            (float n * Math.Log(likelihood / float n)) + 2. * (float nParams) 
    
        ///
        let private optimizeWaveletFits (origData:(float*float) []) (groupedPeaks:(WaveletPeak*Gaussian)list) = 
            match groupedPeaks with 
            | (wp,g)::[] ->
                let min = g.XLoc - 2. * g.Stdev
                let max = g.XLoc + 2. * g.Stdev
                let data = origData |> Array.filter (fun (x,y) -> x >= min  && x <= max )
                let yHat = origData |> Array.map (fst >> g.Function)
                {Start=min;End=max;Data=data;Fits=[g];SumTrace=yHat;Convolved=false}
            | _ -> 
                let min  = groupedPeaks |> List.map (fun (wp,g) -> g.XLoc - 2. * g.Stdev) |> List.min
                let max  = groupedPeaks |> List.map (fun (wp,g) -> g.XLoc + 2. * g.Stdev) |> List.max
                let xy = origData |> Array.filter (fun (x,y) -> x >= min && x <= max) 
                let x,y = xy |> Array.unzip        
                let bestGroup = 
                    groupedPeaks
                    |> List.map (fun (wp,g) -> 
                        let yHat = x |> Array.map g.Function
                        g,yHat
                    )
                    |> FSharpAux.List.powerSetOf
                    |> List.minBy (fun groupedPeaks ->    
                        let m = Matrix.ofJaggedArray (groupedPeaks |> List.map snd |> Array.ofList)
                        let yHat'= m |> Matrix.mapiCols (fun i v -> v |> Vector.sum)
                        let rss = Seq.map2 (fun x y -> (x-y)**2. ) y yHat' |> Seq.sum
                        aic m.NumCols (m.NumRows*3) rss
                    )
                    |> List.map fst
                let sumTrace = 
                    let m = Matrix.ofJaggedArray (bestGroup |> List.map (fun g ->  origData |> Array.map (fst >> g.Function)) |> Array.ofList )
                    let yHat' = m |> Matrix.mapiCols (fun i v -> v |> Vector.sum)
                    yHat'
                    |> Array.ofSeq            
                {Start=min;End=max;Data=xy;Fits=bestGroup;SumTrace=sumTrace;Convolved=true}
             
        ///        
        let identifyPeaksBy (borderpadding:int option) (borderPadMethod:Padding.BorderPaddingMethod) (internalPaddingMethod:Padding.InternalPaddingMethod) (hugeGapPaddingMethod:Padding.HugeGapPaddingMethod) 
            (maxDistance:float) minPeakLength maxPeakLength minPeakDistance noiseQuantile minSNR (trace:(float*float)[]) =
            let minScale = 
                match minPeakLength with 
                | Option.Some v -> v / 6.
                | Option.None -> 0.
            let maxPeakLength = 
                let min = trace |> Array.minBy fst |> fst
                let max = trace |> Array.maxBy fst |> fst
                let delta = max-min
                Math.Min(delta,maxPeakLength) 
            let maxScale = maxPeakLength / 6.
            let averageDistance = Padding.HelperFunctions.getMedianSpacing trace (-)
            let minDistancePadding = averageDistance / 2. 
            let borderpadding =
                match borderpadding with 
                | Some nPoint -> nPoint
                | None        -> Math.Min(((maxPeakLength / minDistancePadding) * 2.)|> int,trace.Length*3)                 
            let paddedData = Padding.pad trace minDistancePadding maxDistance (-) (+) borderpadding borderPadMethod internalPaddingMethod hugeGapPaddingMethod
            let minPeakDistance = 
                match minPeakDistance with 
                | Option.Some v -> Math.Max(v,minDistancePadding*3.)
                | Option.None -> minDistancePadding*3.
            let scales = getScales minScale maxScale minDistancePadding
            let transformations = 
                scales
                |> Array.map Wavelet.createRicker
                |> Array.map (fun x -> ContinuousWavelet.transform paddedData (-) borderpadding x |> Array.unzip)
            let xVals = fst transformations.[0]
            let corrMatrix = 
                transformations
                |> Array.map snd
                |> Matrix.ofJaggedArray
            let noiseLevel = (corrMatrix.Row 0).ToArray() |> FSharp.Stats.Quantile.compute noiseQuantile
            let peakMatrix = 
                corrMatrix
                |> Matrix.mapi (fun m n x -> 
                    if n < 2 || n >= corrMatrix.NumCols-3 then 0. 
                    elif x > corrMatrix.[m,n-1] && x > corrMatrix.[m,n-2] && x > corrMatrix.[m,n+1] && x > corrMatrix.[m,n+2] && x >= (minSNR*noiseLevel) then 
                        x
                    else 0.
                )
            let maxCorrScale = 
                peakMatrix 
                |> Matrix.mapiCols (fun i x ->
                    let mutable maxScaleIdx,maxScale,maxCorr = 0,0.,0.
                    for scale = 1 to x.Length-1 do 
                        if x.[scale] > maxCorr then 
                            maxScaleIdx <- scale 
                            maxScale    <- scales.[scale]
                            maxCorr     <- x.[scale]
                    {XAxisIndex=i ;ScaleIndex=maxScaleIdx; Scale=maxScale; Correlation=maxCorr}   
                )
                |> Array.ofSeq
            let mergedPeaks = 
                let finalPeaks = ResizeArray<WaveletPeak>(100)
                for i = 2 to maxCorrScale.Length-3 do    
                    let currentPeak = maxCorrScale.[i]
                    let currentCorr = currentPeak.Correlation
                    if currentCorr < maxCorrScale.[i-1].Correlation ||
                        currentCorr < maxCorrScale.[i-2].Correlation ||
                        currentCorr < maxCorrScale.[i+1].Correlation ||
                        currentCorr < maxCorrScale.[i+2].Correlation ||
                        currentCorr = 0. then ()
                    else 
                        if currentPeak.ScaleIndex > 0 then 
                            let lowBound  = getColLowBound xVals i (minPeakDistance/2.)
                            let highBound = getColHighBound xVals i (minPeakDistance/2.)
                            let mutable maxCorr = 0.0
                            let mutable maxCol = 0
                            for j = lowBound to highBound do 
                                let scale = maxCorrScale.[j].ScaleIndex
                                if corrMatrix.[scale,j] > maxCorr && scale > 0  then
                                    maxCorr <- corrMatrix.[scale,j]
                                    maxCol <- j 
                            let refinedPeak = maxCorrScale.[maxCol] 
                            if finalPeaks.Count = 0 then 
                                finalPeaks.Add refinedPeak
                            else 
                                let prevPeak = finalPeaks.[finalPeaks.Count-1]
                                let mzDiff = xVals.[refinedPeak.XAxisIndex] - xVals.[prevPeak.XAxisIndex]
                                if mzDiff > minPeakDistance then
                                    finalPeaks.Add refinedPeak
                                elif refinedPeak.Correlation > prevPeak.Correlation then
                                    finalPeaks.[finalPeaks.Count-1] <- refinedPeak
                finalPeaks.ToArray()
            let refinedPeaks = 
                mergedPeaks 
                |> Array.map (fun wp -> 
                    let loc,apex = 
                            let x,y = 
                                paddedData 
                                |> Array.filter (fun (x,y) -> 
                                    let waveletX = xVals.[wp.XAxisIndex]
                                    x > (waveletX-wp.Scale*0.25) && x < (waveletX+wp.Scale*0.25)
                                    )
                                |> Array.unzip
                            let xW = weightedMean y x 
                            let yW  = weightedMean y y  
                            xW,yW
                    let f = FSharp.Stats.Fitting.NonLinearRegression.Table.gaussModel.GetFunctionValue ([|apex;loc;wp.Scale|] |> vector)
                    let peakSim = 
                        trace 
                        |> Array.map fst 
                        |> Array.map (fun x -> x,f x)
                    let gaussian = 
                        {
                            Amplitude   = apex
                            XLoc        = loc
                            Stdev       = wp.Scale
                            Trace       = peakSim
                            Function    = f
                        }
                    wp,gaussian 
                )
                |> Array.sortBy (fun (wp,g) -> g.XLoc)
                |> List.ofArray
                |> List.filter (fun (wp,g) -> g.Stdev > 0. && g.Amplitude > 0.)
            let groupedPeaks = 
                let isOverlapping (currentPeak:(WaveletPeak*Gaussian)) (candidate:(WaveletPeak*Gaussian)) = 
                    (((snd candidate).XLoc < ((snd currentPeak).XLoc + ((snd currentPeak).Stdev * 2.))) && ((snd candidate).XLoc > ((snd currentPeak).XLoc - ((snd currentPeak).Stdev * 2.))) ) ||
                    ( ((snd currentPeak).XLoc < ((snd candidate).XLoc + ((snd candidate).Stdev * 2.))) && ((snd currentPeak).XLoc > ((snd candidate).XLoc - ((snd candidate).Stdev * 2.)))     ) 
                let overlappingPeaks currentFamily peaksLeft = 
                    let rec findOverlappingPeaks  (overlappingPeaks:(WaveletPeak*Gaussian) list) (nonOverlappingPeaks:(WaveletPeak*Gaussian) list) (currentPeak:(WaveletPeak*Gaussian)) (peaksLeft:(WaveletPeak*Gaussian) list) =
                        match peaksLeft with 
                        | [] -> overlappingPeaks,nonOverlappingPeaks
                        | h::t ->
                            if isOverlapping currentPeak h then 
                                findOverlappingPeaks (h::overlappingPeaks) nonOverlappingPeaks currentPeak t 
                            else 
                                findOverlappingPeaks overlappingPeaks (h::nonOverlappingPeaks) currentPeak t 
                              
                    let rec loop newPeakFamily currentFamily peaksLeft =
                        match currentFamily with 
                        | [] -> newPeakFamily,peaksLeft
                        | cp::t -> 
                            let overlappingPeaks,nonOverlappingPeaks = findOverlappingPeaks [] [] cp peaksLeft
                            loop (overlappingPeaks@newPeakFamily) t nonOverlappingPeaks
                    loop [] currentFamily peaksLeft             
                let rec loop (acc:ResizeArray<(WaveletPeak*Gaussian)list>) currentFam peaksLeft =
                    match peaksLeft with 
                    | [] -> 
                        currentFam
                        |> List.chunkBySize 15
                        |> List.iter acc.Add 
                        acc |> Seq.toList
                    | cP::t ->
                        match currentFam with 
                        | [] -> 
                            loop acc (cP::currentFam) (t)                     
                        | cf  ->
                            let overlapping = overlappingPeaks cf (cP::t) 
                            match overlapping with 
                            | [],rest   -> 
                                cf 
                                |> List.chunkBySize 15
                                |> List.iter acc.Add 
                                loop acc [] rest
                            | newPeaks,rest -> 
                                let cf' = newPeaks@cf
                                loop acc cf' rest 
                loop (ResizeArray<(WaveletPeak*Gaussian)list>()) [] refinedPeaks                     
            let origData = 
                xVals 
                |> Array.map (fun x -> 
                    paddedData |> Array.minBy (fun (xx,yy) -> abs(xx-x) )
                )
            let fits = 
                groupedPeaks
                |> List.filter  (List.isEmpty >> not)
                |> List.map (fun x -> x)
                |> List.map (optimizeWaveletFits origData)
                |> List.sortBy (fun x -> x.Start)
            // We could maybe skip the concatenation to keep close peaks together and refine them using 
            // e.g. EM techniques.
            let concatedFits =
                fits
                |> List.map (fun x -> x.Fits)
                |> List.concat
                //|> List.map (fun x-> x.Trace)
            let final =
                let sumTraceMatrix = 
                    concatedFits 
                    |> List.map (fun x -> x.Trace) 
                    |> Array.ofList
                let sumTraceFiltered = 
                    sumTraceMatrix
                    |> JaggedArray.transpose
                    |> Array.map (fun v -> 
                        let max = v |> Array.maxBy snd |> snd
                        v 
                        |> Array.map (fun (x,y) -> if y < max then x,0. else x,y)
                    )
                    |> JaggedArray.transpose
                concatedFits
                |> List.mapi (fun i fit ->
                    let filteredTrace =  
                        sumTraceFiltered.[i]
                        |> Array.map (fun (x,y) -> if abs(x-fit.XLoc) > 2.*fit.Stdev then x, 0. else x,y)
                        |> FSharpAux.Seq.groupWhen (fun (x,y) -> y = 0.)
                        |> Seq.filter (Seq.isEmpty >> not)
                        |> Seq.tryFind (fun signals -> 
                            signals |> Seq.item 0 |> fst < fit.XLoc &&
                            signals |> Seq.last   |> fst > fit.XLoc 
                            )                
                    match filteredTrace with 
                    | Option.Some trace ->
                        let origData' =
                            trace 
                            |> Seq.choose (fun (x',yhat) -> 
                                origData
                                |> Array.tryFind (fun (x,y) -> x=x')
                                )
                        let corr = FSharp.Stats.Correlation.Seq.pearson (Seq.map snd trace) (Seq.map snd origData') 
                        if corr > 0.5 then 
                            Option.Some (origData', fit)
                        else 
                            Option.None
                    |  _ -> Option.None 
                    )
            final        


        ///        
        let identifyPeaks (parameters:Parameters) (trace:(float*float)[]) =
            identifyPeaksBy parameters.Borderpadding parameters.BorderPadMethod parameters.InternalPaddingMethod parameters.HugeGapPaddingMethod
                parameters.HugeGapPaddingDistance parameters.MinPeakLength parameters.MaxPeakLength parameters.MinPeakDistance parameters.NoiseQuantile parameters.MinSNR trace
    
        ///
        let toIdentifiedPeak ((data,peak):(seq<float*float>*Gaussian)) =
                let xData, yData = data|> Array.ofSeq |> Array.unzip
                let apex = Signal.PeakDetection.createPeakFeature -1 peak.XLoc peak.Amplitude
                Signal.PeakDetection.createIdentifiedPeak 
                    apex 
                    None
                    (Signal.PeakDetection.createPeakFeature -1 xData.[0] yData.[0])
                    None
                    (Signal.PeakDetection.createPeakFeature -1 (Array.last xData) (Array.last yData))
                    false
                    false
                    xData
                    yData
    
        ///
        let p = 
            {
                Borderpadding           = None    
                BorderPadMethod         = Padding.BorderPaddingMethod.Random 
                InternalPaddingMethod   = Padding.InternalPaddingMethod.LinearInterpolation 
                HugeGapPaddingMethod    = Padding.HugeGapPaddingMethod.Zero
                HugeGapPaddingDistance  = 100.
                MinPeakDistance         = None
                MinPeakLength           = Some 0.1
                MaxPeakLength           = 1.5 
                NoiseQuantile           = 0.01 
                MinSNR                  = 0.01  
            }

        ///
        let identify parameters (x:float[]) (y:float[]) =
            let trace = (Array.zip x y)
            identifyPeaks p trace
            |> List.choose id 
            |> List.map toIdentifiedPeak
            |> Array.ofList




module SearchDB' =

    module DB' =

        module SQLiteQuery' =

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// Prepares statement to select a ModSequence entry by ModSequenceID
            let prepareSelectModsequenceBySequenceAndGMod (cn:SQLiteConnection) =
                let querystring = "SELECT * FROM ModSequence WHERE Sequence=@sequence AND GlobalMod=@globalMod"
                let cmd = new SQLiteCommand(querystring, cn)
                cmd.Parameters.Add("@sequence", Data.DbType.String) |> ignore
                cmd.Parameters.Add("@globalMod", Data.DbType.Int32) |> ignore
                fun (sequence:string) (globalMod:int) ->
                    cmd.Parameters.["@sequence"].Value  <- sequence
                    cmd.Parameters.["@globalMod"].Value <- globalMod
                    use reader = cmd.ExecuteReader()
                    match reader.Read() with
                        | true ->  (reader.GetInt32(0), reader.GetInt32(1),reader.GetDouble(2), reader.GetInt64(3), reader.GetString(4), reader.GetInt32(5))
                        | false -> -1,-1,nan,-1L,"",-1

    /// Returns a LookUpResult list
    let getThreadSafePeptideLookUpFromFileBySequenceAndGMod (cn:SQLiteConnection) sdbParams = 
        let parseAAString = initOfModAminoAcidString sdbParams.IsotopicMod (sdbParams.FixedMods@sdbParams.VariableMods)
        let selectModsequenceByID = DB'.SQLiteQuery'.prepareSelectModsequenceBySequenceAndGMod cn 
        (fun sequence globalMod -> 
            selectModsequenceByID sequence globalMod 
            |> (createLookUpResultBy parseAAString))   

module ProteinInference' =

    open BioFSharp.PeptideClassification
    open BioFSharp.FileFormats
    open BioFSharp.FileFormats.GFF3
    open BioFSharp.IO.GFF3
    open FSharpAux.IO.SchemaReader.Attribute
    open FSharp.Stats.Interpolation
    open Plotly.NET

    /// For a group of proteins, contains information about all peptides that might be used for its quantification and score calculated for it.
    type InferredProteinClassItemScored =
        {
            GroupOfProteinIDs: string
            PeptideSequence  : string []
            Class            : PeptideEvidenceClass
            TargetScore      : float
            DecoyScore       : float
            Decoy            : bool
            DecoyBigger      : bool
            FoundInDB        : bool
        }

    let createInferredProteinClassItemScored proteinIDs evidenceClass peptideSequences targetScore decoyScore isDecoy decoyBigger foundInDB =
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence   = peptideSequences
            Class             = evidenceClass
            TargetScore       = targetScore
            DecoyScore        = decoyScore
            Decoy             = isDecoy
            DecoyBigger       = decoyBigger
            FoundInDB         = foundInDB
        }

    /// For a group of proteins, contains information about all peptides that might be used for its quantification and score / q-value calculated for it.
    type InferredProteinClassItemQValue =
        {
            GroupOfProteinIDs: string
            PeptideSequence  : string []
            Class            : PeptideEvidenceClass
            TargetScore      : float
            DecoyScore       : float
            QValue           : float
            Decoy            : bool
            DecoyBigger      : bool
            FoundInDB        : bool
        }

    let createInferredProteinClassItemQValue proteinIDs evidenceClass peptideSequences targetScore decoyScore qValue isDecoy decoyBigger foundInDB =
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence   = peptideSequences
            Class             = evidenceClass
            TargetScore       = targetScore
            DecoyScore        = decoyScore
            QValue            = qValue
            Decoy             = isDecoy
            DecoyBigger       = decoyBigger
            FoundInDB         = foundInDB
        }

    /// For a group of proteins, contains information about all peptides that are put into the output file.
    type InferredProteinClassItemOut =
        {
            GroupOfProteinIDs: string
            PeptideSequence  : string
            Class            : PeptideEvidenceClass
            TargetScore      : float
            DecoyScore       : float
            QValue           : float
        }

    let createInferredProteinClassItemOut proteinIDs evidenceClass peptideSequences targetScore decoyScore qValue =
        {
            GroupOfProteinIDs = proteinIDs
            PeptideSequence   = peptideSequences
            Class             = evidenceClass
            TargetScore       = targetScore
            DecoyScore        = decoyScore
            QValue            = qValue
        }

    type PSMInput =
        {
            [<FieldAttribute("PepSequenceID")>]
            PepSequenceID   : int
            [<FieldAttribute("StringSequence")>]
            Seq             :string
            [<FieldAttribute("ModelScore")>]
            ModelScore : float
        }

    ///checks if GFF line describes gene
    let isGene (item: GFFLine<'a>) =
        match item with
        | GFFEntryLine x -> x.Feature = "gene"
        | _ -> false

    ///checks if GFF line describes rna
    let isRNA (item: GFFLine<'a>) =
        match item with
        | GFFEntryLine x -> if x.Feature = "mRNA" then Some x else None
        | _ -> None

    let removeModification pepSeq =
        String.filter (fun c -> System.Char.IsLower c |> not && c <> '[' && c <> ']') pepSeq

    let proteinGroupToString (proteinGroup:string[]) =
        Array.reduce (fun x y ->  x + ";" + y) proteinGroup

    /// Reads geographical information about protein from gff entry and builds the modelinfo of it
    /// This function takes an RNA gff3 entry and therefore will contain the splice variant id of the gene in its result.
    /// This splice variant id should be the same in the given FastA-file.
    let createProteinModelInfoFromEntry i locus (entry:GFFEntry) =
        let attributes = entry.Attributes
        /// Same as in FastA file
        let spliceVariantID =
            match Map.tryFind "ID" attributes with
            | Some res ->
                res.Head
            | None ->
                failwithf "could not find spliceVariantId for locus %s" locus
        let chromosomeID = entry.Seqid

        let direction =
            match entry.Strand with
            |'+' -> PeptideClassification.StrandDirection.Forward
            |'-' -> PeptideClassification.StrandDirection.Reverse
            | _  -> PeptideClassification.StrandDirection.Forward

        PeptideClassification.createProteinModelInfo spliceVariantID chromosomeID direction locus i Seq.empty Seq.empty

    /// By reading GFF creates the protein models (relationships of proteins to each other) which basically means grouping the rnas over the gene loci
    /// TODO: Don't group over order but rather group over id
    let assignTranscriptsToGenes tryParseProteinID (gffLines: seq<GFF3.GFFLine<'a>>)  =
        gffLines
        // transcripts are grouped by the gene they originate from
        |> Seq.groupWhen isGene
        |> Seq.map (fun group ->
            match Seq.head group with
            | GFFEntryLine x ->
                let locus = x.Attributes.["ID"].Head // =genename, this value is used to assign mRNAs of the same gene together
                group
                |> Seq.choose isRNA //the transcripts of the gene are chosen
                |> Seq.mapi (fun i element ->
                    // every transcript of gene gets its own number i and other info is collected from element and used for info of protein
                    let modelInfo = createProteinModelInfoFromEntry i locus element
                    let r =  tryParseProteinID modelInfo.Id
                    // the gff3 id has to be matched with the sequence in the fasta file. therefore the regexpattern is used
                    match r with 
                    | Some v -> 
                        v,
                        modelInfo
                    | None -> 
                        failwithf "could not match gff3 entry id %s with regex pattern. Either gff3 file is corrupt or regexpattern is not correct" modelInfo.Id 
                )

            | _ -> Seq.empty
        )
        |> Seq.concat
        |> Map.ofSeq

    /// Creates a lookup data base to assign peptides to the proteins they are contained in
    let createPeptideProteinRelation (protModels:seq<ProteinModel<'id,'chromosomeId,'geneLocus,'sequence list> option>) =
        let ppRelation = BidirectionalDictionary<'sequence,ProteinModelInfo<'id,'chromosomeId,'geneLocus>>()
        protModels
        |> Seq.iter (fun prot ->
                        // insert peptide-protein relationship
                        // Todo: change type of proteinID in digest
                        match prot with
                        | Some proteinModel ->
                            proteinModel.Sequence
                            |> Seq.iter (fun pepSequence -> ppRelation.Add pepSequence proteinModel.ProteinModelInfo)
                        | None -> ()
        )
        ppRelation

    // Creates a Map of peptides with their highest found score
    let createPeptideScoreMap (psmInputs: PSMInput list list) =
        psmInputs
        |> List.concat
        |> List.groupBy (fun psm -> psm.Seq)
        |> List.map (fun (sequence, psmList) ->
            // This sequence may contain modifications.
            // Depending on the type of lookup this map is used for, the modifications have to be removed.
            sequence,
            psmList
            |> List.maxBy (fun psm -> psm.ModelScore)
            |> fun psm -> psm.ModelScore
        )
        |> Map.ofList

    // Assigns a score to each protein with reverse digested peptides based on the peptides obtained in psm.
    let createReverseProteinScores (reverseProteins: (string*string[])[]) (peptideScoreMap: Map<string,float>) =
        // Remove modifications from map since protein inference was also done using unmodified peptides
        let scoreMapWithoutMods =
            peptideScoreMap
            |> Map.toArray
            |> Array.map (fun (seq, score) ->
                removeModification seq, score
            )
            |> Map.ofArray
        reverseProteins
        |> Array.map (fun (protein, peptides) ->
        protein,
            (
                peptides
                // looks wether the peptides resulting from the reverse digest appear in the peptides from the psm
                |> Array.map (fun pep ->
                    scoreMapWithoutMods.TryFind pep
                    |> (fun x ->
                        match x with
                        | Some score -> score
                        | None -> 0.
                    )
                )
                |> Array.sum,
                peptides
            )
        )
        // only hits are relevant
        |> Array.filter (fun (protein, (score, peptides)) -> score <> 0.)
        |> Map.ofArray

    /// Sums up score of all peptide sequences
    let assignPeptideScores (peptideSequences : string []) (peptideScoreMap : Map<string,float>) =
        peptideSequences
        |> Array.map (fun sequence -> peptideScoreMap.Item sequence)
        |> Array.sum

    /// Looks if the given protein accession is present in a map of identified decoy proteins and assigns its score when found.
    let assignDecoyScoreToTargetScore (proteins: string) (decoyScores: Map<string,(float * string[])>) =
        let prots = proteins |> String.split ';'
        prots
        |> Array.map (fun protein ->
            decoyScores.TryFind protein
            |> fun protOpt ->
                match protOpt with
                // peptides which pointed to the decoy version of this protein are discarded here, they can be included if needed
                | Some (score,peptides) -> score
                | None -> 0.
        )
        |> Array.max

module FDRControl' =

    open Plotly.NET
    open FSharp.Stats.Interpolation
    open FSharpAux.IO

    /// for given data, creates a logistic regression model and returns a mapping function for this model
    let getLogisticRegressionFunction (x:vector) (y:vector) epsilon =
        let alpha =
            match FSharp.Stats.Fitting.LogisticRegression.Univariable.estimateAlpha epsilon x y with
            | Some a -> a
            | None -> failwith "Could not find an alpha for logistic regression of fdr data"
        let weight = FSharp.Stats.Fitting.LogisticRegression.Univariable.coefficient epsilon alpha x y
        FSharp.Stats.Fitting.LogisticRegression.Univariable.predict weight

    ///// Calculates q value mapping funtion for target/decoy dataset
    //let getQValueFunc pi0 bw (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) =
    //    let (scores,_,q) = FDRControl.binningFunction bw pi0 scoreF isDecoyF data
    //    FDRControl.getLogisticRegressionFunction scores q 0.0000001

    ///// Calculates q values for target/decoy dataset
    //let getQValues pi0 (scoreF: 'A -> float) (isDecoyF: 'A -> bool) (data:'A[]) =
    //    let f = getQValueFunc pi0 0.01 scoreF isDecoyF data
    //    Array.map (scoreF >> f) data

    // FDR estimation using MAYU
    // Code form 'stirlingLogFactorial' to 'estimatePi0HG' translated from percolator 'ProteinFDREstimator.cpp'

    let private stirlingLogFacorial (n: float) =
        log(sqrt(2. * Ops.pi  *n)) + n * log(n) - n

    let private exactLogFactorial (n: float) =
        let rec loop i log_fact =
            if i > n then
                log_fact
            else
                let new_log_fact = log_fact + (log i)
                loop (i + 1.) new_log_fact
        loop 2. 0.

    let private logFactorial (n: float) =
        if n < 1000. then
            exactLogFactorial n
        else
            stirlingLogFacorial n

    let private logBinomial (n: float) (k: float) =
        (logFactorial n) - (logFactorial k) - (logFactorial (n - k))

    let private hypergeometric (x: float) (n: float) (w: float) (d: float) =
        //natural logarithm of the probability
        if (d > 0.) then
            exp((logBinomial w x) + (logBinomial (n - w) (d - x)) - (logBinomial n d))
        else 0.

    /// Estimates the false positives given the total number of entries, the number of target hits and the number of decoy hits
    let estimatePi0HG (n: float) (targets: float) (cf: float) =
        let rec loop (fp: float) (logprob: float list) =
            if fp > cf then
                logprob |> List.rev
            else
                let tp = targets - fp
                let w = n - tp
                let prob = hypergeometric fp n w cf
                loop (fp + 1.) (prob::logprob)
        let logprob = loop 0. []
        let sum = logprob |> List.sum
        let logprob_Norm =
            logprob
            |> List.map (fun x ->
                x / sum
            )
        // MAYU rounds here to first decimal
        let expectation_value_FP_PID =
            logprob_Norm
            |> List.foldi (fun i acc x ->
                acc + x * (float i)
            ) 0.
        if (Ops.isNan expectation_value_FP_PID) || (Ops.isInf expectation_value_FP_PID) then
            0.
        else
            expectation_value_FP_PID

    /// Returns proteins sorted into bins according to their size.
    /// proteins are the proteins which were found in either the reverse or forward proteininference, proteinsFromDB are the proteins with peptide sequence
    /// on which the inference was performed.
    let binProteinsLength (proteins: ProteinInference'.InferredProteinClassItemScored []) (proteinsFromDB: (string*string)[]) binCount =
        // need to bin all proteins in the db, not only those with hit?
        let proteinsNotMatchedDB =
            let proteinsMatched =
                proteins
                |> Array.collect (fun protein ->
                    protein.GroupOfProteinIDs
                    |> String.split ';'
                )
                |> Set.ofArray
            // Creates an entry for every protein that is present in the search db and wasn't inferred
            proteinsFromDB
            |> Array.choose (fun (proteinName, peptideSequence) ->
                if Set.contains proteinName proteinsMatched then
                    None
                else
                    Some ((ProteinInference'.createInferredProteinClassItemScored proteinName BioFSharp.PeptideClassification.PeptideEvidenceClass.Unknown [|peptideSequence|] (-1.) (-1.) false false false),
                         float peptideSequence.Length)
            )

        // Adds the length of the peptide sequence to every protein, since they should be binned according to it
        let addedSequenceLength =
            let proteinLengthMap =
                proteinsFromDB
                |> Array.map (fun (protein,sequence) ->
                    protein, sequence.Length
                )
                |> Map.ofArray
            proteins
            |> Array.map (fun ipcis ->
                ipcis,
                let groupOfProteins =
                    ipcis.GroupOfProteinIDs
                    |> String.split ';'
                groupOfProteins
                |> Array.map (fun protein ->
                    match proteinLengthMap.TryFind protein with
                    | None -> failwith "A protein where you are trying to get the length from isn't present in the database"
                    | Some length -> float length
                )
                // lengths are averaged for protein groups
                |> Array.average
            )
        let combined = Array.append addedSequenceLength proteinsNotMatchedDB
        // sorting treats protein groups as one protein with average length. They are also treated as one protein for total and target counts.
        let sortLength =
            combined
            |> Array.sortBy snd
            |> Array.map fst
        let bins =
            let binSize =
                ceil (float sortLength.Length / binCount)
            sortLength
            |> Array.chunkBySize (int binSize)
        bins

    // The original paper of Mayu describes a protein as:
    // FP = if only if all its peptides with q <= threshold are decoy
    // TP = at least one of its peptides with q <= threshold is target
    // However, the way mayu estimates it on the program is like this:
    // FP = any protein that contains a decoy psm with q <= threshold
    // TP = any protein that contains a target psm with q <= threshold
    // They do not consider protein containing both decoy and target psms.
    // ProteomIQon currently uses the picked target decoy approach with the following definitions:
    // FP = protein where the score of decoy hits is higher than score of target hits
    // TP = protein where the score of target hits is higher than score of decoy hits
    // Also, mayu estimates q as the empirical (target-decoy) q value.
    // Percolator estimates q as the empirical (target-decoy) q value and adjusted by pi0
    // Mayu extracts the list of TP and FP proteins from PSM level whereas percolator
    // extract the list of TP and FP proteins from peptide level, this avoids redundancy and
    // gives a better calibration since peptide level q values are re-adjusted in percolator.
    // ProteomIQon extracts TP and FP proteins from the result of the picked target decoy approach.
    // This creates sometimes a difference in the number of TP and FP proteins between percolator, Mayu, and ProteomIQon,
    // which causes a slight difference in the estimated protein FDR.

    /// Calculates the expected false positives for every protein bin and sums them up.
    let expectedFP (proteinBin: ProteinInference'.InferredProteinClassItemScored []) =
        let numberTarget =
            proteinBin
            |> Array.sumBy (fun protein ->
                match not protein.DecoyBigger && protein.FoundInDB with
                | true -> 1.
                | false -> 0.
            )
        let numberDecoy =
            proteinBin
             |> Array.sumBy (fun protein ->
                 match protein.DecoyBigger && protein.FoundInDB with
                 | true -> 1.
                 | false -> 0.
             )
        let total =
            let notFound =
                proteinBin
                |> Array.sumBy (fun protein ->
                    match not protein.FoundInDB with
                    | true -> 1.
                    | false -> 0.
                )
            notFound + numberTarget + numberDecoy
        // MAYU rounds the number of expected false positives for every bin to the first decimal
        let fpInBin = estimatePi0HG total numberTarget numberDecoy
        fpInBin

    /// Calculates the fdr of the data using the MAYU method. The proteinsFromDB is the DB that was used for the inference.
    let calculateFDRwithMAYU (data: ProteinInference'.InferredProteinClassItemScored []) (proteinsFromDB: (string*string)[]) =
        let proteinBins = binProteinsLength data proteinsFromDB 10.
        let estimatedFP =
            proteinBins
            |> Array.fold (fun acc proteinBin -> acc + expectedFP proteinBin) 0.
        let targetCount =
            data
            |> Array.sumBy (fun x ->
            match x.DecoyBigger with
            | true -> 0.
            | false -> 1.
            )
        let fdr =
            if (Ops.isNan estimatedFP) || (Ops.isInf estimatedFP) || estimatedFP = 0. then
                1.
            elif (estimatedFP / targetCount < 0.) || (estimatedFP / targetCount > 1.) then
                1.
            else
                estimatedFP / targetCount
        fdr

    /// Calculates Decoy/Target ratio
    let calculateFDRwithDecoyTargetRatio (data: ProteinInference'.InferredProteinClassItemScored []) =
        // Should decoy Hits be doubled?: Target-decoy search strategy for increasedconfidence in large-scale proteinidentifications by mass spectrometry
        let decoyCount  =
            data
            |> Array.sumBy (fun x ->
                match x.DecoyBigger with
                | true -> 1.
                | false -> 0.
            )
        let targetCount =
            data
            |> Array.sumBy (fun x ->
                match x.DecoyBigger with
                | true -> 0.
                | false -> 1.
            )
        decoyCount/targetCount

    // Assigns a q value to an InferredProteinClassItemScored
    let assignQValueToIPCIS (qValueF: float -> float) (item: ProteinInference'.InferredProteinClassItemScored) =
        if item.Decoy then
            ProteinInference'.createInferredProteinClassItemQValue item.GroupOfProteinIDs item.Class item.PeptideSequence item.TargetScore item.DecoyScore (qValueF item.DecoyScore) item.Decoy item.DecoyBigger true
        else
            ProteinInference'.createInferredProteinClassItemQValue item.GroupOfProteinIDs item.Class item.PeptideSequence item.TargetScore item.DecoyScore (qValueF item.TargetScore) item.Decoy item.DecoyBigger true

    /// Creates a Histogram based on a given score of a target/decoy dataset. Each bin contains the information of the total count, the decoy count and the median score.
    /// (Bin, Count, DecoyCount, Median Score)
    let createTargetDecoyHis bandwidth (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) (data: 'a[]) =
        let halfBw = bandwidth / 2.0
        let scoreDecoyInfo =
            data
            |> Array.map (fun x -> 
                if isDecoy x then
                    {|Score = decoyScoreF x; Decoy = true|}
                else
                    {|Score = targetScoreF x; Decoy = false|}
            )
        scoreDecoyInfo
        |> Array.groupBy (fun x ->
            floor (x.Score / bandwidth)) 
        |> Array.map (fun (k,values) -> 
            let count = (Array.length(values))
            let decoyCount = (values |> Array.filter (fun x -> x.Decoy = true) |> Array.length)
            let medianScore = values |> Array.map (fun x -> x.Score) |> Array.median
            // first part of the tuple only needed for debugging
            if k < 0. then
                ((k  * bandwidth) + halfBw, count, decoyCount, medianScore)
            else
                ((k + 1.) * bandwidth - halfBw, count, decoyCount, medianScore)
        )

    /// Calculates the PEP value based on the ratio of Decoys to targets at a given score
    let calculatePEPValues (totalCountF: 'a -> float) (decoyCountF: 'a -> float) (scoreF: 'a -> float) (dataFreq: 'a[]) = 
        dataFreq
        |> Array.map (fun x -> 
            scoreF x,(decoyCountF x)/(totalCountF x)
        )
        |> Array.sortBy fst
        |> Array.toList

    /// Logit transforms pep values (log10)
    let logitTransformPepValues score pepVal  =
        Array.zip score pepVal
        // 0 and 1 are + and - infinity
        |> Array.filter (fun (y,x) -> x <> 0. && x <> 1.)
        |> Array.map (fun (score,pep) ->
            score,
            log10 (pep/(1.-pep))
        )
        |> Array.unzip

    /// Calculates monotonized PEP values for a target/decoy dataset based on the decoy/target ratio. Entries are binned with a given bandwidth as intital estiamtor based on the scores. 
    /// Returns a function which maps from score to PEP value based on a fit of a linear function using linear regression. The linear regression is performed on the logit transformed 
    /// pep values. The fit focuses on the pep values centered aound the middle of the score distribution
    let initCalculateLin (logger: NLog.Logger) bandwidth (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) (data: 'a[]) =
        let lowerScore, upperScore =
            let decoy = 
                data
                |> Array.filter isDecoy
                |> Array.map decoyScoreF
                |> Array.filter (fun x -> x < 0.)
                |> Array.median
            let target =
                data
                |> Array.filter (isDecoy >> not)
                |> Array.map targetScoreF
                |> Array.filter (fun x -> x > 0.)
                |> Array.median
            decoy, target
        logger.Trace (sprintf "Lower Score: %f; Upper Score: %f" lowerScore upperScore)
        let filteredData =
            data
            |> Array.filter (fun entry ->
                if isDecoy entry then
                    let score = decoyScoreF entry
                    score >= lowerScore && score <= upperScore
                else
                    let score = targetScoreF entry
                    score >= lowerScore && score <= upperScore
            )
        logger.Trace(sprintf "Initial Bandwidth: %f" bandwidth)
        let fittingFunction, score, pep =
            let xPointRange =
                let min = Math.Min((Array.minBy targetScoreF filteredData) |> targetScoreF, (Array.minBy decoyScoreF filteredData) |> decoyScoreF)
                let max = Math.Max((Array.maxBy targetScoreF filteredData) |> targetScoreF, (Array.maxBy decoyScoreF filteredData) |> decoyScoreF)
                max-min
            let upperBW = Math.Min(10., xPointRange/10.)
            [|bandwidth .. 0.1 .. upperBW|]
            |> Array.choose (fun bw ->
                let targetDecoyHis = createTargetDecoyHis bw (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) (filteredData: 'a[])
                let score',pep' = 
                    calculatePEPValues (fun (_,count,_,_) -> float count) (fun (_,_,decoyCount,_) -> float decoyCount) (fun (_,_,_,medianScore) -> medianScore) targetDecoyHis
                    |> Array.ofList
                    |> Array.unzip
                let logitScore, logitPEPVal = logitTransformPepValues score' pep'
                let coeff = Fitting.LinearRegression.OLS.Linear.Univariable.fit (vector logitScore) (vector logitPEPVal)
                let fittingFunction' = (Fitting.LinearRegression.OLS.Linear.Univariable.predict coeff) >> (fun x -> 10.**(x)/(1.+10.**(x)))
                let sos = FSharp.Stats.Fitting.GoodnessOfFit.calculateSumOfSquares fittingFunction' score' pep'
                if coeff.[1] < 0. then
                    Some (sos.Error/sos.Count, fittingFunction', score', pep', bw)
                else
                    None
            )
            |> Array.minBy (fun (error,_,_,_,_) -> error)
            |> fun (error, fit,s,p,bw) ->logger.Trace(sprintf "Chosen Bandwidth: %f" bw); fit,s,p
        fittingFunction



module Drafo = 
    
    module Core =
        
        open Deedle
        open DynamicObj
        open System.Collections.Generic

        type Key() = 
            inherit DynamicObj ()
           
            member this.addCol(columns:(string*'a)) =
                this.SetProperty columns
                this

            override this.ToString() = 
                let sb = new System.Text.StringBuilder()
                // sb.Append
                (this.GetProperties true)        
                |> Seq.iter (fun x -> 
                    let value = x.Value.ToString() 
                    sb.AppendLine(x.Key+": "+value)
                    |> ignore
                    )
                (sb.ToString())

            override this.Equals(b) =
                match b with
                | :? Key as p -> 
                    let propA = (this.GetProperties true) |> Seq.map (fun x -> x.Key,x.Value) 
                    let propB = (p.GetProperties true) |> Seq.map (fun x -> x.Key,x.Value) 
                    Seq.map2 (fun (x1,x2) (y1,y2) -> x1 = unbox y1 && x2 = unbox y2) propA propB 
                    |> Seq.contains false
                    |> not
                | _ -> false

            override this.GetHashCode() = 
                let sb = new System.Text.StringBuilder()
                // sb.Append
                (this.GetProperties true)        
                |> Seq.iter (fun x -> 
                    let value = x.Value.ToString() 
                    sb.Append(value)
                    |> ignore
                    )
                (sb.ToString())
                |> hash   

        ///
        let indexWithColumnValues (keyCols:(seq<string>) ) (f:Frame<_,string>) :Frame<Key,_>= 
            f
            |> Frame.indexRowsUsing (fun s -> 
                    keyCols
                    |> Seq.fold (fun key x -> 
                        let value = s.GetAs<string>(x) 
                        key.addCol (x,value)
                        ) (Key())
                )  

        ///
        let readFrame (fp:string) = Frame.ReadCsv(fp,hasHeaders=true,inferTypes=false,separators="\t")
          
        ///
        let readAndIndexFrame keyCols fp = 
            readFrame fp
            |> indexWithColumnValues keyCols

        ///
        let inline getColumn<'a> column (f:Frame<Key, string> )  =
            f.GetColumn<'a>(column)

        ///
        let inline seriesToFrame (s: Series<Key, 'a>) =
            s
            |> Series.map (fun k s -> 
                (k.GetProperties true) 
                |> Seq.map (fun x -> x.Key,x.Value)
                |> Series.ofObservations
            )    
            |> Frame.ofRows
            |> Frame.addCol "Value" s
            |> Frame.indexRowsOrdinally

        ///
        let inline rowKeyToColumns (f: Frame<Key, string>) =
            let rowKeysAsColumns = 
                f
                |> Frame.mapRows (fun k s -> 
                    (k.GetProperties true) 
                    |> Seq.map (fun x -> x.Key,x.Value)
                    |> Series.ofObservations
                )    
                |> Frame.ofRows
            Frame.join JoinKind.Inner rowKeysAsColumns f 
            |> Frame.indexRowsOrdinally

        ///
        let createFilter (op : 'a -> bool) (s: Series<'KeyType, 'a>) = 
            s
            |> Series.mapValues op

        ///
        let transform (op : 'a -> 'b) (s: Series<'KeyType, 'a>) = 
            s
            |> Series.mapValues op

        ///
        let zip (op : 'a -> 'a -> 'b) (s1: Series<'KeyType, 'a>) (s2: Series<'KeyType, 'a>) = 
            Series.zipInner s1 s2
            |> Series.mapValues (fun (s1,s2) -> 
                op s1 s2
                )

        ///
        let dropAllKeyColumnsBut (keyColumns:seq<string>) (key:Key) = 
            let newK = Key()
            key.GetProperties true
            |> Seq.filter (fun x -> keyColumns |> Seq.contains x.Key )
            |> Seq.fold (fun (k:Key) x -> k.addCol (x.Key,x.Value) ) newK

        ///
        let dropKeyColumns (keyColumns:seq<string>) (key:Key) = 
            let newK = Key()
            key.GetProperties true
            |> Seq.filter (fun x -> keyColumns |> Seq.contains x.Key |> not)
            |> Seq.fold (fun (k:Key) x -> k.addCol (x.Key,x.Value) ) newK

        ///
        let groupTransform (op :'a [] -> 'a -> 'b) modifyKeyColumns (keyColumns:seq<string>) (s: Series<'KeyType, 'a>) =
            s
            |> Series.groupBy (fun k v -> modifyKeyColumns keyColumns k )
            |> Series.mapValues (fun valueCol -> 
                let fInit = valueCol |> Series.values |> Array.ofSeq |> op
                let filterCol = 
                    valueCol
                    |> Series.mapValues fInit
                filterCol 
            )
            |> Series.values
            |> Series.mergeAll

        ///
        let createGroupFilter (op :'a [] -> 'a -> bool) modifyKeyColumns (keyColumns:seq<string>) (s: Series<'KeyType, 'a>) =
            groupTransform op modifyKeyColumns keyColumns s

        ///
        let aggregate (op : seq<'A> -> 'C) modifyKeyColumns (keyColumns:seq<string>) (filters:seq<Series<Key,bool>>) (col:Series<Key,'A>)  :Series<Key,'C> = 
            let filtered = 
                filters
                |> Seq.map (fun s -> 
                    System.Guid.NewGuid(),
                    s |> Series.filterValues id)
                |> Series.ofObservations
                |> Frame.ofColumns
                |> Frame.dropSparseRows
            let colID = (System.Guid.NewGuid())
            filtered
            |> Frame.addCol colID col
            |> Frame.dropSparseRows
            |> Frame.getCol colID
            |> Series.applyLevel (modifyKeyColumns keyColumns) (Series.values >> op)

        ///
        let assemble (cols:seq<(string * #ISeries<Key>)>) =
            Frame.ofColumns cols

        ///
        let pivot (pivotCol:string) (assembeledFrame:Frame<Key,string>) =
            assembeledFrame
            |> Frame.nestBy (fun k -> 
                let value: string option = k.TryGetTypedPropertyValue pivotCol
                value.Value
                )
            |> Series.map (fun k f ->
                    f
                    |> Frame.mapColKeys (fun ck -> sprintf "%s.%s" k ck)
                    |> Frame.mapRowKeys (dropKeyColumns [pivotCol])
                )
            |> Series.values
            |> Frame.mergeAll