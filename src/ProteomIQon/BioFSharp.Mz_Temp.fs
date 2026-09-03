namespace ProteomIQon

open FSharp.Stats
open FSharpAux
open System
open FSharp.Stats.Signal

// TODO upstream to FSharp.Stats (generic dynamic time warping; not BioFSharp.Mz domain)
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
    
    

// TODO upstream to FSharpAux.IO (generic CSV writer incl. empty-result header fix)
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
    
// TODO upstream to FSharp.Stats (wavelet peak detection fork; unify with BioFSharp.Mz FSharp.Stats.fs fork there)
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




// TODO upstream to Drafo (data-frame tooling; not BioFSharp.Mz domain)
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