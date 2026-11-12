namespace ProteomIQon

open FSharp.Stats
open FSharpAux
open System
open FSharp.Stats.Signal

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
    
        /// Same values as ContinuousWavelet.transform with (-) as getDiff, where that one does not throw; beyond the padded trace the
        /// kernel sees zeros. The constants of the Ricker function are hoisted and the
        /// kernel value of a pair of points is evaluated once, the right side of the left point and the left side of the right point
        /// share it (the same difference of the same two x values, the kernel is even). The sums run in the same order. A kernel value
        /// is a factor of the intensities of its two points only, with both intensities zero (or not an output point) it is not evaluated.
        let private transformFast (kernelValues: float [] ref) (paddedData: (float*float)[]) (borderpadding: int) (wavelet: Wavelet.Ricker) =
            let n = paddedData.Length
            let xs = paddedData |> Array.map fst
            let ys = paddedData |> Array.map snd
            let padd = wavelet.PaddingArea
            let scale = wavelet.Scale
            let ss = pown scale 2
            let twoSs = 2. * ss
            let fakA = 2. / (sqrt(3. * scale) * 1.331335364)
            let ricker (x: float) =
                let xx = x * x
                fakA * (1. - (xx/ss)) * exp (-xx / twoSs)
            let stepSize = xs.[1] - xs.[0]
            let last = n - borderpadding - 1
            // kernel values of all computed rows in one buffer, which the scales of a trace share
            let used = ref 0
            let rowOff : int [] = Array.zeroCreate n
            let rowStart : int [] = Array.zeroCreate n
            let rowLen : int [] = Array.zeroCreate n
            let computeRow i start =
                rowOff.[i] <- used.Value
                rowStart.[i] <- start
                let xi = xs.[i]
                let iOut = i >= borderpadding && i <= last
                let mutable j = start
                let mutable go = true
                while go do
                    // the kernel stops at the end of the padded trace, which is shorter than the kernel for short traces
                    if j >= n then go <- false
                    else
                    let diff = xs.[j] - xi
                    if diff > padd then go <- false
                    else
                        if used.Value = kernelValues.Value.Length then
                            let bigger = Array.zeroCreate (2 * kernelValues.Value.Length)
                            Array.blit kernelValues.Value 0 bigger 0 used.Value
                            kernelValues.Value <- bigger
                        let jOut = j >= borderpadding && j <= last
                        kernelValues.Value.[used.Value] <- if (iOut && ys.[j] <> 0.) || (jOut && ys.[i] <> 0.) then ricker diff else 0.
                        used.Value <- used.Value + 1
                        j <- j + 1
                rowLen.[i] <- j - start
            // the points left of the first output point that have it within the kernel
            let mutable i0 = borderpadding - 1
            while i0 >= 0 && xs.[borderpadding] - xs.[i0] <= padd do
                computeRow i0 borderpadding
                i0 <- i0 - 1
            let center = ricker 0.
            for j in borderpadding .. last do
                computeRow j (j + 1)
            let buf = kernelValues.Value
            Array.init (last - borderpadding + 1) (fun o ->
                let j = borderpadding + o
                let mutable right = 0.
                let off = rowOff.[j]
                for k in 0 .. rowLen.[j] - 1 do
                    right <- right + buf.[off + k] * ys.[j + 1 + k]
                let mutable left = 0.
                let mutable i = j - 1
                let mutable go = true
                while go do
                    if i < 0 then go <- false
                    else
                    let diff = xs.[j] - xs.[i]
                    if diff > padd then go <- false
                    else
                        left <- left + buf.[rowOff.[i] + j - rowStart.[i]] * ys.[i]
                        i <- i - 1
                let correlationValue = right + left + center * ys.[j]
                xs.[j], correlationValue / (Math.Sqrt (Math.Abs scale)) * stepSize
                )

        /// Bit mask (bit j = j-th peak) of the subset of the peaks with the smallest AIC, the subset List.minBy chooses of
        /// FSharpAux.List.powerSetOf: the first of equal values, and the first subset if it has no AIC (nan). The sums of the fits
        /// are added in the order of the peaks, as Vector.sum does, along a depth first walk over the subsets that reuses the sum
        /// of the prefix.
        let private bestSubsetFast (y: float []) (yHats: float [] []) =
            let k = yHats.Length
            let cols = y.Length
            let buffers = Array.init (k + 1) (fun _ -> Array.zeroCreate<float> cols)
            let best = ref 0
            let bestAic = ref nan
            let firstIsNan = ref false
            let rec visit (start: int) (sum: float []) (mask: int) (count: int) =
                for j in start .. k - 1 do
                    let row = yHats.[j]
                    let mask' = mask ||| (1 <<< j)
                    let sum' =
                        if count = 0 then row
                        else
                            let target = buffers.[count]
                            for c in 0 .. cols - 1 do
                                target.[c] <- sum.[c] + row.[c]
                            target
                    let mutable rss = 0.
                    for c in 0 .. cols - 1 do
                        rss <- rss + (y.[c] - sum'.[c]) ** 2.
                    let a = aic cols ((count + 1) * 3) rss
                    if mask' = 1 && Double.IsNaN a then firstIsNan.Value <- true
                    if (a < bestAic.Value) || (a = bestAic.Value && mask' < best.Value) || (best.Value = 0 && not (Double.IsNaN a)) then
                        best.Value <- mask'
                        bestAic.Value <- a
                    visit (j + 1) sum' mask' (count + 1)
            visit 0 null 0 0
            if firstIsNan.Value then 1 else best.Value

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
                    let withHats = groupedPeaks |> List.map (fun (wp,g) -> g, x |> Array.map g.Function)
                    // the groups hold at most 15 peaks (List.chunkBySize 15 in identifyPeaksBy), within the 30 peaks that fit the bit mask of bestSubsetFast
                    let subset = bestSubsetFast y (withHats |> List.map snd |> Array.ofList)
                    withHats |> List.mapi (fun j (g,_) -> j,g) |> List.filter (fun (j,_) -> (subset >>> j) &&& 1 = 1) |> List.map snd
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
            let kernelValues = ref (Array.zeroCreate<float> (1 <<< 12))
            let transformations = 
                scales
                |> Array.map Wavelet.createRicker
                |> Array.map (fun x -> transformFast kernelValues paddedData borderpadding x |> Array.unzip)
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
            identifyPeaks parameters trace
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