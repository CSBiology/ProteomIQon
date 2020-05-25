namespace ProteomIQon

open BioFSharp
open System
open System.IO
open FSharp.Plotly
open FSharp.Stats
open FSharp.Stats.Signal
open FSharp.Stats.Fitting.NonLinearRegression


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
    let identify (x:float[]) (y:float[]) =
        let trace = (Array.zip x y)
        identifyPeaks p trace
        |> List.choose id 
        |> List.map toIdentifiedPeak
        |> Array.ofList








