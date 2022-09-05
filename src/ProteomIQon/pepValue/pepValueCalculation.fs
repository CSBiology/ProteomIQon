(*******************************************************************************
This is derivative work based on Percolator (https://github.com/percolator/percolator)

Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*******************************************************************************)

namespace ProteomIQon

open FSharp.Stats
open FSharpAux
open System
open FDRControl'

module PepValueCalculation =

    let initCalculatePEPValueIRLS (logger: NLog.Logger) bandwidth (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) (data: 'a[]) =
        
        let targetDecoyHis = createTargetDecoyHis bandwidth (isDecoy: 'a -> bool) (decoyScoreF: 'a -> float) (targetScoreF: 'a -> float) (data: 'a[])

        let binSize, negativeCounts, scores =
            targetDecoyHis
            |> Array.map (fun (_,count,decoyCount,medianScore) ->
                float count,float decoyCount,medianScore
            )
            |> Array.sortBy (fun (_,_,score) -> score)
            |> Array.unzip3

        let mutable x = Vector.zeroCreate 1
        let mutable y = Vector.zeroCreate 1
        let mutable m = Vector.zeroCreate 1
        let mutable transf: float -> float = fun x -> x

        let transform deltaLow deltaHigh doLogit doLog (xx: float) =
            if not doLogit && not doLog then
                xx
            elif deltaLow > 0. || deltaHigh > 0. then
                if doLogit then
                    log ((xx * (1. - deltaHigh - deltaLow)) + deltaLow)
                else
                    log (xx + deltaLow)
            elif doLogit then
                log (xx / (1. - xx))
            else
                failwith "unexpected transform input"

        let setData (xx: Vector<float>) =
            x <- Vector.zeroCreate xx.Length
            let minV = xx |> Vector.toArray |> Array.min
            let maxV = xx |> Vector.toArray |> Array.max
            if minV >= 0. && maxV <= 1. then
                transf <-
                    transform
                        (
                            if minV > 0. then 0. else 1e-20
                        )
                        (
                            if maxV < 1. then 0. else 1e-10
                        )
                        true
                        false
            if minV >= 0. then
                transf <-
                    transform
                        (
                            if minV > 0. then 0. else 1e-20
                        )
                        0.
                        false
                        true
            x <- xx |> Vector.map transf

        let lrSetData (xx: Vector<float>) (yy: Vector<float>) (mm: Vector<float>) =
            y <- yy
            m <- mm
            setData xx

        lrSetData 
            (scores |> Vector.ofArray)
            (negativeCounts |> Vector.ofArray)
            (binSize |> Vector.ofArray)

        logger.Trace (
            scores
            |> Array.map (fun x -> string x)
            |> String.concat ";"
        )

        logger.Trace (
            negativeCounts
            |> Array.map (fun x -> string x)
            |> String.concat ";"
        )

        logger.Trace (
            binSize
            |> Array.map (fun x -> string x)
            |> String.concat ";"
        )

        let convergeEpsilon = 1e-4
        let stepEpsilon = 1e-8
        let weightSlope = 1e1
        let scaleAlpha = 1.
        let tao = 2. / (1. + sqrt(5.0))

        let mutable g = Vector.zeroCreate x.Length
        let mutable gNew = Vector.zeroCreate x.Length
        // Change back to zerocreate
        let mutable w = Vector.init x.Length (fun x -> 1.)
        let mutable z = Vector.init x.Length (fun x -> 0.5)
        let mutable gamma = Vector.zeroCreate (x.Length - 2)

        let mutable Q: Matrix<float> = Matrix.zero 1 1
        let mutable Qt: Matrix<float> = Matrix.zero 1 1
        let mutable R: Matrix<float> = Matrix.zero 1 1
        let mutable dx: Vector<float> = Vector.zeroCreate 1

        let mutable p = Array.zeroCreate 1

        let gRange = 35.

        let initg () =
            let n = x.Length
            g <- Vector.zeroCreate n
            gNew <- Vector.zeroCreate n
            w <- Vector.zeroCreate n
            z <- Vector.init n (fun x -> 0.5)
            gamma <- Vector.zeroCreate (n - 2)

        let invlogit (input: float) =
            let e = exp input
            e / (1. + e)

        let logit p =
            log(p / (1. - p))

        let limitg () =
            for ix = (gNew.Length - 1) downto 0 do
                gNew.[ix] <- Math.Min(gRange, Math.Max(-gRange,gNew.[ix]))

        let limitgamma () =
            for ix = (gamma.Length - 1) downto 0 do
                gamma.[ix] <- Math.Min(gRange, Math.Max(-gRange,gamma.[ix]))

        let calcPZW () =
            for ix = (z.Length - 1) downto 0 do
                let e = exp(g.[ix])
                let epsilon = 1e-15
                p.[ix] <- Math.Min(Math.Max(e / (1. + e),epsilon), 1. - epsilon)
                w.[ix] <- Math.Max(m.[ix] * p.[ix] * (1. - p.[ix]), epsilon)
                z.[ix] <- Math.Min(gRange, Math.Max(-gRange, g.[ix] + (y.[ix] - p.[ix] * m.[ix]) / w.[ix]))

        let lrInitg () =
            initg()
            let n = x.Length
            Array.Resize(&p,n)
            gNew <- Vector.zeroCreate n
            for ix = (g.Length - 1) downto 0 do
            let p = (y.[ix] + 0.05) / (m.[ix] + 0.1)
            gNew.[ix] <- log(p / (1. - p))

        let initiateQR () =
            let n = x.Length
            let dx' = Vector.zeroCreate (n - 1)
            for ix = 0 to (n - 2)do
                dx'.[ix] <- x.[ix + 1] - x.[ix]
                if not (dx'.[ix] > 0.) then
                    failwith "value muste be > 0"
            let Q' = Matrix.zero n (n - 2)
            let R' = Matrix.zero (n - 2) (n - 2)
            // Fill Q
            Q'.[0,0] <- 1. / dx'.[0]
            Q'.[1,0] <- -1. / dx'.[0] - 1. / dx'.[1]
            Q'.[1,1] <- 1. / dx'.[1]
            for j = 2 to (n - 3) do
                Q'.[j,j - 2] <- 1. / dx'.[j - 1]
                Q'.[j,j - 1] <- -1. / dx'.[j - 1] - 1. / dx'.[j]
                Q'.[j,j] <- 1. / dx'.[j]
            Q'.[n - 2,n - 4] <- 1. / dx'.[n - 3]
            Q'.[n - 2,n - 3] <- -1. / dx'.[n - 3] - 1. / dx'.[n - 2]
            Q'.[n - 1,n - 3] <- 1. / dx'.[n - 2]
            // Fill R
            for i = 0 to n - 4 do
                R'.[i,i] <- (dx'.[i] + dx'.[i + 1]) / 3.
                R'.[i,i + 1] <- dx'.[i + 1] / 6.
                R'.[i + 1,i] <- dx'.[i + 1] / 6.
            R'.[n - 3, n - 3] <- (dx'.[n - 3] + dx'.[n - 2]) / 3.
            let Qt' = Q' |> Matrix.transpose
            Q <- Q'
            Qt <- Qt'
            R <- R'
            dx <- dx'
                
        let solveInPlace (mat: Matrix<float>) (res: Vector<float>) =
            Algebra.LinearAlgebra.SolveLinearSystem mat res

        
        let splineEval (xx: float) =
            let xxLogit = transf xx
            let n = x.Length
            let right =
                x
                |> Vector.toArray
                |> Array.tryFindIndex (fun e -> e >= xxLogit)
            if right.IsNone then
                let derl = (g.[n - 1] - g.[n - 2]) / (x.[n - 1] - x.[n - 2]) + (x.[n - 1] - x.[n - 2]) / 6. * gamma.[n - 3]
                let gx = g.[n - 1] + (xx - x.[n - 1]) * derl
                gx
            elif x.[right.Value] = xx then
                g.[right.Value]
            elif right.Value > 0 then
                let left = right.Value - 1
                let dr = x.[right.Value] - xx
                let dl = xx - x.[left]
                let gamr = 
                    if right.Value < (n - 1) then
                        gamma.[right.Value - 1]
                    else
                        0.
                let gaml = 
                    if right.Value > 1 then
                        gamma.[right.Value - 1 - 1]
                    else
                        0.
                let h = x.[right.Value] - x.[left]
                let gx = (dl * g.[right.Value] + dr * g.[right.Value - 1]) / h - dl * dr / 6. * ((1.0 + dl / h) * gamr + (1.0 + dr / h) * gaml)
                gx
            else
                let derr = (g.[1] - g.[0]) / (x.[1] - x.[0]) - (x.[1] - x.[0]) / 6. * gamma.[0]
                let gx = g.[0] - (x.[0] - xx) * derr
                gx

        let iterativeReweightedLeastSquares (alpha: float) =
            let mutable step = 0.
            let mutable iter = 0
            // do .. while
            let mutable init = true
            let n = x.Length
            while init || ((step > stepEpsilon || step < 0.) && iter < 20) do
                init <- false
                iter <- iter + 1
                g <- gNew
                calcPZW()
                // strange vector division again
                let diag = ((Matrix.diag (Vector.map2 (fun x y -> x / y) (Vector.init (n) (fun x -> 1.)) w)) * alpha)
                let aWiQ = diag * Q
                let M = R + (Qt * aWiQ)
                gamma <- Qt * z
                gamma <- solveInPlace M gamma
                gNew <- z - (aWiQ*gamma)
                limitg()
                let difference = g - gNew
                step <- (Vector.norm difference) / float n
            g <- gNew

        let evaluateSlope (alpha: float) =
            iterativeReweightedLeastSquares(alpha)
            let n = g.Length
            let mutable mixg = 1
            let mutable maxg = g.[mixg]
            for ix = mixg to (n - 2) do
                //assert(ix=g.index(ix)); //This should be a filled vector
                if g.[ix] >= maxg then
                    maxg <- g.[ix]
                    mixg <- ix
            let mutable maxSlope = -10e6
            let mutable slopeix = -1
            for ix = (mixg + 1) to (n - 3) do
                let slope = g.[ix - 1] - g.[ix]
                if slope > maxSlope then
                    maxSlope <- slope
                    slopeix <- ix
            maxSlope * weightSlope + alpha

        let rec alphaLinearSearchBA (min_p': float) (max_p': float) (p1': float) (p2': float) (cv1': float) (cv2': float) =
            // Minimize Slope score
            // Use neg log of 0<p<1 so that we allow for searches 0<alpha<inf
            let mutable oldCV = 0.
            let mutable min_p = min_p'
            let mutable max_p = max_p'
            let mutable p1 = p1'
            let mutable p2 = p2'
            let mutable cv1 = cv1'
            let mutable cv2 = cv2'
            if cv2 < cv1 then
                // keep point 2
                min_p <- p1
                p1 <- p2
                p2 <- min_p + tao * (max_p - min_p)
                oldCV <- cv1
                cv1 <- cv2
                cv2 <- evaluateSlope(-scaleAlpha*log(p2))
            else
                // keep point 1
                max_p <- p2
                p2 <- p1
                p1 <- min_p + (1. - tao) * (max_p - min_p)
                oldCV <- cv2
                cv2 <- cv1
                cv1 <- evaluateSlope(-scaleAlpha*log(p1))
            if ((oldCV - (min cv1 cv2)) / oldCV < 1e-5 || (abs(p2 - p1) < 1e-10)) then
                if cv1 < cv2 then
                    -scaleAlpha*log(p1)
                else
                    -scaleAlpha*log(p2)
            else
                alphaLinearSearchBA min_p max_p p1 p2 cv1 cv2

        let roughnessPenaltyIRLS () =
            initiateQR()
            initg()
            let mutable p1 = 1. - tao
            let mutable p2 = tao
            let alpha =
                let mutable min_p = 0.
                let mutable max_p = 1.
                let mutable cv1 = (evaluateSlope(-scaleAlpha*log(p1)))
                let mutable cv2 = (evaluateSlope(-scaleAlpha*log(p2)))
                alphaLinearSearchBA 
                    min_p
                    max_p
                    p1
                    p2
                    cv1
                    cv2
            iterativeReweightedLeastSquares(alpha)

        let revLogitAndPi01 (peps: float[])=
            let top = Math.Min(1., Math.Exp(peps |> Array.max))
            let mutable crap = false
            peps
            |> Array.map (fun x ->
                if crap then
                    top
                else
                    let temp = Math.Exp x
                    if temp >= top then
                        crap <- true
                        top
                    else
                        temp
            )

        lrInitg()
        limitg()
        limitgamma()
        roughnessPenaltyIRLS()
        data
        |> Array.filter (isDecoy >> not)
        |> Array.sortByDescending targetScoreF
        |> fun filteredSorted ->
            filteredSorted
            |> Array.map (fun x -> 
                x
                |> targetScoreF,
                x
                |> targetScoreF
                |> splineEval
            )
        |> Array.unzip
        |> fun (x, y) ->
            x,
            y
            |> revLogitAndPi01
            |> fun arr ->
                let head::tail = arr |> Array.rev |> Array.toList
                tail
                |> List.fold (fun (acc: float list) newPEPValue ->
                    let pepValue = acc.Head
                    if newPEPValue > pepValue then
                        pepValue::acc
                    else
                        newPEPValue::acc
                )[head]
                |> Array.ofList
        |> fun (x, y) ->
            let coeff = FSharp.Stats.Interpolation.LinearSpline.initInterpolate x y
            let fitLinSp = Interpolation.LinearSpline.interpolate coeff
            fitLinSp