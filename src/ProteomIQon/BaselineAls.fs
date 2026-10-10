namespace ProteomIQon

// Asymmetric least squares baseline of the PSM based quantification tools, only for speed. It is identical to
// FSharp.Stats Signal.Baseline.baselineAls' (same arithmetic, same values) and a candidate for upstreaming to FSharp.Stats,
// after which the tools can call FSharp.Stats again. Do not change its arithmetic.
module BaselineAls =

    /// Asymmetric least squares baseline with the arithmetic of FSharp.Stats.Signal.Baseline.baselineAls'
    /// (Givens QR of FSharp.Stats.Algebra.LinearAlgebraManaged.leastSquares, maxIterations + 1 solves), on the
    /// nonzero entries of the pentadiagonal system only instead of n x n matrices.
    let baselineAls (maxIterations: int) (lambdaExponent: int) (p: float) (y: float []) =
        let n = y.Length
        if n < 3 then FSharp.Stats.Signal.Baseline.baselineAls' maxIterations lambdaExponent p y |> Array.ofSeq else
        let lambda = 10. ** float lambdaExponent
        // diagonals 0, 1 and 2 of lambda D'D, D the second difference matrix (column r has 1, -2, 1 at rows r, r + 1, r + 2)
        let d0 = Array.zeroCreate n
        let d1 = Array.zeroCreate (n - 1)
        let d2 = Array.zeroCreate (n - 2)
        for r in 0 .. n - 3 do
            d0.[r] <- d0.[r] + 1.
            d0.[r + 1] <- d0.[r + 1] + 4.
            d0.[r + 2] <- d0.[r + 2] + 1.
            d1.[r] <- d1.[r] - 2.
            d1.[r + 1] <- d1.[r + 1] - 2.
            d2.[r] <- d2.[r] + 1.
        let d0 = d0 |> Array.map (fun v -> lambda * v)
        let d1 = d1 |> Array.map (fun v -> lambda * v)
        let d2 = d2 |> Array.map (fun v -> lambda * v)
        // Givens rotation of the sparse rows x and y (ApplyGivens of QRSparse), returns the encoded rotation
        let applyGivens (x: ResizeArray<int*float>) (y: ResizeArray<int*float>) =
            let a = snd x.[0]
            let b = snd y.[0]
            let mutable c = -1.0
            let mutable s = -1.0
            if b = 0.0 then
                c <- 1.0
                s <- 0.0
            elif abs b > abs a then
                let tau = -a / b
                s <- 1.0 / sqrt (1.0 + tau * tau)
                c <- s * tau
                if c < 0.0 then
                    c <- -c
                    s <- -s
            else
                let tau = -b / a
                c <- 1.0 / sqrt (1.0 + tau * tau)
                s <- c * tau
            x.[0] <- (fst x.[0], c * a - s * b)
            y.RemoveAt 0
            let mutable p = 1
            let mutable q = 0
            while p < x.Count && q < y.Count do
                if fst x.[p] = fst y.[q] then
                    let xNew = c * snd x.[p] - s * snd y.[q]
                    let yNew = s * snd x.[p] + c * snd y.[q]
                    if xNew <> 0.0 then
                        x.[p] <- (fst x.[p], xNew)
                        p <- p + 1
                    else
                        x.RemoveAt p
                    if yNew <> 0.0 then
                        y.[q] <- (fst y.[q], yNew)
                        q <- q + 1
                    else
                        y.RemoveAt q
                elif fst x.[p] < fst y.[q] then
                    let k = fst x.[p]
                    let xNew = c * snd x.[p]
                    let yNew = s * snd x.[p]
                    x.[p] <- (k, xNew)
                    p <- p + 1
                    y.Insert(q, (k, yNew))
                    q <- q + 1
                else
                    let k = fst y.[q]
                    let xNew = -s * snd y.[q]
                    let yNew = c * snd y.[q]
                    x.Insert(p, (k, xNew))
                    p <- p + 1
                    y.[q] <- (k, yNew)
                    q <- q + 1
            if p < x.Count then
                while p < x.Count do
                    let k = fst x.[p]
                    let xNew = c * snd x.[p]
                    let yNew = s * snd x.[p]
                    x.[p] <- (k, xNew)
                    p <- p + 1
                    y.Insert(q, (k, yNew))
                    q <- q + 1
            elif q < y.Count then
                while q < y.Count do
                    let k = fst y.[q]
                    let xNew = -s * snd y.[q]
                    let yNew = c * snd y.[q]
                    x.Insert(p, (k, xNew))
                    p <- p + 1
                    y.[q] <- (k, yNew)
                    q <- q + 1
            if abs s < c then s
            elif s >= 0.0 then abs (1.0 / c)
            elif 1.0 / c < 0.0 then 1.0 / c
            else (1.0 / c) * -1.0
        let solve (w: float []) =
            // QR of W + lambda D'D, row by row (QRSparse)
            let r = Array.init n (fun _ -> ResizeArray<int*float>())
            let qs = Array.init n (fun _ -> ResizeArray<int*float>())
            for a in 0 .. n - 1 do
                let mutable row = ResizeArray<int*float>()
                for j in max 0 (a - 2) .. min (n - 1) (a + 2) do
                    let v =
                        match j - a with
                        | 0 -> w.[a] + d0.[a]
                        | 1 -> d1.[a]
                        | -1 -> d1.[j]
                        | 2 -> d2.[a]
                        | _ -> d2.[j]
                    if v <> 0.0 then row.Add((j, v))
                let mutable q = 0
                while row.Count > 0 && fst row.[0] < a do
                    let j = fst row.[0]
                    if r.[j].Count <= 0 || fst r.[j].[0] > j then
                        let temp = row
                        row <- ResizeArray<int*float>(r.[a])
                        r.[a] <- ResizeArray<int*float>(temp)
                        qs.[a].Insert(q, (j, 1.0))
                    else
                        let ret = applyGivens r.[j] row
                        qs.[a].Insert(q, (j, ret))
                    q <- q + 1
                r.[a] <- ResizeArray<int*float>(row)
            // Q' (w .* y) (MultiplyVectorWithTransposeGivensMatrix)
            let v = Array.map2 (*) w y
            for i in 0 .. n - 1 do
                for (k, tau) in qs.[i] do
                    if tau = 1.0 then
                        let temp = v.[i]
                        v.[i] <- v.[k]
                        v.[k] <- temp
                    else
                        let mutable c = -1.0
                        let mutable s = -1.0
                        if abs tau < 1.0 then
                            s <- tau
                            c <- sqrt (1.0 - s * s)
                        else
                            c <- 1.0 / tau
                            s <- sqrt (1.0 - c * c)
                            if c < 0.0 then
                                c <- -c
                                s <- -s
                        let newxk = c * v.[k] - s * v.[i]
                        v.[i] <- s * v.[k] + c * v.[i]
                        v.[k] <- newxk
            // back substitution on R (SolveTriangularLinearSystems), columns from the last to the first
            let z = Array.zeroCreate n
            for i in n - 1 .. -1 .. 0 do
                let mutable s = v.[i]
                let mutable diag = 0.0
                for e in r.[i].Count - 1 .. -1 .. 0 do
                    let (j, value) = r.[i].[e]
                    if j > i then s <- s - value * z.[j]
                    elif j = i then diag <- value
                z.[i] <- s / diag
            z
        let mutable w = Array.create n 1.
        for _ in 1 .. maxIterations do
            let z = solve w
            w <- Array.map2 (fun yi zi -> if yi > zi then p else 1. - p) y z
        solve w
