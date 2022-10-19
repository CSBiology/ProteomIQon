namespace ProteomIQon


open System.IO
open System
open FSharpAux.Colors.Table.StatisticalGraphics24
open FSharp.Stats
open FSharpAux.IO.SchemaReader
open Plotly.NET
open BioFSharp
open Microsoft
// open Microsoft.ML
// open Microsoft.ML.Data
// open Microsoft.ML.AutoML   
open Dto
open Dto.QuantificationResult
open Core.InputPaths


// The Spline module was implemented according to a julia implementation available at: https://github.com/nignatiadis/SmoothingSplines.jl/blob/master/LICENSE.md 
// Copyright (c) 2016: Nikolaos Ignatiadis.
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
module Spline =

    open MKLNET


    ///Helper functions
    ///  
    type ReinschQ(values:float[]) =    
        member m.h = values
        member m.Item(i,j) = 
                let h = values 
                if (i = j) then
                    1. / h.[i]
                elif (i = j+1) then
                    -1. / h.[j] - 1. / h.[j+1]
                elif (i = j+2) then
                    1. / h.[j+1]
                else
                    0.
        member m.Size = 
            let n = values.Length
            n+1,n-1


    //tests
    let rq = ReinschQ([|5.;5.;5.|])
    // should be 0.2
    rq.[0,0]
    rq.[1,1]
    // (should be -0.4)
    rq.[2,1]
    // (should be 4,2)
    rq.Size

    let multiplyInplace (out:float[]) (q:ReinschQ) (g:float[]) = 
        let n = out.Length
        if n <> (q.Size |> fst)  then failwith "DimensionMismatch"
        elif g.Length <> (q.Size |> snd) then failwith "DimensionMismatch"
        else 
            for i = 0 to n-1 do 
                out.[i] <- 0.
                for j = (max 0 (i-2)) to (min i (n-3)) do
                    out.[i] <- out.[i] + g.[j]*q.[i,j]
            out

    //tests
    let testmul1 = 
        let rq = ReinschQ([|5.;5.;5.|])
        let out = Array.init 4 (fun x -> 0.)
        let gTest = [|2.;2.|]
        // should be [|0.4; -0.4; -0.4; 0.4|]
        multiplyInplace out rq gTest    

    let testmul2 = 
        let rq = ReinschQ([|1. .. 10.|])
        let out = Array.init 11 (fun x -> 0.)
        let gTest = [|for i = 0 to 8 do yield 2.|]
        // should be [|2.0; -2.0; 1.110223025e-16; 1.110223025e-16; 0.0; -5.551115123e-17;
        // -5.551115123e-17; 0.0; 0.0; -0.2; 0.2|]
        multiplyInplace out rq gTest    


    let multiplyInplaceTransposedQ (out:float[]) (q:ReinschQ) (g:float[]) = 
        let n = out.Length
        if n <> (q.Size |> snd)  then failwith "DimensionMismatch"
        elif g.Length <> (q.Size |> fst) then failwith "DimensionMismatch"
        else 
            let h = q.h
            let mutable deltaGP1 = (g.[1] - g.[0])/h.[0]
            for i = 0 to out.Length-1 do 
                let deltaG = deltaGP1
                deltaGP1 <- (g.[i+2] - g.[i+1])/h.[i+1]
                out.[i] <- deltaGP1 - deltaG
            out 

    let testmultrans2 = 
        let rq = ReinschQ([|1. .. 10.|])
        let out = Array.init 9 (fun x -> 0.)
        let gTest = [|0. .. 10.|] |> Array.rev
        // should be   [|0.5; 0.1666666667; 0.08333333333; 0.05; 0.03333333333; 0.02380952381;
        // 0.01785714286; 0.01388888889; 0.01111111111|]
        multiplyInplaceTransposedQ out rq gTest    

    ///Helper functions
    ///  
    type ReinschR(values:float[]) =    
        member m.h = values
        member m.Item(i,j) = 
            let h = values 
            if (i=j) then
                (h.[i] + h.[i+1]) / 3.
            elif abs(i-j) = 1 then
                h.[max i j] / 6.
            else
                0.
        member m.Size = 
            let n = values.Length
            n-1,n-1

    //tests
    let reinschrTests =
        let rr = ReinschR([|5.;5.;5.|])
        // should be 0.2
        rr.[0,0]
        rr.[1,1]
        // (should be -0.4)
        rr.[2,1]
        // (should be 4,2)
        rr.Size

    ///
    let QtQpR (h:float[]) (alpha:float) (w:float []) = 
        let n = h.Length-1
        let Q = ReinschQ h
        let R = ReinschR h 
        let out = Matrix.zero 3 n
        // main diagonal
        for i=0 to n-1 do
            let fstTerm = (Q.[i+2,i] / w.[i+2] - Q.[i+1,i] / w.[i+1]) /  h.[i+1]
            let sndTerm = (Q.[i+1,i] / w.[i+1] - Q.[i,i] / w.[i]) / h.[i]
            out.[2,i] <- alpha * (fstTerm - sndTerm) + R.[i,i]
        // 1st superdiagonal
        for i=0 to n-2 do
            let fstTerm = (Q.[i+2,i+1]  / w.[i+2] - Q.[i+1,i+1] / w.[i+1]) / h.[i+1]
            let sndTerm =  Q.[i+1, i+1] / w.[i+1] / h.[i]
            out.[1,i+1] <- alpha * (fstTerm - sndTerm) + R.[i, i+1]
        // 2nd superdiagonal
        for i=0 to n-3 do
            let fstTerm = Q.[i+2,i+2] / w.[i+2] / h.[i+1]
            out.[0,i+2] <- alpha * fstTerm + R.[i,i+2]   
        out


    let QtQpRTests =
        let rq = [|1. .. 10.|]
        let out = Array.init (rq.Length+1) (fun x -> 1.)
        // should be 
            // matrix [[0.0; 0.0; 1.666666667; 0.8333333333; 0.5; 0.3333333333;
            //    0.2380952381; 0.1785714286; 0.1388888889]
            //   [0.0; -11.33333333; -4.222222222; -1.916666667; -0.8; -0.126984127;
            //    0.3418367347; 0.7033730159; 1.00308642]
            //   [36.0; 12.22222222; 7.472222222; 6.05; 5.688888889; 5.77324263;
            //    6.077806122; 6.503858025; 7.002469136]]
        QtQpR rq 10. out


    let QtQpRTests2 =
        let rq = [|1. .. 10.|]
        let out = Array.init (rq.Length+1) (fun x -> 1.)
        // should be 
            // matrix [[0.0; 0.0; 1.666666667; 0.8333333333; 0.5; 0.3333333333;
            //    0.2380952381; 0.1785714286; 0.1388888889]
            //   [0.0; -11.33333333; -4.222222222; -1.916666667; -0.8; -0.126984127;
            //    0.3418367347; 0.7033730159; 1.00308642]
            //   [36.0; 12.22222222; 7.472222222; 6.05; 5.688888889; 5.77324263;
            //    6.077806122; 6.503858025; 7.002469136]]
        QtQpR rq 10. out

    let diff (x:float []) = 
        Array.init (x.Length-1) (fun i -> x.[i+1] - x.[i])

    let diffTests = 
        let t = [|1.;20.;3.;4.;5.|] 
        diff t
      
    let initFullSparseRpαQtQ (RpαQtQ:matrix) =
        let indexSequence = 
            [|
                for m = 0 to 2 do
                    for n = 0 to RpαQtQ.NumCols-1 do 
                        if m = 2 then n,n,RpαQtQ.[m,n]
                        
                        if m = 1 && n > 0 then 
                            n-1,n,RpαQtQ.[m,n] 
                        if m = 0 && n > 1 then 
                            n-2,n,RpαQtQ.[m,n]
                        
                        if m = 1 && n > 0 then 
                            n,n-1,RpαQtQ.[m,n]                    
                        if m = 0 && n > 1 then 
                            n,n-2,RpαQtQ.[m,n] 
            |]
        Matrix.initSparse RpαQtQ.NumCols RpαQtQ.NumCols indexSequence
                        

    let initFullSparseRpαQtQTests =
        let rq = [|1. .. 10.|]
        let out = Array.init (rq.Length+1) (fun x -> 1.)
        // should be:
        // matrix [[36.0; -11.33333333; 1.666666667; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
        //           [0.0; 12.22222222; -4.222222222; 0.8333333333; 0.0; 0.0; 0.0; 0.0;
        //            0.0]
        //           [0.0; 0.0; 7.472222222; -1.916666667; 0.5; 0.0; 0.0; 0.0; 0.0]
        //           [0.0; 0.0; 0.0; 6.05; -0.8; 0.3333333333; 0.0; 0.0; 0.0]
        //           [0.0; 0.0; 0.0; 0.0; 5.688888889; -0.126984127; 0.2380952381; 0.0;
        //            0.0]
        //           [0.0; 0.0; 0.0; 0.0; 0.0; 5.77324263; 0.3418367347; 0.1785714286;
        //            0.0]
        //           [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 6.077806122; 0.7033730159;
        //            0.1388888889]
        //           [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 6.503858025; 1.00308642]
        //           [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 7.002469136]]
        QtQpR rq 10. out
        |> initFullSparseRpαQtQ


    // expects row major matrix
    let solveLeastSquares (m:matrix) (y:vector) = 
        let mFlat = 
            let flat = m |> Matrix.toJaggedArray |> Array.concat
            let tmp = System.GC.AllocateArray<float>(flat.Length,pinned=true) 
            flat |> Array.iteri (fun i x -> tmp.[i] <- x)
            tmp
        // Computes the Cholesky factorization of a symmetric (Hermitian) positive-definite band matrix.
        let res = MKLNET.Lapack.pbtrf(Layout.RowMajor,UpLoChar.Upper,m.NumCols,2,mFlat,max 1 (m.NumCols))
        let yFlat =
            let tmp = System.GC.AllocateArray<float>(y.Length,pinned=true) 
            y |> Seq.iteri (fun i x -> tmp.[i] <- x)
            tmp 
        let mFlat' = 
            let tmp = System.GC.AllocateArray<float>(mFlat.Length,pinned=true) 
            mFlat |> Array.iteri (fun i x -> tmp.[i] <- x)
            tmp
        let resT = MKLNET.Lapack.pbtrs(Layout.RowMajor,UpLoChar.Upper,m.NumCols, 2, 1,mFlat',m.NumCols,yFlat,1)
        yFlat

        
    let fitSplineSparseLapack (x:float []) (y:float []) lambda = 
        let ws = Array.create x.Length 1.
        let diffx = diff x 
        let RpαQtQ = QtQpR diffx lambda ws
        let h = diff x
        let Q = ReinschQ(h)
        // fitSpline xd yd 2. should yield: gamma = [|-0.1; 0.1; 0.4; -0.1; -0.4; 0.2; 0.2; -0.4|]
        let gamma = 
            multiplyInplaceTransposedQ (Array.zeroCreate (x.Length-2)) Q y
            |> vector 
            |> solveLeastSquares RpαQtQ
        let g = 
            multiplyInplace (Array.zeroCreate (x.Length)) Q (gamma) 
            |> Array.map2 (fun ws g -> g / ws) ws
            |> Array.map (fun g -> g * lambda)
            |> Array.map2 (fun y g -> y - g ) y
        g


    let initPredict (gamma:float[]) (g:float[]) (x:float[]) (xToPred:float) = 
        let n = x.Length-1
        let idxL = 
            match Array.tryFindIndexBack (fun v -> v <= xToPred) x with
            | Some i -> i
            | None -> -1 
        let idxR = idxL + 1
        if idxL = -1 then 
            let gl = g.[0]
            let gr = g.[1]
            let gamma' = gamma.[0]
            let xl = x.[0]
            let xr = x.[1]
            let gprime = (gr-gl) / (xr-xl) - 1. / 6. * (xr-xl) * gamma'
            gl - (xl-xToPred) * gprime
        elif idxL = n then  
            let gl = g.[n-1]
            let gr = g.[n]
            let gamma' = gamma.[n-2]
            let xl = x.[n-1]
            let xr = x.[n]
            let gprime = (gr-gl) / (xr-xl) + 1./ 6. * (xr-xl) * gamma'
            gr + (xToPred - xr)*gprime
        else 
            let xl = x.[idxL]
            let xr = x.[idxR]
            let gammaL = if idxL = 0 then 0. else gamma.[idxL-1]
            let gammaR = if idxL = n-1 then 0. else gamma.[idxR-1]
            let gl = g.[idxL]
            let gr = g.[idxR]
            let h = xr-xl
            let tmp = ((xToPred-xl) * gr + (xr-xToPred) * gl) / h
            tmp - (1. / 6. * (xToPred-xl) * (xr-xToPred) * ((1. + (xToPred-xl) / h) * gammaR + (1.+ (xr-xToPred) / h) * gammaL))

    
    let fitSplineSparseLapack' (x:float []) (y:float []) lambda = 
        let ws = Array.create x.Length 1.
        let diffx = diff x 
        let RpαQtQ = QtQpR diffx lambda ws
        let h = diff x
        let Q = ReinschQ(h)
        // fitSpline xd yd 2. should yield: gamma = [|-0.1; 0.1; 0.4; -0.1; -0.4; 0.2; 0.2; -0.4|]
        let gamma = 
            multiplyInplaceTransposedQ (Array.zeroCreate (x.Length-2)) Q y
            |> vector 
            |> solveLeastSquares RpαQtQ
        let g = 
            multiplyInplace (Array.zeroCreate (x.Length)) Q (gamma) 
            |> Array.map2 (fun ws g -> g / ws) ws
            |> Array.map (fun g -> g * lambda)
            |> Array.map2 (fun y g -> y - g ) y
        initPredict gamma g x 

module QuantBasedAlignment = 

    ///
    let paletteArray =
        [|
            Blue1     
            Blue2     
            Blue3     
            Red1      
            Red2      
            Red3      
            LightBlue1
            LightBlue2
            LightBlue3
            LightRed1 
            LightRed2 
            LightRed3 
            Green1    
            Green2    
            Green3    
            Orange1   
            Orange2   
            Orange3   
            Cyan1     
            Cyan2     
            Cyan3     
            Magenta1  
            Magenta2  
            Magenta3  
        |]

    /// Define the random number generator outside of a potential loop.
    let getRandomColor (rnd: System.Random) =
        let index = rnd.Next(0,23)
        paletteArray.[index]

    // ///
    // let downcastPipeline (x : IEstimator<_>) = 
    //     match x with 
    //     | :? IEstimator<ITransformer> as y -> y
    //     | _ -> failwith "downcastPipeline: expecting a IEstimator<ITransformer>"

    ///
    type PeptideIon = 
        {
            Sequence             : string
            GlobalMod            : int
            Charge               : int        
        }
    
    ///
    let toPeptideIon (qp:QuantificationResult) = 
        {
            Sequence             = qp.StringSequence
            GlobalMod            = qp.GlobalMod
            Charge               = qp.Charge
        }

    ///
    type AlignmentFile = {
        FileName                   : string
        QuantifiedPeptides         : Map<PeptideIon,QuantificationResult> 
        MissingPeptides            : PeptideIon []
        GainedPeptides             : AlignmentResult []
        }
    

    [<CLIMutable>]
    type PeptideForLearning = 
        {
            // [<ColumnName("Sequence")>]
            Sequence                     : string
            // [<ColumnName("GlobalMod")>]
            GlobalMod                    : string
            // [<ColumnName("Charge")>]
            Charge                       : string
            // [<ColumnName("PepSequenceID")>]
            PepSequenceID                : string
            // [<ColumnName("ModSequenceID")>]
            ModSequenceID                : string
            // [<ColumnName("SourceScanTime")>]
            SourceScanTime               : float32
            // [<ColumnName("TargetScanTime")>]
            TargetScanTime               : float32
            // [<ColumnName("ScanTimeDifference")>]
            ScanTimeDifference           : float32
        }
    
    [<CLIMutable>]
    type PeptideComplete = 
        {
            // [<ColumnName("Sequence")>]
            Sequence                     : string
            // [<ColumnName("GlobalMod")>]
            GlobalMod                    : string
            // [<ColumnName("Charge")>]
            Charge                       : string
            // [<ColumnName("PepSequenceID")>]
            PepSequenceID                : string
            // [<ColumnName("ModSequenceID")>]
            ModSequenceID                : string
            // [<ColumnName("SourceScanTime")>]
            SourceScanTime               : float32
            // [<ColumnName("SourceIntensity")>]
            SourceIntensity              : float32
            // [<ColumnName("SourceStabw")>]
            SourceStabw                  : float32
            // [<ColumnName("TargetScanTime")>]
            TargetScanTime               : float32
            // [<ColumnName("TargetIntensity")>]
            TargetIntensity              : float32
            // [<ColumnName("RtTrace_SourceFile")>]
            RtTrace_SourceFile                              : float [] 
            // [<ColumnName("IntensityTrace_SourceFile")>]
            IntensityTrace_SourceFile                       : float []
            // [<ColumnName("RtTrace_TargetFile")>]
            RtTrace_TargetFile                              : float []
            // [<ColumnName("IntensityTrace_TargetFile")>]
            IntensityTrace_TargetFile                       : float []          
            // [<ColumnName("IsotopicPatternMz_SourceFile")>]
            IsotopicPatternMz_SourceFile                    : float []          
            // [<ColumnName("IsotopicPatternIntensity_Observed_SourceFile")>]
            IsotopicPatternIntensity_Observed_SourceFile    : float []         
            // [<ColumnName("IsotopicPatternMz_TargetFile")>]
            IsotopicPatternMz_TargetFile                    : float []         
            // [<ColumnName("IsotopicPatternIntensity_Observed_TargetFile")>]
            IsotopicPatternIntensity_Observed_TargetFile    : float []
        }

    /////
    //let formatString s = String.filter (fun x -> Char.IsUpper x) s

    ///
    let toPeptideForLearning (targetPep:QuantificationResult option) (sourcePep:QuantificationResult) = 
        match targetPep with 
        | Some tP -> 
            {
                Sequence                                        = (*formatString*) sourcePep.StringSequence
                GlobalMod                                       = sourcePep.GlobalMod          |> string      
                Charge                                          = sourcePep.Charge             |> string
                PepSequenceID                                   = sourcePep.PepSequenceID      |> string
                ModSequenceID                                   = sourcePep.ModSequenceID      |> string
                SourceScanTime                                  = getTargetScanTime sourcePep  |> float32
                TargetScanTime                                  = getTargetScanTime tP         |> float32
                ScanTimeDifference                              = ((getTargetScanTime tP) - (getTargetScanTime sourcePep)) |> float32
            },
            {
                Sequence                                        = (*formatString*) sourcePep.StringSequence
                GlobalMod                                       = sourcePep.GlobalMod          |> string      
                Charge                                          = sourcePep.Charge             |> string
                PepSequenceID                                   = sourcePep.PepSequenceID      |> string
                ModSequenceID                                   = sourcePep.ModSequenceID      |> string
                SourceScanTime                                  = getTargetScanTime sourcePep  |> float32
                SourceIntensity                                 = getTargetIntensity sourcePep |> float32
                SourceStabw                                     = getTargetStabw sourcePep     |> float32
                TargetScanTime                                  = getTargetScanTime tP         |> float32
                TargetIntensity                                 = getTargetIntensity tP        |> float32
                RtTrace_SourceFile                              = getTargetRtTrace sourcePep
                IntensityTrace_SourceFile                       = getTargetIntensityTrace sourcePep
                RtTrace_TargetFile                              = getTargetRtTrace tP
                IntensityTrace_TargetFile                       = getTargetIntensityTrace tP
                IsotopicPatternMz_SourceFile                    = getIsotopicPatternMz sourcePep
                IsotopicPatternIntensity_Observed_SourceFile    = getIsotopicPatternIntensity_Observed sourcePep
                IsotopicPatternMz_TargetFile                    = getIsotopicPatternMz tP  
                IsotopicPatternIntensity_Observed_TargetFile    = getIsotopicPatternIntensity_Observed tP
            }
        | None -> 
            {
                Sequence                                        = (*formatString*) sourcePep.StringSequence
                GlobalMod                                       = sourcePep.GlobalMod          |> string   
                Charge                                          = sourcePep.Charge             |> string
                PepSequenceID                                   = sourcePep.PepSequenceID      |> string
                ModSequenceID                                   = sourcePep.ModSequenceID      |> string
                SourceScanTime                                  = getTargetScanTime sourcePep  |> float32
                TargetScanTime                                  = nan                          |> float32
                ScanTimeDifference                              = nan                          |> float32
            },
            {
                Sequence                                        = (*formatString*) sourcePep.StringSequence
                GlobalMod                                       = sourcePep.GlobalMod          |> string   
                Charge                                          = sourcePep.Charge             |> string
                PepSequenceID                                   = sourcePep.PepSequenceID      |> string
                ModSequenceID                                   = sourcePep.ModSequenceID      |> string
                SourceScanTime                                  = getTargetScanTime sourcePep  |> float32
                SourceIntensity                                 = getTargetIntensity sourcePep |> float32
                SourceStabw                                     = getTargetStabw sourcePep     |> float32
                TargetScanTime                                  = nan                          |> float32
                TargetIntensity                                 = nan                          |> float32
                RtTrace_SourceFile                              = getTargetRtTrace sourcePep
                IntensityTrace_SourceFile                       = getTargetIntensityTrace sourcePep
                RtTrace_TargetFile                              = [||]
                IntensityTrace_TargetFile                       = [||]
                IsotopicPatternMz_SourceFile                    = getIsotopicPatternMz sourcePep
                IsotopicPatternIntensity_Observed_SourceFile    = getIsotopicPatternIntensity_Observed sourcePep
                IsotopicPatternMz_TargetFile                    = [||]
                IsotopicPatternIntensity_Observed_TargetFile    = [||]
            }

    ///
    [<CLIMutable>]
    type ScanTimePrediction = 
        {
            // [<ColumnName("Score")>]
            TargetScanTime : float32
        }

    ///
    let getQuantifiedPeptides (quantFilePath:string) = 
        ///
        let peptides =
            Csv.CsvReader<QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(quantFilePath,'\t',false,1)
            |> Array.ofSeq
        let filteredPeptides =
            peptides
            |> Array.filter (fun qp -> 
                // Filter for peptides where fit assures a good estimation of the ScanTime
                (qp.GlobalMod = 0 && qp.Params_Light |> Array.isEmpty |> not) || (qp.GlobalMod = 1 && qp.Params_Heavy |> Array.isEmpty |> not)
                )
            |> Array.filter (fun qp -> tryTargetGetScanTime qp |> Option.isSome)
            |> Array.filter (fun qp -> (getTargetScanTimeDifference qp |> abs) / (getTargetStabw qp) < 2. )          
        filteredPeptides 
    
    /// Reads quant files, filters for quality quantifications and creates AlignmentFiles.
    let createAlignmentFiles (sourceFiles :string []) (targetFile:string) = 
        let toAlignmentFile (allPepIons:PeptideIon []) (filePath:string) peptides = 
            let presentPeptides = 
                peptides 
                |> Array.map (fun qp -> toPeptideIon qp, qp)
                |> Map.ofArray
            let missingPeptides = 
                allPepIons
                |> Array.filter (fun pep -> 
                    let pepUnMod = {pep with GlobalMod = 0}
                    let pepMod = {pep with GlobalMod = 1}
                    (presentPeptides.ContainsKey pepMod || presentPeptides.ContainsKey pepUnMod) 
                    |> not
                    )
            {
                FileName                   = Path.GetFileNameWithoutExtension(filePath)
                QuantifiedPeptides         = presentPeptides 
                MissingPeptides            = missingPeptides 
                GainedPeptides             = [||]
            }
        let targetPeptides = 
            targetFile 
            |> getQuantifiedPeptides
        let sourcePeptides = 
            sourceFiles
            |> Array.map getQuantifiedPeptides
        let allPeptideIons =
            sourcePeptides
            |> Array.append [|targetPeptides|] 
            |> Array.concat
            |> Array.map (fun qp -> toPeptideIon qp)
            |> Array.distinct
        let targetAlignmentFile = toAlignmentFile allPeptideIons targetFile targetPeptides
        let sourceAlignmentFiles = Array.map2 (toAlignmentFile allPeptideIons) sourceFiles sourcePeptides
        sourceAlignmentFiles, targetAlignmentFile 
    
    ///
    let createAlignmentResult (quantifiedPeptide:QuantificationResult) (scanTimePrediction:ScanTimePrediction) = 
        {
            StringSequence                                  = quantifiedPeptide.StringSequence
            GlobalMod                                       = quantifiedPeptide.GlobalMod
            Charge                                          = quantifiedPeptide.Charge
            PepSequenceID                                   = quantifiedPeptide.PepSequenceID
            ModSequenceID                                   = quantifiedPeptide.ModSequenceID
            Mz                                              = Mass.toMZ (quantifiedPeptide.TheoMass) (float quantifiedPeptide.Charge) 
            ProteinNames                                    = quantifiedPeptide.ProteinNames
            PredictedScanTime                               = float scanTimePrediction.TargetScanTime
            ScanTime_SourceFile                             = getTargetScanTime quantifiedPeptide
            ApexIntensity_SourceFile                        =
                if quantifiedPeptide.GlobalMod = 1 then
                    quantifiedPeptide.MeasuredApex_Heavy
                else
                    quantifiedPeptide.MeasuredApex_Light
            Quant_SourceFile                                =
                if quantifiedPeptide.GlobalMod = 1 then
                    quantifiedPeptide.Quant_Heavy
                else
                    quantifiedPeptide.Quant_Light
            RtTrace_SourceFile                              = getTargetRtTrace quantifiedPeptide
            IntensityTrace_SourceFile                       = getTargetIntensityTrace quantifiedPeptide
            IsotopicPatternMz_SourceFile                    = getIsotopicPatternMz quantifiedPeptide       
            IsotopicPatternIntensity_Observed_SourceFile    = getIsotopicPatternIntensity_Observed quantifiedPeptide       
        }   


    ///
    type ModelMetrics = 
        {
        RSquared                             : float
        Sequence                            : string []
        GlobalMod                           : int []
        Charge                              : int []
        PepSequenceID                       : int []
        ModSequenceID                       : int []
        X_Intensities                       : float []
        X_Stabw                             : float []        
        X_Test                              : float []
        X_IsotopicPatternMz                 : float [][]
        X_IsotopicPatternIntensity_Observed : float [][]
        X_RtTrace                           : float [][]
        X_IntensityTrace                    : float [][]   
        Y_Test                              : float []
        YHat_Test                           : float []
        YHat_Refined_Test                   : float []
        Y_Intensities                       : float []
        Y_IsotopicPatternMz                 : float [][]
        Y_IsotopicPatternIntensity_Observed : float [][]
        Y_RtTrace                           : float [][]
        Y_IntensityTrace                    : float [][]
        DtwDistanceBefore                   : float []
        }



    ///
    let createMetricsChart fileName (*stabwMedian*) (rnd:System.Random) (metrics:ModelMetrics) = 
        ///
        let traceName = ""//(sprintf "#TestPeptides:%i Rsquared:%f RMS:%f, %s " metrics.X_Test.Length metrics.Metrics.RSquared metrics.Metrics.MeanAbsoluteError fileName)
        let color = getRandomColor rnd |> FSharpAux.Colors.toWebColor        
        let xVsY = 
            [
            Chart.Point(metrics.X_Test,metrics.Y_Test) |> Chart.withMarkerStyle(Color = color)
            Chart.Line(Array.zip metrics.X_Test metrics.YHat_Test |> Array.sortBy fst)
            ]
            |> Chart.Combine
            |> Chart.withX_AxisStyle("source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("measured (dotted) and predicted (Line) target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let yVsYHat = 
            Chart.Point(metrics.Y_Test,metrics.YHat_Test)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("target ScanTimes (Y_Test)")
            |> Chart.withY_AxisStyle("predicted target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let xVsDifferenceYandYHat = 
            Chart.Point(metrics.X_Test, Array.map2 (fun y yHat -> y - yHat) metrics.Y_Test metrics.YHat_Test)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("target ScanTimes (Y_Test) - predicted target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let xVsDifferenceYandYHatNormed = 
            Chart.Point(metrics.X_Test, Array.map3 (fun y yHat stabwMedian -> (y - yHat) / stabwMedian) metrics.Y_Test metrics.YHat_Test metrics.X_Stabw)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("normed target ScanTimes (Y_Test) - predicted target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let xVsDifferenceYandYHat_refined = 
            Chart.Point(metrics.X_Test, Array.map2 (fun y yHat -> y - yHat) metrics.Y_Test metrics.YHat_Refined_Test)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("target ScanTimes (Y_Test) - refined target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let xVsDifferenceYandYHatNormed_refined = 
            Chart.Point(metrics.X_Test, Array.map3 (fun y yHat stabwMedian -> (y - yHat) / stabwMedian) metrics.Y_Test metrics.YHat_Refined_Test metrics.X_Stabw)
            |> Chart.withMarkerStyle(Color = color)
            |> Chart.withX_AxisStyle("Source ScanTimes (X_Test)")
            |> Chart.withY_AxisStyle("normed target ScanTimes (Y_Test) - refined target ScanTimes (YHat_Test)")
            |> Chart.withTraceName traceName
        let makeBar yHat_Test = 
            let metric =  Array.map3 (fun y yHat stabwMedian -> (y - yHat) / stabwMedian) metrics.Y_Test yHat_Test metrics.X_Stabw
            let UpToOne   = metric |> Array.filter (fun x -> abs x >= 0. && abs x < 1.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToTwo   = metric |> Array.filter (fun x -> abs x >= 1. && abs x < 2.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToThree = metric |> Array.filter (fun x -> abs x >= 2. && abs x < 3.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToFour  = metric |> Array.filter (fun x -> abs x >= 3. && abs x < 4.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToFive  = metric |> Array.filter (fun x -> abs x >= 4. && abs x < 5.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToSix   = metric |> Array.filter (fun x -> abs x >= 5. && abs x < 6.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let UpToSeven = metric |> Array.filter (fun x -> abs x >= 6. && abs x < 7.)   |> Array.length |> float |> fun x -> x / float metric.Length
            let OutOfSeven  = metric |> Array.filter (fun x -> abs x > 7.)                |> Array.length |> float |> fun x -> x / float metric.Length           
            [
                Chart.StackedColumn(values=[UpToOne],keys=[traceName],Name="UpToOne")       |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Green1)
                Chart.StackedColumn(values=[UpToTwo],keys=[traceName],Name="UpToTwo")       |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Green2)
                Chart.StackedColumn(values=[UpToThree],keys=[traceName],Name="UpToThree")   |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Green3)
                Chart.StackedColumn(values=[UpToFour],keys=[traceName],Name="UpToFour")     |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Orange1)
                Chart.StackedColumn(values=[UpToFive],keys=[traceName],Name="UpToFive")     |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Orange2)
                Chart.StackedColumn(values=[UpToSix],keys=[traceName],Name="UpToSix")       |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Orange3)
                Chart.StackedColumn(values=[UpToSeven],keys=[traceName],Name="UpToSeven")   |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Red1)
                Chart.StackedColumn(values=[OutOfSeven],keys=[traceName],Name="OutOfSeven") |> Chart.withMarkerStyle(Color=FSharpAux.Colors.toWebColor FSharpAux.Colors.Table.StatisticalGraphics24.Red2)
            ]
            |> Chart.Combine
            |> Chart.withY_AxisStyle("")
        let bar = makeBar metrics.YHat_Test
        let barRef = makeBar metrics.YHat_Refined_Test
        [
        xVsY
        yVsYHat
        xVsDifferenceYandYHat
        xVsDifferenceYandYHatNormed
        xVsDifferenceYandYHat_refined
        xVsDifferenceYandYHatNormed_refined        
        bar
        barRef
        ]
        |> Chart.Stack 2
        |> Chart.withSize(2000.,2500.)


    ///
    let saveMetrics outputDir targetFileName sourceFileName (metrics:ModelMetrics) = 
        let fileName = (targetFileName ) + ".alignmetric"
        let outFilePath = Path.Combine [|outputDir;fileName|]
        let data = 
            [|
            for i = 1 to metrics.X_Intensities.Length-1 do 
                yield
                    {           
                        Sequence                             = metrics.Sequence.[i]
                        GlobalMod                            = metrics.GlobalMod.[i]
                        Charge                               = metrics.Charge.[i]
                        PepSequenceID                        = metrics.PepSequenceID.[i]
                        ModSequenceID                        = metrics.ModSequenceID.[i]
                        X_FileName                           = sourceFileName 
                        X_Intensities                        = metrics.X_Intensities.[i] 
                        X_Stabw                              = metrics.X_Stabw.[i]
                        X_Test                               = metrics.X_Test.[i] 
                        X_IsotopicPatternMz                  = metrics.X_IsotopicPatternMz.[i]
                        X_IsotopicPatternIntensity_Observed  = metrics.X_IsotopicPatternIntensity_Observed.[i]
                        X_RtTrace                            = metrics.X_RtTrace.[i]
                        X_IntensityTrace                     = metrics.X_IntensityTrace.[i]
                        Y_Test                               = metrics.Y_Test.[i] 
                        YHat_Test                            = metrics.YHat_Test.[i] 
                        YHat_Refined_Test                    = metrics.YHat_Refined_Test.[i]
                        Y_Intensity                          = metrics.Y_Intensities.[i]
                        Y_IsotopicPatternMz                  = metrics.Y_IsotopicPatternMz.[i]
                        Y_IsotopicPatternIntensity_Observed  = metrics.Y_IsotopicPatternIntensity_Observed.[i]
                        Y_RtTrace                            = metrics.Y_RtTrace.[i]
                        Y_IntensityTrace                     = metrics.Y_IntensityTrace.[i]   
                        DtwDistanceBefore                    = metrics.DtwDistanceBefore.[i]
                    }
            |]
        if System.IO.File.Exists outFilePath then 
            data
            |> SeqIO'.csv "\t" false false
            |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        else
            data
            |> SeqIO'.csv "\t" true false
            |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)

    ///
    // let initAlign (logger:NLog.Logger) (ctx:MLContext) (pepsForLearning: (PeptideForLearning*PeptideComplete) []) = 
    let initAlign (logger:NLog.Logger) (pepsForLearning: (PeptideForLearning*PeptideComplete) []) = 
        logger.Trace ("Sampling test and train data")
        let pepsForLearning = pepsForLearning |> Array.distinctBy (fun (pepLearning,pepComp) -> pepLearning.SourceScanTime)
        let train,test = 
            let getBinIdx width scantime = int ((scantime / width))    
            let _,train,test = 
                pepsForLearning
                |> Array.groupBy (fun (pl,_) -> getBinIdx 10. (pl.SourceScanTime |> float))
                |> Array.map (fun (binIdx,ions) -> 
                    let dataS = ions |> Array.shuffleFisherYates
                    let knotX = 
                        dataS 
                        |> Array.map fst
                        |> Array.maxBy (fun x -> x.SourceScanTime)
                    let train,test =
                        let nTest = 
                            (float dataS.Length) * 0.9
                            |> int
                        dataS.[.. nTest], dataS.[nTest+1 ..] 
                    float knotX.SourceScanTime, train, test
                    )
                |> Array.unzip3
            train |> Array.concat |> Array.sortBy (fun (pepLearning,pepComp) -> pepLearning.SourceScanTime), 
            test |> Array.concat |> Array.sortBy (fun (pepLearning,pepComp) -> pepLearning.SourceScanTime)
        let trainLearn,trainComp = train |> Array.unzip
        let testLearn,testComp = test |> Array.unzip    
        logger.Trace ("Training Spline")
        let trainer lambda = 
            let fit = Spline.fitSplineSparseLapack' (train |> Array.map (fun (x,y) -> float x.SourceScanTime)) (train |> Array.map (fun (x,y) -> float x.TargetScanTime)) lambda 
            let x,y,yHat = 
                test
                |> Array.map (fun (x,y) -> 
                    x, x.TargetScanTime, fit (float x.SourceScanTime)
                    )
                |> Array.unzip3
            let rS = FSharp.Stats.Fitting.GoodnessOfFit.calculateDeterminationFromValue yHat (Array.map float y)
            lambda, rS, fit                  
        logger.Trace ("Optimizing Lambdas")
        let models =      
            [|1. .. 1. .. 200.|]
            |> Array.map trainer
        let lambda,rSquared,model = 
            models
            |> Array.maxBy (fun (x,y,z) -> y )                  
        logger.Trace (sprintf "Optimizing Lambdas:Finished, selected lamda:%f" lambda)
        logger.Trace ("Training Spline:Finished")
        // [
        //     test
        //     |> Array.map (fun (x,y) -> 
        //             x.SourceScanTime, (float x.TargetScanTime)
        //             )
        //     |> Chart.Point    
        //     models
        //     |> Array.map (fun (x,y,z) -> 
        //         test
        //         |> Array.map (fun (x,y) -> 
        //                 x.SourceScanTime, z (float x.SourceScanTime)
        //                 )
        //         |> Chart.Line
        //         |> Chart.withTraceName (sprintf "lambda %f" x)
        //         )
        //     |> Chart.Combine
        // ]
        // |> Chart.Combine
        // // |> Chart.Show
        // models
        // |> Array.map (fun (x,y,z) -> x,y)
        // |> Chart.Point
        // // |> Chart.Show
        let metrics =             
            let rSquared = rSquared
            let yHat          = testComp |> Seq.map (fun x -> model (float x.SourceScanTime))|> Array.ofSeq
            let sequence      = testComp |> Seq.map (fun x -> x.Sequence)       |> Array.ofSeq
            let globalMod     = testComp |> Seq.map (fun x -> x.GlobalMod)      |> Array.ofSeq
            let charge        = testComp |> Seq.map (fun x -> x.Charge)         |> Array.ofSeq
            let pepSequenceID = testComp |> Seq.map (fun x -> x.PepSequenceID)  |> Array.ofSeq
            let modSequenceID = testComp |> Seq.map (fun x -> x.ModSequenceID)  |> Array.ofSeq
            let i             =  testComp |> Seq.map (fun x -> float x.SourceIntensity)  |> Array.ofSeq
            let std           =  testComp |> Seq.map (fun x -> float x.SourceStabw)      |> Array.ofSeq
            let x             =  testComp |> Seq.map (fun x -> float x.SourceScanTime)   |> Array.ofSeq
            let y             =  testComp |> Seq.map (fun x -> float x.TargetScanTime)   |> Array.ofSeq
            let yIntensities  =  testComp |> Seq.map (fun x -> float x.TargetIntensity)  |> Array.ofSeq
            let xSource = testComp |> Seq.map (fun x -> Array.ofSeq x.RtTrace_SourceFile)         |> Array.ofSeq
            let ySource = testComp |> Seq.map (fun x -> Array.ofSeq x.IntensityTrace_SourceFile)  |> Array.ofSeq
            let xTarget = testComp |> Seq.map (fun x -> Array.ofSeq x.RtTrace_TargetFile)         |> Array.ofSeq
            let yTarget = testComp |> Seq.map (fun x -> Array.ofSeq x.IntensityTrace_TargetFile)  |> Array.ofSeq     
            let yHatAfterRefinement,dtwDistanceBefore,dtwDistanceAfter = 
                [|
                    for i = 0 to xSource.Length-1 do                         
                        let target = Array.zip xTarget.[i] (DTW'.zNorm yTarget.[i])
                        let source = Array.zip xSource.[i] (DTW'.zNorm ySource.[i])
                        let yRefined = 
                            DTW'.align' target source x.[i] 
                            |> snd
                        let alignment = 
                            DTW'.align target source 
                            |> Array.ofList
                        let dtwDistanceBefore = 
                            DTW'.distance None None None None None None (Array.map snd target)  (Array.map snd source)
                        let dtwDistanceAfter = 
                            DTW'.distance None None None None None None (Array.map snd target)  (Array.map snd alignment)
                        yield yRefined,dtwDistanceBefore,dtwDistanceAfter
                |]
                |> Array.unzip3

            let x_IsotopicPatternMz                 = testComp |> Seq.map (fun x -> Array.ofSeq x.IsotopicPatternMz_SourceFile)                 |> Array.ofSeq
            let x_IsotopicPatternIntensity_Observed = testComp |> Seq.map (fun x -> Array.ofSeq x.IsotopicPatternIntensity_Observed_SourceFile) |> Array.ofSeq
            let y_IsotopicPatternMz                 = testComp |> Seq.map (fun x -> Array.ofSeq x.IsotopicPatternMz_TargetFile)                 |> Array.ofSeq
            let y_IsotopicPatternIntensity_Observed = testComp |> Seq.map (fun x -> Array.ofSeq x.IsotopicPatternIntensity_Observed_TargetFile) |> Array.ofSeq       
            {
                RSquared                            = rSquared
                Sequence                            = sequence                           
                GlobalMod                           = globalMod     |> Array.map int                     
                Charge                              = charge        |> Array.map int                     
                PepSequenceID                       = pepSequenceID |> Array.map int                     
                ModSequenceID                       = modSequenceID |> Array.map int                     
                X_Intensities                       = i
                X_Stabw                             = std      
                X_Test                              = x
                X_IsotopicPatternMz                 = x_IsotopicPatternMz
                X_IsotopicPatternIntensity_Observed = x_IsotopicPatternIntensity_Observed
                X_RtTrace                           = xSource
                X_IntensityTrace                    = ySource
                Y_Test                              = y
                YHat_Test                           = yHat
                YHat_Refined_Test                   = yHatAfterRefinement
                Y_Intensities                       = yIntensities
                Y_IsotopicPatternMz                 = y_IsotopicPatternMz
                Y_IsotopicPatternIntensity_Observed = y_IsotopicPatternIntensity_Observed
                Y_RtTrace                           = xTarget
                Y_IntensityTrace                    = yTarget
                DtwDistanceBefore                   = dtwDistanceBefore
            }

        // Spline
        let predict quantifiedPeptide = 
            quantifiedPeptide
            |> toPeptideForLearning None
            |> fun (pepToLearn,pepComp) -> 
                let v = model (float pepToLearn.SourceScanTime)
                {
                TargetScanTime = float32 v
                }
            |> createAlignmentResult quantifiedPeptide

        metrics, predict

    type Alignment = {
        Metrics      : ModelMetrics
        MetricsChart : GenericChart.GenericChart
        AlignFunc    : (QuantificationResult->AlignmentResult)
        SourceFile   : AlignmentFile
        }
    ///
    let performAlignment outDir rnd align (target: AlignmentFile) (source: AlignmentFile) =    
        ///
        let peptidesForLearning = 
            target.QuantifiedPeptides
            |> Seq.choose (fun tarQP -> 
                match Map.tryFind tarQP.Key source.QuantifiedPeptides with 
                | Some sourceQP ->
                    toPeptideForLearning (Some tarQP.Value) sourceQP
                    |> Some
                | None -> None
                ) 
            |> Array.ofSeq
            |> Array.groupBy (fun (x,y) -> x.SourceScanTime)
            |> Array.map (fun (_,(d)) -> d |> Array.maxBy (fun (x,y) -> y.SourceIntensity))
        ///
        let metrics,model: ModelMetrics*(QuantificationResult->AlignmentResult) = 
            align peptidesForLearning

        {
        Metrics      = metrics
        MetricsChart = createMetricsChart source.FileName (*stabwMedian*) rnd metrics
        AlignFunc    = model
        SourceFile   = source
        }
        
    ///
    let alignFiles diagCharts (processParams:AlignmentParams) (outputDir:string) (sourceFiles:string []) (targetFile:string) = 
        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension targetFile)
        logger.Trace (sprintf "target file: %A" targetFile)
        logger.Trace (sprintf "source files: %A" sourceFiles)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" processParams)
        let getPlotFilePathFilePath (plotName:string) (fileName:string) =
            let fileName = (Path.GetFileNameWithoutExtension fileName) + "_" + plotName 
            Path.Combine [|outputDir;fileName|]
        logger.Trace "Init align function"
        // let ctx = ML.MLContext()
        let rnd = System.Random()
        let align = initAlign logger 
        logger.Trace "Init align function: finished"
         
        logger.Trace "Reading and preparing .quant files for alignment"
        let sourceAlignmentFiles, targetAlignmentFile  = createAlignmentFiles sourceFiles targetFile
        logger.Trace "Reading and preparing .quant files for alignment: finished"

        logger.Trace "Performing Alignments"
        let alignmentsSortedByQuality = 
            let alignments =
                sourceAlignmentFiles
                |> Array.map (fun (sourceFile)  -> 
                    logger.Trace (sprintf "Performing Alignment %s vs %s" targetAlignmentFile.FileName sourceFile.FileName)
                    let alignResult = performAlignment outputDir rnd align targetAlignmentFile sourceFile
                    saveMetrics outputDir targetAlignmentFile.FileName sourceFile.FileName alignResult.Metrics
                    logger.Trace (sprintf "Performing Alignment %s vs %s: finished" targetAlignmentFile.FileName sourceFile.FileName)
                    alignResult
                    ) 
            let sortedByQuality = 
                alignments
                |> Array.sortByDescending (fun x -> x.Metrics.RSquared)
            if diagCharts then 
                sortedByQuality
                |> Array.map (fun x -> x.MetricsChart)
                |> Chart.Combine
                |> Chart.withTitle(targetAlignmentFile.FileName)
                |> Chart.SaveHtmlAs(getPlotFilePathFilePath "Metrics" targetAlignmentFile.FileName)                    
            sortedByQuality 
        logger.Trace "Performing Alignments: finished"
        logger.Trace "Transfer identifications"
        let result = 
            alignmentsSortedByQuality
            |> Array.fold (fun target alignment -> 
                    let peptideIonsToTransfer,peptideIonsStillMissing =
                        Array.append target.MissingPeptides (target.QuantifiedPeptides |> Map.toArray |> Array.map fst)
                        |> Array.fold (fun (pepsToTransfer,stillMissingPeps) missingPep -> 
                            match Map.tryFind missingPep alignment.SourceFile.QuantifiedPeptides with 
                            | Some pepToTransfer -> (pepToTransfer::pepsToTransfer,stillMissingPeps)
                            | None -> (pepsToTransfer,missingPep::stillMissingPeps)
                            ) ([],[])
                    let alignmentResults : AlignmentResult [] = 
                        peptideIonsToTransfer 
                        |> List.toArray 
                        |> Array.map alignment.AlignFunc 
                    {target with MissingPeptides = peptideIonsStillMissing |> Array.ofList; GainedPeptides = Array.append target.GainedPeptides alignmentResults}
                ) targetAlignmentFile
        logger.Trace "Writing Results"
        let outFilePath =
            let fileName = (result.FileName) + ".align"
            Path.Combine [|outputDir;fileName|]
        logger.Trace (sprintf "outFilePath:%s" outFilePath)
        result.GainedPeptides
        |> SeqIO'.csv "\t" true false
        |> FSharpAux.IO.SeqIO.Seq.writeOrAppend (outFilePath)
        logger.Trace "Writing Results:finished"
                
    open CLIArgumentParsing
    open Argu 
    
    let execute argv =
        let errorHandler = ProcessExiter(colorizer = function ErrorCode.HelpText -> None | _ -> Some System.ConsoleColor.Red)
        let parser = ArgumentParser.Create<CLIArguments>(programName =  (System.Reflection.Assembly.GetExecutingAssembly().GetName().Name),errorHandler=errorHandler)     
        let directory = Environment.CurrentDirectory
        let getPathRelativeToDir = getRelativePath directory
        let results = parser.Parse argv
        let i   = results.GetResult TargetFiles |> List.map getPathRelativeToDir
        let ii  = results.GetResult SourceFiles |> List.map getPathRelativeToDir 
        let o   = results.GetResult OutputDirectory    |> getPathRelativeToDir
        // let p   = results.GetResult ParamFile          |> getPathRelativeToDir
        let dc  = results.Contains DiagnosticCharts
        Logging.generateConfig o
        let logger = Logging.createLogger "QuantBasedAlignment"
        logger.Info (sprintf "InputFilePath -i = %A" i)
        logger.Info (sprintf "OutputFilePath -o = %s" o)
        // logger.Info (sprintf "ParamFilePath -p = %s" p)
        logger.Trace (sprintf "CLIArguments: %A" results)
        Directory.CreateDirectory(o) |> ignore
        //let p =
        //    Json.ReadAndDeserialize<Dto.QuantificationParams> p
        //    |> Dto.QuantificationParams.toDomain
        let targetFiles = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.quant"))) i
            |> Array.ofSeq
        let sourceFiles = 
            parsePaths (fun path -> Directory.GetFiles(path,("*.quant"))) ii
            |> Array.ofSeq
        if targetFiles.Length = 1 then
            logger.Info "single file detected"
            alignFiles dc {Placeholder=true} o sourceFiles targetFiles.[0]
        else
            logger.Info "directory found"
            logger.Trace (sprintf ".quant files : %A" targetFiles)
            let c =
                match results.TryGetResult Parallelism_Level with
                | Some c    -> c
                | None      -> 1
            logger.Trace (sprintf "Program is running on %i cores" c)
            targetFiles
            |> FSharpAux.PSeq.withDegreeOfParallelism c
            |> FSharpAux.PSeq.iter (fun (t) -> 
                let sourceFiles' = 
                    sourceFiles
                    |> Array.filter (fun x -> x <> t)
                alignFiles dc {Placeholder=true} o sourceFiles' t)
            |> ignore
        logger.Info "Done"




































