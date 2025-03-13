namespace ProteomIQon

open ProteomIQon
open FSharp.Stats
open Deedle
open FSharpAux

module RatioLFQ =

    let equalityConstrainedLS sumOfB (m:matrix) (b:vector) = 
        let mTm = (Matrix.scale 2. m.Transpose) * m  
        let leftM = 
            let tmp = Matrix.init (mTm.NumRows+1) (mTm.NumCols+1) (fun m n  -> 1.)
            tmp
            |> Matrix.mapi (fun m n v -> 
                if m < mTm.NumRows && n < mTm.NumCols then
                    mTm.[m,n]
                elif m = tmp.NumRows-1 && n = tmp.NumCols-1 then 
                    0.
                else 
                    v
                )
        let rightM = 
            let tmp = (m.Transpose |> Matrix.scale 2.) * b
            Vector.init (tmp.Length+1) (fun i -> if i < tmp.Length then tmp.[i] else sumOfB)
        FSharp.Stats.Algebra.LinearAlgebra.LeastSquares leftM rightM     


    let computeLFQProt (proteinName: string) (protMedianMap: Map<string,float>) allPepsOfProtein = 
        let loggedA = 
            allPepsOfProtein
        
        let loggedAFiltered =
            loggedA
            |> Array.choose id
        if loggedAFiltered.Length = 1 then None
        else
            let indexArray =
                loggedA
                |> Array.foldi (fun i acc x -> 
                    match x with
                    | Some _ -> i::acc
                    | None -> acc
                ) []
                |> List.rev
                |> Array.ofList
            // Matrix for peptides. Since proteins are only one row in matrix, can be reworked for lists. Should work either way.
            let ratioMatrix: Matrix<float> = 
                [|loggedAFiltered|]
                |> Matrix.ofJaggedArray
                |> FSharp.Stats.Correlation.Matrix.columnWiseCorrelationMatrix 
                    (fun peps1 peps2 -> 
                        Seq.zip peps1 peps2
                        |> Seq.map (fun (pep1,pep2) -> 
                            pep2 - pep1
                        )
                        |> Seq.median
                    )

            // let A,b :matrix*vector= 
            let numberOfComparisons = (FSharp.Stats.SpecialFunctions.Binomial.coeffcient (ratioMatrix.NumCols) 2) |> round 0 |> int
            let A = Matrix.create numberOfComparisons ratioMatrix.NumCols 0.
            let b = Vector.create numberOfComparisons 0.
            let mutable nB = 0
            for i = 1 to ratioMatrix.NumRows-1 do
                for j = 0 to i-1 do 
                    let rIJ = ratioMatrix.[i,j]
                    // hier
                    b.[nB] <- rIJ
                    A.[nB,j] <- -1.
                    A.[nB,i] <- 1.
                    nB <- nB + 1
            // don't take avg of ratios. Take median of medianOfRatios normalized N15 quants and median norm N15 quants.
            // Compare output matrix with input matrix. Should only deviate slightly.
            // average of runs where the protein was quantified multiplied by the number of runs
            let avgSumI =  
                match protMedianMap.TryFind proteinName with
                | Some name -> protMedianMap.[proteinName] * (float loggedAFiltered.Length)
                | None -> failwith proteinName

            let myLFQsRel = equalityConstrainedLS avgSumI A b
    
            let myLFQsRelExp = 
                myLFQsRel
                |> Vector.map (fun x -> 2.**x)
                |> Vector.toArray
                |> Array.take (ratioMatrix.NumCols)

            let finalRatioMatrix = 
                myLFQsRelExp
                |> Array.map (fun i ->  
                    myLFQsRelExp
                    |> Array.map (fun j -> i / j)
                    )
                |> Matrix.ofJaggedArray
        
            let myLFQsRelExpFilled =
                let arr = Array.zeroCreate loggedA.Length
                Array.map2 (fun i value ->
                    arr.[i] <- value
                ) indexArray myLFQsRelExp
                |> ignore
                arr
    
            (myLFQsRelExpFilled,finalRatioMatrix,avgSumI)
            |> Some

    let readQuantResult (path: string) (proteinColumn: string) = 
        let data: Frame<_,string> =
            Frame.ReadCsv(path, separators="\t")
            |> Frame.indexRowsUsing (fun s ->
                s.GetAs<string>(proteinColumn)
            )
        data

    let combinedRatio (path: string) (proteinColumn: string) (ratioColumnEnding: string) =
        readQuantResult path proteinColumn
        |> Frame.filterCols (fun ck os -> ck.EndsWith(ratioColumnEnding))
        |> Frame.fillMissingWith nan
        |> Frame.mapValues log2
        |> Frame.filterRows (fun rk s -> isNan (Stats.mean s) |> not)

    let combinedOriginal (path: string) (proteinColumn: string) (referenceColumnEnding: string) =
        readQuantResult path proteinColumn
        |> Frame.filterCols (fun ck os -> ck.EndsWith(referenceColumnEnding))
        |> Frame.fillMissingWith nan
        |> Frame.mapValues log2
        |> Frame.mapRows (fun rk s -> Stats.mean s)
        |> Series.observations
        |> Map.ofSeq

    let lfq (inputPath: string) (outputPath: string) (proteinColumn: string) (ratioColumnEnding: string) (referenceColumnEnding: string) =
        let combinedRatioRead = combinedRatio inputPath proteinColumn ratioColumnEnding
        let combinedOriginalRead = combinedOriginal inputPath proteinColumn referenceColumnEnding
        let ck = 
            combinedRatioRead.ColumnKeys
            |> Array.ofSeq
            |> Array.map (fun ck -> ck + "_LFQ")
        combinedRatioRead
        |> Frame.mapRows(fun k f -> 
            k,
            f 
            |> Series.valuesAll
            |> Array.ofSeq
            |> Array.map (fun x -> 
                match x with
                | Some value -> Some (value :?> float)
                | None -> None
            )
            |> computeLFQProt k combinedOriginalRead
        )
        |> Series.values
        |> Array.ofSeq
        |> Array.map (fun (prot,result) ->
            match result with
            | Some (lfq,ratioMatrix,scalingFactor) -> Some (prot,Array.zip ck lfq |> Series.ofObservations)
            | None -> None
        )
        |> Array.choose id
        |> Frame.ofRows
        |> fun f -> f.SaveCsv(outputPath, true, separator = '\t')