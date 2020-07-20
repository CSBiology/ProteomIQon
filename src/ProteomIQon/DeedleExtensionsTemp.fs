namespace Deedle

open Deedle
open Deedle.Vectors
open Deedle.Indices
open Deedle.Indices.Linear
open FSharpAux
//open FSharpAux.Collections

[<AutoOpen>]
module Frame =

    /// A generic transformation that works when at most one value is defined
    let private atMostOne = 
        { new IBinaryTransform with
            member vt.GetFunction<'R>() = (fun (l:OptionalValue<'R>) (r:OptionalValue<'R>) -> 
              if l.HasValue && r.HasValue then invalidOp "Combining vectors failed - both vectors have a value."
              if l.HasValue then l else r)
            member vt.IsMissingUnit = true }
        |> VectorListTransform.Binary

    let private transformColumn (vectorBuilder:IVectorBuilder) scheme rowCmd (vector:IVector) = 
      { new VectorCallSite<IVector> with
          override x.Invoke<'T>(col:IVector<'T>) = 
            vectorBuilder.Build<'T>(scheme, rowCmd, [| col |]) :> IVector }
      |> vector.Invoke

    /// Reorder elements in the index to match with another index ordering after applying the given function to the first index
    let private reindexBy (keyFunc : 'RowKey1 -> 'RowKey2) (index1:IIndex<'RowKey1>) (index2:IIndex<'RowKey2>) vector = 
        let relocations = 
                seq {  
                    for KeyValue(key, newAddress) in index1.Mappings do
                    let key = keyFunc key
                    let oldAddress = index2.Locate(key)
                    if oldAddress <> Addressing.Address.invalid then 
                        yield newAddress, oldAddress }
        Vectors.Relocate(vector, index1.KeyCount, relocations)

    ///Returns all possible combinations of second and first frame keys sharing the same keyfunction result
    let private combineIndexBy (keyFunc1 : 'RowKey1 -> 'SharedKey) (keyFunc2 : 'RowKey2 -> 'SharedKey) (index1:IIndex<'RowKey1>) (index2:IIndex<'RowKey2>) =
        // Group by the shared key
        let g1 = index1.Keys |> Seq.groupBy keyFunc1
        let g2 = index2.Keys |> Seq.groupBy keyFunc2
    
        //Intersect over the shared Keys
        let s1,s2 = g1 |> Seq.map fst |> Set.ofSeq, g2 |> Seq.map fst |> Set.ofSeq    
        let keyIntersect =  Set.intersect s1 s2 

        //For each shared key, create all possible combinations of index1 and index2 keys
        let m1,m2 = g1 |> Map.ofSeq,g2 |> Map.ofSeq
        let newKeys = 
            keyIntersect
            |> Seq.collect (fun key -> 
                    Map.find key m1
                    |> Seq.collect(fun k1 -> 
                        Map.find key m2
                        |> Seq.map (fun k2 -> 
                            key,k1,k2)
                        )
                )

        //Return new index
        LinearIndexBuilder.Instance.Create(newKeys,None)

    /// Create transformation on indices/vectors representing the align operation
    let private createAlignTransformation (keyFunc1 : 'RowKey1 -> 'SharedKey) (keyFunc2 : 'RowKey2 -> 'SharedKey) (thisIndex:IIndex<_>) (otherIndex:IIndex<_>) thisVector otherVector =
        let combinedIndex = combineIndexBy keyFunc1 keyFunc2 thisIndex otherIndex
        let rowCmd1 = reindexBy (fun (a,b,c) -> b) combinedIndex thisIndex thisVector
        let rowCmd2 = reindexBy (fun (a,b,c) -> c) combinedIndex otherIndex otherVector
        combinedIndex,rowCmd1,rowCmd2

    /// Create transformation on indices/vectors representing the expand operation
    let private createExpandTransformation (expandF : 'R -> 'RS seq) (index : IIndex<'R>) (vector : VectorConstruction) = 
        ///new Keys * old Keys 
        let newKeys = 
            index.Keys
            |> Seq.collect (fun r -> expandF r |> Seq.map (fun rs -> rs,r))
        let keyFunc = 
            let m = newKeys |> Map.ofSeq
            fun rs -> m.[rs]
        let newIndex = LinearIndexBuilder.Instance.Create(newKeys |> Seq.map fst,None)
        let rowCmd = reindexBy keyFunc newIndex index vector
        newIndex,rowCmd
        
    
    let distinctRowValues colName (df:Frame<int,_>) = 
        df
        |> Frame.groupRowsByString colName
        |> Frame.applyLevel fst (fun os -> os.FirstValue())
        |> Frame.indexRowsOrdinally


    // Appends a given series by a key it's value 
    let append (s:Series<'key,'value>) key value =
        s
        |> Series.observations
        |> Seq.append [(key,value)]
        |> Series.ofObservations



    let composeRowsBy f (keyColName:'a) (valueColName:'a) (table:Frame<_,'a>)=
        let keys =
            table.GetColumn<'a>(keyColName)
            |> Series.values
        let values = 
            table.GetColumn<'a>(valueColName)
            |> Series.values
    
        Seq.zip keys values
        |> Map.compose
        |> Map.map (fun k v -> f v )
        |> Map.toSeq
        |> Series.ofObservations


    let decomposeRowsBy f colName (df:Frame<int,_>) =
        df
        |> Frame.groupRowsByString colName
        |> Frame.rows
        |> Series.observations
        |> Seq.collect (fun (k,os) -> f k
                                      |> Seq.mapi (fun i v -> let k',i' = k
                                                              (k',i'+ i),append os "decomposed" (v))) // * (box v)
        |> Frame.ofRows
        |> Frame.indexRowsOrdinally
    
 
     // Frame toJaggedArray
    let Frame_ofJaggedArrayCol (rowNames:'rKey seq) (colNames:'cKey seq) (colJarray:'value array array) =       
        //TODO:
        // If length
    
        colJarray
        |> Seq.map2 (fun colKey arr -> colKey,Series(rowNames, arr)  ) colNames 
        |> frame 


    /// Align two data frames by a shared key received through mapping functions. 
    ///
    /// The columns of the joined frames must not overlap and their rows are aligned and multiplied
    /// by the shared key. The keyFuncs are used to map the rowKeys of the two frames to a shared key. 
    /// The resulting keys will result from the intersection of the shared Keys.
    ///
    /// The key of the resulting frame will be a triplet of shared key and the two input keys.
    let align (keyFunc1 : 'RowKey1 -> 'SharedKey) (keyFunc2 : 'RowKey2 -> 'SharedKey) (frame1 : Frame<'RowKey1, 'TColumnKey>) (frame2 : Frame<'RowKey2, 'TColumnKey>) =  
        //Get needed transformation objects and data form the Frame
        let index1 = frame1.RowIndex
        let index2 = frame2.RowIndex
        let indexBuilder = LinearIndexBuilder.Instance
        let vectorbuilder = ArrayVector.ArrayVectorBuilder.Instance
        let data1 = frame1.GetFrameData().Columns |> Seq.map snd |> ``F# Vector extensions``.Vector.ofValues
        let data2 = frame2.GetFrameData().Columns |> Seq.map snd |> ``F# Vector extensions``.Vector.ofValues
        // Intersect mapped row indices and get transformations to apply to input vectors
        let newRowIndex, rowCmd1, rowCmd2 = 
          createAlignTransformation keyFunc1 keyFunc2 index1 index2 (Vectors.Return 0) (Vectors.Return 0)
        // Append the column indices and get transformation to combine them
        let newColumnIndex, colCmd = 
            indexBuilder.Merge( [(frame1.ColumnIndex, Vectors.Return 0); (frame2.ColumnIndex, Vectors.Return 1) ], atMostOne, true)
        // Apply transformation to both data vectors
        let newData1 = data1.Select(transformColumn vectorbuilder newRowIndex.AddressingScheme rowCmd1)
        let newData2 = data2.Select(transformColumn vectorbuilder newRowIndex.AddressingScheme rowCmd2)
        // Combine column vectors a single vector & return results
        let newData = vectorbuilder.Build(newColumnIndex.AddressingScheme, colCmd, [| newData1; newData2 |])
        Frame(newRowIndex, newColumnIndex, newData, indexBuilder, vectorbuilder)

    /// Applies the function f to the rowKeys of the frame and adds the result as a new column to the frame
    let columnOfRowKeysBy (columnKey : 'C) (f : 'R -> 'T) (frame : Frame<'R,'C>) : Frame<'R,'C> =
        let newColumn = 
            frame.RowKeys
            |> Seq.map f
            |> Seq.zip frame.RowKeys
            |> series
        Frame.addCol columnKey newColumn frame

    /// Adds the rowKeys as a new column to the frame
    let columnOfRowKeys (columnKey : 'C) (frame : Frame<'R,'C>) : Frame<'R,'C> =
        columnOfRowKeysBy columnKey id frame
    
    /// If the predicate returns false for a value, replaces the value with missing
    let filter (predicate : 'R -> 'C -> 'a -> bool) (frame : Frame<'R,'C>) : Frame<'R,'C> =
        frame
        |> Frame.map (fun r c a -> 
            if predicate r c a then 
                OptionalValue a
            else 
                OptionalValue.Missing)

    /// If the predicate returns false for a value, replaces the value with missing
    let filterValues (predicate : 'a -> bool) (frame : Frame<'R,'C>) : Frame<'R,'C> =
        frame
        |> Frame.mapValues (fun a -> 
            if predicate a then 
                OptionalValue a
            else 
                OptionalValue.Missing)
    
    /// Creates a new data frame that contains only those columns of the original 
    /// data frame which contain at least one value.
    let dropEmptyRows (frame:Frame<'R, 'C>) = 
        //Get needed transformation objects and data form the Frame
        let indexBuilder =  LinearIndexBuilder.Instance
        let vectorbuilder = ArrayVector.ArrayVectorBuilder.Instance
        let data = frame.GetFrameData().Columns |> Seq.map snd |> ``F# Vector extensions``.Vector.ofValues

        // Create a combined vector that has 'true' for rows which have some values
        let hasSomeFlagVector = 
            Frame.rows frame
            |> Series.map (fun _ s -> Series.valuesAll s |> Seq.exists (fun opt -> opt.IsSome))
            |> Series.values
            |> ``F# Vector extensions``.Vector.ofValues

        // Collect all rows that have at least some values
        let newRowIndex, cmd = 
            indexBuilder.Search( (frame.RowIndex, Vectors.Return 0), hasSomeFlagVector, true)
        let newData = data.Select(transformColumn vectorbuilder newRowIndex.AddressingScheme cmd)
        Frame<_, _>(newRowIndex, frame.ColumnIndex, newData, indexBuilder, vectorbuilder)        

    /// Creates a new data frame that contains only those columns of the original 
    /// data frame which contain at least one value.
    let dropEmptyCols (frame:Frame<'R, 'C>) = 
        //Get needed transformation objects and data form the Frame
        let indexBuilder =  LinearIndexBuilder.Instance
        let vectorbuilder = ArrayVector.ArrayVectorBuilder.Instance
        let data = frame.GetFrameData().Columns |> Seq.map snd |> ``F# Vector extensions``.Vector.ofValues

        let newColKeys, newData =
            [| for KeyValue(colKey, addr) in frame.ColumnIndex.Mappings do
                match data.GetValue(addr) with
                | OptionalValue.Present(vec) when vec.ObjectSequence |> Seq.exists (fun opt -> opt.HasValue) ->
                    yield colKey, vec :> IVector
                | _ -> () |] |> Array.unzip
        let colIndex = indexBuilder.Create(Deedle.Internal.ReadOnlyCollection.ofArray newColKeys, None)
        Frame(frame.RowIndex, colIndex, vectorbuilder.Create(newData), indexBuilder, vectorbuilder )
            
    /// Creates a Frame where each row is mapped to multiple rows based on the input function.
    let expandRowsByKey (expandF : 'R -> 'RS seq) (frame : Frame<'R,'C>) : Frame<'RS,'C> =
        //Get needed transformation objects and data form the Frame
        let index = frame.RowIndex
        let indexBuilder = LinearIndexBuilder.Instance
        let vectorbuilder = ArrayVector.ArrayVectorBuilder.Instance
        let data = frame.GetFrameData().Columns |> Seq.map snd |> ``F# Vector extensions``.Vector.ofValues
        //expand rows via collection of function results
        let newRowIndex, rowCmd = 
          createExpandTransformation expandF index (Vectors.Return 0)

        // Apply transformation to both data vectors
        let newData = data.Select(transformColumn vectorbuilder newRowIndex.AddressingScheme rowCmd)
        // Combine column vectors a single vector & return results
        Frame(newRowIndex, frame.ColumnIndex, newData, indexBuilder, vectorbuilder)
 
    /// Creates a Frame where each row is mapped to multiple rows based on the input function. The input function takes the rowkey and the value of the given column at this rowkey and returns new keys.
    let expandRowsByColumn (column : 'C) (expandF : 'R -> 'V -> 'RS seq) (frame : Frame<'R,'C>) : Frame<'RS,'C> =
        //Get needed transformation objects and data form the Frame
        let index = frame.RowIndex
        let indexBuilder = LinearIndexBuilder.Instance
        let vectorbuilder = ArrayVector.ArrayVectorBuilder.Instance
        let data = frame.GetFrameData().Columns |> Seq.map snd |> ``F# Vector extensions``.Vector.ofValues
        /// The column over which the rows are expanded
        let column = frame.GetColumn column
        let expandF r = expandF r (column.Get r)
                   
        //expand rows via collection of function results
        let newRowIndex, rowCmd = 
          createExpandTransformation expandF index (Vectors.Return 0)

        // Apply transformation to both data vectors
        let newData = data.Select(transformColumn vectorbuilder newRowIndex.AddressingScheme rowCmd)
        // Combine column vectors a single vector & return results
        Frame(newRowIndex, frame.ColumnIndex, newData, indexBuilder, vectorbuilder)

