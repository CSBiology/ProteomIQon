namespace ProteomIQon

open System

// The peaks of a spectrum sorted by m/z, read by the XIC extraction of the quantification tools.
module SortedPeaks =
    /// The peaks of a spectrum sorted by m/z, ties in file order. Index is the position of the peak in the file,
    /// peaks without a m/z are left out. IonMobility is empty when the spectrum has none.
    type MzSortedPeaks =
        {
            Mz          : float[]
            Intensity   : float[]
            IonMobility : float[]
            Index       : int[]
        }

    /// Sorts the peaks of a spectrum, given as arrays in file order, by m/z. TIMS peaks are stored per mobility
    /// scan and are not sorted. The ion mobilities may be empty.
    let sortPeakArraysByMz (mzs: float[]) (intensities: float[]) (ionMobilities: float[]) =
        let order =
            Array.init mzs.Length id
            |> Array.filter (fun i -> not (Double.IsNaN mzs.[i]))
        let mz = order |> Array.map (fun i -> mzs.[i])
        Array.Sort(mz, order)
        // equal m/z in file order
        let mutable first = 0
        while first < mz.Length do
            let mutable last = first + 1
            while last < mz.Length && mz.[last] = mz.[first] do
                last <- last + 1
            if last - first > 1 then Array.Sort(order, first, last - first)
            first <- last
        {
            Mz          = mz
            Intensity   = order |> Array.map (fun i -> intensities.[i])
            IonMobility = if ionMobilities.Length > 0 then order |> Array.map (fun i -> ionMobilities.[i]) else [||]
            Index       = order
        }
