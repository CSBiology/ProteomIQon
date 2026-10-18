namespace ProteomIQon

open System
open System.Data.SQLite
open MzIO.Binary
open MzIO.Model

// Direct decoding of the peak data of TIMS spectra and the m/z sort of their peaks, shared by PSMBasedQuantificationTIMs and AlignmentBasedQuantification_v2.
module TimsPeaks =
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


    /// Reads the peaks of a spectrum. The layout of TIMS spectra (ZLib, m/z, intensity and ion mobility as float64) is decoded
    /// from the peak data directly, without peak objects, and given to ofArrays as m/z, intensity and ion mobility in file order;
    /// other layouts are read by MzIO and given to fallback. Returns the read function and its command, which the caller
    /// disposes when the file is done.
    let initReadPeakArrays (reader: MzIO.MzSQL.MzSQL) (fallback: Peak1DArray -> 'T) (ofArrays: float[] -> float[] -> float[] -> 'T) =
        let cmd = new SQLiteCommand("SELECT PeakArray, PeakData FROM Spectrum WHERE SpectrumID = @spectrumID", reader.Connection)
        let spectrumIDParameter = cmd.Parameters.Add("@spectrumID", Data.DbType.String)
        let inflate (data: byte[]) (offset: int) (length: int) =
            use input = new IO.MemoryStream(data, offset, length)
            use deflate = new IO.Compression.DeflateStream(input, IO.Compression.CompressionMode.Decompress)
            use output = new IO.MemoryStream()
            deflate.CopyTo output
            let values = Array.zeroCreate<float> (int output.Length / 8)
            Buffer.BlockCopy(output.GetBuffer(), 0, values, 0, values.Length * 8)
            values
        let read (spectrumID: string) =
            spectrumIDParameter.Value <- spectrumID
            // the last row of the spectrum, as MzIO
            let mutable row = None
            do
                use rows = cmd.ExecuteReader()
                while rows.Read() do
                    row <- Some (MzIO.Json.MzIOJson.FromJson<Peak1DArray>(rows.GetString 0), rows.GetValue 1 :?> byte[])
            match row with
            | Some (peakArray, data)
                when
                    peakArray.CompressionType = BinaryDataCompressionType.ZLib
                    && peakArray.IntensityDataType = BinaryDataType.Float64
                    && peakArray.MzDataType = BinaryDataType.Float64
                    && peakArray.IonMobilityDataType = Some BinaryDataType.Float64
                ->
                // intensities, m/z and ion mobilities, each as length and deflated bytes
                let mutable position = 0
                let next () =
                    let length = BitConverter.ToInt32(data, position)
                    let values = inflate data (position + 4) length
                    position <- position + 4 + length
                    values
                let intensities = next ()
                let mzs = next ()
                let ionMobilities = next ()
                // a blob of another layout (not ending after the three arrays, arrays of different lengths) is read by MzIO; a
                // truncated or corrupt blob throws
                if position = data.Length && mzs.Length = intensities.Length && mzs.Length = ionMobilities.Length then
                    ofArrays mzs intensities ionMobilities
                else
                    fallback (reader.ReadSpectrumPeaks spectrumID)
            | _ ->
                fallback (reader.ReadSpectrumPeaks spectrumID)
        read, cmd
