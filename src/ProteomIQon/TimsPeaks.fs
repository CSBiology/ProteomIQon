namespace ProteomIQon

open System
open System.Data.SQLite
open MzIO.Binary
open MzIO.Model

// Direct decoding of the peak data of TIMS spectra, shared by PSMBasedQuantificationTIMs and AlignmentBasedQuantification_v2.
module TimsPeaks =

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
                ofArrays mzs intensities ionMobilities
            | _ ->
                fallback (reader.ReadSpectrumPeaks spectrumID)
        read, cmd
