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
                // a blob of another layout (not ending after the three arrays, arrays of different lengths) is read by MzIO
                if position = data.Length && mzs.Length = intensities.Length && mzs.Length = ionMobilities.Length then
                    ofArrays mzs intensities ionMobilities
                else
                    fallback (reader.ReadSpectrumPeaks spectrumID)
            | _ ->
                fallback (reader.ReadSpectrumPeaks spectrumID)
        read, cmd

    /// Encodes the peak data of TIMS spectra: the intensities, the m/z and the ion mobilities, each as length and deflated bytes
    /// of the float64 values. These are the bytes MzIO's encoder writes for ZLib compression and a Peak1DArray with ion mobility
    /// (the layout initReadPeakArrays decodes), made without peak objects and with buffers that are reused from spectrum to spectrum.
    type PeakDataEncoder() =
        let output = new IO.MemoryStream(1 <<< 20)
        let lengthPlaceholder : byte[] = Array.zeroCreate 4
        let mutable bytes : byte[] = Array.zeroCreate (1 <<< 16)

        let write (values: float[]) (count: int) =
            let lengthPosition = int output.Position
            output.Write(lengthPlaceholder, 0, 4)
            if bytes.Length < count * 8 then bytes <- Array.zeroCreate (max (count * 8) (bytes.Length * 2))
            Buffer.BlockCopy(values, 0, bytes, 0, count * 8)
            do
                use deflate = new IO.Compression.DeflateStream(output, IO.Compression.CompressionMode.Compress, true)
                deflate.Write(bytes, 0, count * 8)
            let length = int output.Position - lengthPosition - 4
            let buffer = output.GetBuffer()
            buffer.[lengthPosition] <- byte length
            buffer.[lengthPosition + 1] <- byte (length >>> 8)
            buffer.[lengthPosition + 2] <- byte (length >>> 16)
            buffer.[lengthPosition + 3] <- byte (length >>> 24)

        /// The peak data of the first count values of the three arrays.
        member _.Encode(intensities: float[], mzs: float[], ionMobilities: float[], count: int) : byte[] =
            output.SetLength 0L
            write intensities count
            write mzs count
            write ionMobilities count
            output.ToArray()
