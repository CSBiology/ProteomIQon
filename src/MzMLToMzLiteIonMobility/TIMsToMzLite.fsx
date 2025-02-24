#r "nuget: ProteomIQon, 0.0.8"
#r "nuget: MzIO, 0.1.1"
#r "nuget: MzIO.Processing, 0.1.2"
#r "nuget: MzIO.SQL, 0.1.4"
#r "nuget: MzIO.MzML, 0.1.6"
#r "nuget: System.Data.SQLite.Core, 1.0.118"

open ProteomIQon
open ProteomIQon.Domain
open ProteomIQon.Core
open System.IO
open System.Collections.Generic
open BioFSharp.Mz
open MzIO
open MzIO.Binary
open MzIO.Model
open MzIO.Processing
open MzIO.IO
open MzIO.MzSQL
open MzIO.IO.MzML
open ProteomIQon.Core.MzIO.Processing
open ProteomIQon.Dto
open ProteomIQon.Domain
open MzIO.MetaData
open MzIO.MetaData.ParamEditExtension
open MzIO.MetaData.PSIMSExtension
open MzIO.MetaData.UO
open MzIO.MetaData.UO.UO
open MzIO.Commons.Arrays



open System
open System.Runtime.InteropServices
open System.Runtime.Serialization
open System.Text
open System.IO
open System.Data.SQLite

[<Struct>]
type ChromatogramJob =
    val mutable id : int64
    val mutable time_begin : float
    val mutable time_end : float
    val mutable mz_min : float
    val mutable mz_max : float
    val mutable ook0_min : float
    val mutable ook0_max : float

[<Struct>]
type MSMS_SPECTRUM_FUNCTOR_RESULT =
    {
        mutable precursorId: int64
        mutable numPeaks: uint32
        mutable mzValues: float array
        mutable areaValues: float array
    }

[<UnmanagedFunctionPointer(CallingConvention.Cdecl)>]
type MSMS_SPECTRUM_FUNCTOR = delegate of MSMS_SPECTRUM_FUNCTOR_RESULT -> unit

[<Struct>]
type MSMS_PROFILE_SPECTRUM_FUNCTOR_RESULT =
    {
        mutable precursorId: int64
        mutable numPoints: uint32
        mutable intensityValues: float array
        }

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint64 tims_open_v2(string analysis_directory, uint32 use_recalibrated_state, uint32 pressure_compensation_strategy)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern void tims_close(uint64 handle)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_get_last_error_string(StringBuilder error_string, uint32 error_string_length)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_has_recalibrated_state(uint64 handle)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_read_scans_v2(uint64 handle, int64 frame_id, uint32 scan_begin, uint32 scan_end, IntPtr buffer, uint32 buffer_length)

//[<UnmanagedFunctionPointer(CallingConvention.Cdecl)>]
//type MSMS_SPECTRUM_FUNCTOR = delegate of int64 * uint32 * float[] * float[] -> unit

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_read_pasef_msms(uint64 handle, IntPtr precursor_list, uint32 num_precursors, MSMS_SPECTRUM_FUNCTOR callback)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_read_pasef_msms_for_frame(uint64 handle, int64 frame_id, MSMS_SPECTRUM_FUNCTOR callback)

[<UnmanagedFunctionPointer(CallingConvention.Cdecl)>]
type MSMS_PROFILE_SPECTRUM_FUNCTOR = delegate of MSMS_PROFILE_SPECTRUM_FUNCTOR_RESULT -> unit

//[<UnmanagedFunctionPointer(CallingConvention.Cdecl)>]
//type MSMS_PROFILE_SPECTRUM_FUNCTOR = delegate of int64 * uint32 * int32[] -> unit
    
[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_read_pasef_profile_msms(uint64 handle, IntPtr precursor_list, uint32 num_precursors, MSMS_PROFILE_SPECTRUM_FUNCTOR callback)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_read_pasef_profile_msms_for_frame(uint64 handle, int64 frame_id, MSMS_PROFILE_SPECTRUM_FUNCTOR callback)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_extract_centroided_spectrum_for_frame_v2(uint64 handle, int64 frame_id, uint32 scan_begin, uint32 scan_end, MSMS_SPECTRUM_FUNCTOR callback, IntPtr context)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_extract_centroided_spectrum_for_frame_ext(uint64 handle, int64 frame_id, uint32 scan_begin, uint32 scan_end, double peak_picker_resolution, MSMS_SPECTRUM_FUNCTOR callback, IntPtr context)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_extract_profile_for_frame(uint64 handle, int64 frame_id, uint32 scan_begin, uint32 scan_end, MSMS_PROFILE_SPECTRUM_FUNCTOR callback, IntPtr context)

[<UnmanagedFunctionPointer(CallingConvention.Cdecl)>]
type CHROMATOGRAM_JOB_GENERATOR = delegate of ChromatogramJob[] * IntPtr -> uint32

[<UnmanagedFunctionPointer(CallingConvention.Cdecl)>]
type CHROMATOGRAM_TRACE_SINK = delegate of int64 * uint32 * int64[] * uint64[] * IntPtr -> uint32

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_extract_chromatograms(uint64 handle, CHROMATOGRAM_JOB_GENERATOR job_generator, CHROMATOGRAM_TRACE_SINK trace_sink, IntPtr user_data)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_index_to_mz(uint64 handle, int64 frame_id, IntPtr indices, IntPtr mzs, uint32 count)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_mz_to_index(uint64 handle, int64 frame_id, IntPtr mzs, IntPtr indices, uint32 count)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_scannum_to_oneoverk0(uint64 handle, int64 frame_id, IntPtr scan_nums, IntPtr mobilities, uint32 count)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_oneoverk0_to_scannum(uint64 handle, int64 frame_id, IntPtr mobilities, IntPtr scan_nums, uint32 count)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_scannum_to_voltage(uint64 handle, int64 frame_id, IntPtr scan_nums, IntPtr voltages, uint32 count)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern uint32 tims_voltage_to_scannum(uint64 handle, int64 frame_id, IntPtr voltages, IntPtr scan_nums, uint32 count)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double tims_oneoverk0_to_ccs_for_mz(double ook0, int32 charge, double mz)

[<DllImport(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\timsdata\win64\timsdata.dll", CallingConvention = CallingConvention.Cdecl)>]
extern double tims_ccs_to_oneoverk0_for_mz(double ccs, int32 charge, double mz)

//let private _throwLastTimsDataError (dllHandle: System.IntPtr) =
//    let len = tims_get_last_error_string(null, 0u)
//    let buf = Marshal.AllocHGlobal(int32 len)
//    tims_get_last_error_string(buf, uint32 len)
//    let errorMsg = Marshal.PtrToStringAnsi(buf)
//    Marshal.FreeHGlobal(buf)
//    raise (new System.Runtime.InteropServices.ExternalException(errorMsg))


type PressureCompensationStrategy =
    | NoPressureCompensation = 0
    | AnalyisGlobalPressureCompensation = 1
    | PerFramePressureCompensation = 2

// Define a helper function to allocate unmanaged memory and copy data to it
let copyToUnmanagedArrayFloat (data: float[]) =
    let ptr = Marshal.AllocHGlobal(data.Length * 8) // 8 bytes for a float
    Marshal.Copy(data,0,ptr,data.Length)
    // for i = 0 to data.Length - 1 do
    //     Marshal.WriteInt64(ptr, i * 8, BitConverter.DoubleToInt64Bits(data.[i]))
    ptr

// Define a helper function to allocate unmanaged memory and copy data to it
let copyToUnmanagedArrayInt64 (data: int64[]) =
    let ptr = Marshal.AllocHGlobal(data.Length * 8) // 8 bytes for a float
    Marshal.Copy(data,0,ptr,data.Length)
    // for i = 0 to data.Length - 1 do
    //     Marshal.WriteInt64(ptr, i * 8, BitConverter.DoubleToInt64Bits(data.[i]))
    ptr

// Define a helper function to read a float from unmanaged memory
let readFloat (ptr: IntPtr, index: int) =
    BitConverter.Int64BitsToDouble(Marshal.ReadInt64(ptr, index * 8))

// Helper function to read a uint32 from unmanaged memory
let readUInt32 (ptr: IntPtr, index: int) =
    uint32(Marshal.ReadInt32(ptr, index * 4)) // 4 bytes for a int32

let executeSQLiteQueryFloat(connection: SQLiteConnection, commandText: string) =
    try
        use command = new SQLiteCommand(commandText, connection)
        let reader = command.ExecuteReader()
        let rec loop (reader:SQLiteDataReader) acc =
            match reader.Read() with
            | true  -> loop reader ((reader.GetFloat(0))::acc)
            | false -> acc |> List.rev
        loop reader []
    with
    | :? SQLiteException as ex ->
        // Handle exceptions specific to SQLite
        failwithf "SQLite Exception: %s" ex.Message
    | ex ->
        // Handle other exceptions
        failwithf "Exception: %s" ex.Message

type TimsData(analysisDirectory: string, ?useRecalibratedState: bool, ?pressureCompensationStrategy: PressureCompensationStrategy) =

    let mutable handle = 0UL
    let mutable conn : SQLiteConnection = null
    let mutable initialFrameBufferSize = 128u

    do
        let useRecalibrated = defaultArg useRecalibratedState false
        let pressureStrategy = defaultArg pressureCompensationStrategy PressureCompensationStrategy.NoPressureCompensation

        handle <- tims_open_v2(analysisDirectory, (if useRecalibrated then 1u else 0u), uint32 pressureStrategy)

        //if handle = 0UL then
        //    _throwLastTimsDataError()

        conn <- new SQLiteConnection(Path.Combine("Data source="+analysisDirectory, "analysis.tdf"))
        conn.Open()

    member this.Handle
        with get() = handle
        and set(value) = handle <- value

    member this.InitialFrameBufferSize
        with get() = initialFrameBufferSize
        and set(value) = initialFrameBufferSize <- value

    member this.Conn
        with get() = conn
        and set(value) = conn <- value

    member this.Close() =
        if handle <> 0UL then
            tims_close(handle)
            handle <- 0UL
        if not (conn = null) then
            conn.Close()
            conn <- null

    member this.CallConversionFunc(frameId: int64, inputArray: float[], func: (uint64 * int64 * IntPtr * IntPtr * uint32) -> uint32) =
        let inArray = copyToUnmanagedArrayFloat inputArray
        let cnt = inputArray.Length
        let outArray = Marshal.AllocHGlobal(int32 (8 * int cnt)) // 8 bytes for a float

        try
            let success =
                func (this.Handle, frameId, inArray, outArray, (uint32 cnt))
            //if success = 0u then
            //    _throwLastTimsDataError this.dll
            let result = Array.init cnt (fun i -> readFloat(outArray, i))
            result
        finally
            Marshal.FreeHGlobal inArray
            Marshal.FreeHGlobal outArray

    member this.IndexToMz (frameId: int64, indices: float array) =
        this.CallConversionFunc(frameId, indices, tims_index_to_mz)
        
    member this.MzToIndex (frameId: int64, mzs: float array) =
        this.CallConversionFunc(frameId, mzs, tims_mz_to_index)

    member this.ScanNumToOneOverK0 (frameId: int64, scanNums: float array) =
        this.CallConversionFunc(frameId, scanNums, tims_scannum_to_oneoverk0)

    member this.OneOverK0ToScanNum (frameId: int64, mobilities: float array) =
        this.CallConversionFunc(frameId, mobilities, tims_oneoverk0_to_scannum)

    member this.ScanNumToVoltage (frameId: int64, scanNums: float array) =
        this.CallConversionFunc(frameId, scanNums, tims_scannum_to_voltage)

    member this.VoltageToScanNum (frameId: int64, voltages: float array) =
        this.CallConversionFunc(frameId, voltages, tims_voltage_to_scannum)


    interface IDisposable with
        member this.Dispose() =
            this.Close()
            
    member this.ReadScansDllBuffer (frameId: int64, scanBegin: uint32, scanEnd: uint32) =
        // Buffer-growing loop
        let rec growBuffer cnt =
            let buf = Marshal.AllocHGlobal(int cnt * 4) // Create an empty array
            let len = cnt * 4u // 4 bytes for a uint32
            let requiredLen = tims_read_scans_v2(this.Handle, frameId, scanBegin, scanEnd, buf, len)
            if requiredLen = 0u then
                failwith "An error occurred."
            if requiredLen > uint32 len then
                growBuffer (requiredLen / 4u + 1u)
            else
                let result = Array.init (int cnt) (fun i -> readUInt32(buf, i))
                Marshal.FreeHGlobal buf // Release unmanaged memory
                result
        growBuffer this.InitialFrameBufferSize

    member this.ReadScans (frameId: int64, scanBegin: uint32, scanEnd: uint32): (uint32[]*uint32[])[]  =
        let buf = this.ReadScansDllBuffer(frameId, scanBegin, scanEnd)
        //printfn "here"
        let mutable result = []
        let mutable d = int (scanEnd - scanBegin)
        for i in [(int scanBegin) .. (int scanEnd) - 1] do
            let nPeaks = buf.[i - int scanBegin]
            let indices = Array.sub buf d (int nPeaks)
            d <- d + int nPeaks
            let intensities = Array.sub buf d (int nPeaks)
            d <- d + int nPeaks
            result <- (indices, intensities) :: result

        result |> List.rev
        |> List.toArray

    
    //member this.ReadPasefMsMs (precursorList: int64 array) =
    //    let precursorsForDll = copyToUnmanagedArrayInt64 precursorList
    //    let result = new System.Collections.Generic.Dictionary<int64, (float array * float array)>()

    //    let callbackForDll(resultStruct: MSMS_SPECTRUM_FUNCTOR_RESULT) =
    //        result.Add(resultStruct.precursorId, (resultStruct.mzValues, resultStruct.areaValues))

    //    let rc = tims_read_pasef_msms(handle, precursorsForDll, uint32 precursorList.Length, callbackForDll)

    //    result

    //member this.ReadPasefMsMsForFrame (frameId: int64) =
    //    let result = new System.Collections.Generic.Dictionary<int64, (float array * float array)>()

    //    let callbackForDll(resultStruct: MSMS_SPECTRUM_FUNCTOR_RESULT) =
    //        result.Add(resultStruct.precursorId, (resultStruct.mzValues, resultStruct.areaValues))

    //    let rc = tims_read_pasef_msms_for_frame(handle, frameId, callbackForDll)

    //    result

    //member this.ReadPasefProfileMsMs (precursorList: int64 array) =
    //    let precursorsForDll = copyToUnmanagedArrayInt64 precursorList
    //    let result = new System.Collections.Generic.Dictionary<int64, float array>()

    //    let callbackForDll(resultStruct: MSMS_PROFILE_SPECTRUM_FUNCTOR_RESULT) =
    //        result.Add(resultStruct.precursorId, resultStruct.intensityValues)

    //    let rc = tims_read_pasef_profile_msms(handle, precursorsForDll, uint32 precursorList.Length, callbackForDll)

    //    result

    //member this.ReadPasefProfileMsMsForFrame (frameId: int64) =
    //    let result = new System.Collections.Generic.Dictionary<int64, float array>()

    //    let callbackForDll(resultStruct: MSMS_PROFILE_SPECTRUM_FUNCTOR_RESULT) =
    //        result.Add(resultStruct.precursorId, resultStruct.intensityValues)

    //    let rc = tims_read_pasef_profile_msms_for_frame(handle, frameId, callbackForDll)

    //    result

    //member this.ExtractCentroidedSpectrumForFrame (frameId: int64, scanBegin: uint32, scanEnd: uint32, peakPickerResolution: double option) =
    //    let mutable result: MSMS_SPECTRUM_FUNCTOR_RESULT = 
    //        {
    //            precursorId = 0L
    //            numPeaks = 0u
    //            mzValues =[||]
    //            areaValues = [||]
    //        }

    //    let callbackForDll(resultStruct: MSMS_SPECTRUM_FUNCTOR_RESULT) =
    //        result.mzValues <- resultStruct.mzValues
    //        result.areaValues <-  resultStruct.areaValues

    //    let rc =
    //        match peakPickerResolution with
    //        | Some resolution ->
    //            tims_extract_centroided_spectrum_for_frame_ext(handle, frameId, scanBegin, scanEnd, resolution, callbackForDll, IntPtr.Zero)
    //        | None ->
    //            tims_extract_centroided_spectrum_for_frame_v2(handle, frameId, scanBegin, scanEnd, callbackForDll, IntPtr.Zero)

    //    result

    //member this.ExtractProfileForFrame (frameId: int64, scanBegin: uint32, scanEnd: uint32) =
    //    let mutable result =
    //        {
    //            precursorId = 0L
    //            numPoints = 0u
    //            intensityValues = [||]
    //        }

    //    let callbackForDll(resultStruct: MSMS_PROFILE_SPECTRUM_FUNCTOR_RESULT) =
    //        result.intensityValues <- resultStruct.intensityValues

    //    let rc = tims_extract_profile_for_frame(handle, frameId, scanBegin, scanEnd, callbackForDll, IntPtr.Zero)

    //    result


    ///
let private initPeakPicking (*(logger: NLog.Logger)*) (peakPickingParams:PeakPicking) (peaks: Peak1DArray) =

    match peakPickingParams with
    | PeakPicking.ProfilePeaks ->
        peaks.Peaks
        |> Core.MzIO.Peaks.unzipIMzliteArray
    | PeakPicking.Centroid (CentroidizationMode.Wavelet waveletParams) ->
        match waveletParams.PaddingParams with
        | Some pParams ->
            let initPaddingParameters yThreshold =
                SignalDetection.Padding.createPaddingParameters
                    yThreshold
                    pParams.MaximumPaddingPoints
                    pParams.Padding_MzTolerance
                    pParams.WindowSize
                    pParams.SpacingPerc

            let initwaveletParameters yThreshold =
                SignalDetection.Wavelet.createWaveletParameters
                    waveletParams.NumberOfScales
                    yThreshold
                    waveletParams.Centroid_MzTolerance
                    waveletParams.SNRS_Percentile
                    waveletParams.MinSNR
                    waveletParams.RefineMZ
                    waveletParams.SumIntensities
            match waveletParams.YThreshold with
            | YThreshold.Fixed yThreshold ->
                let paddingParams = initPaddingParameters yThreshold
                let waveletParameters = initwaveletParameters yThreshold
                let mzData, intensityData =
                    peaks.Peaks
                    |> Core.MzIO.Peaks.unzipIMzliteArray
                let paddedMz,paddedIntensity =
                    SignalDetection.Padding.paddDataBy paddingParams mzData intensityData
                BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters paddedMz paddedIntensity
            | YThreshold.MinSpectrumIntensity ->
                //logger.Trace (sprintf "ms2Centroidization with: %A" waveletParams)
                let mzData, intensityData =
                    peaks.Peaks
                    |> Core.MzIO.Peaks.unzipIMzliteArray
                let yThreshold = Array.min intensityData
                let paddedMz,paddedIntensity =
                    let paddingParams = initPaddingParameters yThreshold
                    SignalDetection.Padding.paddDataBy paddingParams mzData intensityData
                let waveletParameters = initwaveletParameters yThreshold
                BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters paddedMz paddedIntensity
        | None ->
            let initwaveletParameters yThreshold =
                SignalDetection.Wavelet.createWaveletParameters
                    waveletParams.NumberOfScales
                    yThreshold
                    waveletParams.Centroid_MzTolerance
                    waveletParams.SNRS_Percentile
                    waveletParams.MinSNR
                    waveletParams.RefineMZ
                    waveletParams.SumIntensities
            match waveletParams.YThreshold with
            | YThreshold.Fixed yThreshold ->
                //logger.Trace (sprintf "ms1Centroidization with: %A" waveletParams)
                let waveletParameters = initwaveletParameters yThreshold
                let mzData, intensityData =
                    peaks.Peaks
                    |> Core.MzIO.Peaks.unzipIMzliteArray
                BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters mzData intensityData
            | YThreshold.MinSpectrumIntensity ->
                let mzData, intensityData =
                    peaks.Peaks
                    |> Core.MzIO.Peaks.unzipIMzliteArray
                let yThreshold = Array.min intensityData
                let waveletParameters = initwaveletParameters yThreshold
                BioFSharp.Mz.SignalDetection.Wavelet.toCentroidWithRicker2D waveletParameters mzData intensityData
    | PeakPicking.Centroid CentroidizationMode.Manufacturer ->
        failwith "Manufacturer peak picking is only supported for .baf (Bruker) files."
    
let binBy (projection: 'a -> float) bandwidth (data: seq<'a>) =
    if bandwidth = 0. then raise (System.DivideByZeroException("Bandwidth cannot be 0."))
    let halfBw = bandwidth / 2.0
    let decBandwidth = decimal bandwidth
    let tmp = 
        data
        |> Seq.groupBy (fun x -> (decimal (projection x) / decBandwidth) |> float |> floor) 
        |> Seq.map (fun (k,values) -> 
            let count = (Seq.length(values)) |> float
            if k < 0. then
                ((k  * bandwidth) + halfBw, values)   
            else
                ((k + 1.) * bandwidth) - halfBw, values)
        |> Seq.sortBy fst
    tmp

let fixSpectrum (m:MzIO.Model.MassSpectrum) =
    if isNull(m.Precursors) then
        m.Precursors <- new MzIO.Model.PrecursorList()
    if isNull(m.Scans) then
        m.Scans <- new MzIO.Model.ScanList()
    if isNull(m.Products) then
        m.Products <- new MzIO.Model.ProductList()
    m

let createPeak1DArrayCopy (source: MzIO.Binary.Peak1DArray) =
    let pa = MzIO.Binary.Peak1DArray()
    pa.CompressionType <- source.CompressionType
    pa.IntensityDataType <- source.IntensityDataType
    pa.MzDataType <- source.MzDataType
    pa

let createBinnedPeaks (binSize: float) (mz: float[]) (intensity: float[]) (ionMobility: float[]) = 

    let zippedPeaks = (Array.zip mz intensity) |> Seq.zip (ionMobility)

    let binnedPeakData =
        zippedPeaks
        |> binBy (fun (mirim, peak) -> mirim) binSize

    binnedPeakData
    |> Seq.map(fun (bin, binnedData) ->
        let pa = MzIO.Binary.Peak1DArray()
        pa.CompressionType <- BinaryDataCompressionType.NoCompression
        pa.IntensityDataType <- BinaryDataType.Float64
        pa.MzDataType <- BinaryDataType.Float64
        pa.Peaks <-
            MzIO.Commons.Arrays.ArrayWrapper(
                binnedData
                |> Seq.map (fun (_,(mz,intensity)) -> new Peak1D (mz, intensity))
                |> Seq.toArray
            )
        bin,pa
    )

let insertSpectrum (compress:BinaryDataCompressionType) (outReader: MzSQL.MzSQL) (runID:string)
    (ms1PeakPicking: Peak1DArray -> float [] * float []) (ms2PeakPicking: Peak1DArray -> float [] * float [])
        (spectrum: MassSpectrum) (peaks: Peak1DArray) =
    match MassSpectrum.getMsLevel spectrum with
    | 1 ->
        let mzData,intensityData =
            try
            ms1PeakPicking peaks
            with
            | _ -> [||],[||]
        if Array.isEmpty mzData || Array.isEmpty intensityData then ()
        else
            let peaks' = PeakArray.createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData 
            outReader.Insert (runID, spectrum, peaks')
    | 2 ->
        let mzData,intensityData =
            try
            ms2PeakPicking peaks
            with
            | _ -> [||],[||]
        if Array.isEmpty mzData || Array.isEmpty intensityData then ()
        else
            let peaks' = PeakArray.createPeak1DArray compress BinaryDataType.Float64 BinaryDataType.Float64 mzData intensityData
            outReader.Insert (runID, spectrum, peaks')
    | _ ->
        failwith "Only mass spectra of level 1 and 2 are supported."

/// Returns the default runID used by manufacturers
let getDefaultRunID (mzReader:TimsData) = 
    match mzReader with
    | :? TimsData as r                 -> "sample=0"
    //| :? MzML.MzMLReader as r       -> "sample=0"
    //| :? MzMLReaderMIRIM as r       -> "sample=0"

let executeSQLiteQueryInt(connection: SQLiteConnection, commandText: string) =
    try
        use command = new SQLiteCommand(commandText, connection)
        let reader = command.ExecuteReader()
        let rec loop (reader:SQLiteDataReader) acc =
            match reader.Read() with
            | true  -> loop reader ((reader.GetInt32(0))::acc)
            | false -> acc |> List.rev
        loop reader []
    with
    | :? SQLiteException as ex ->
        // Handle exceptions specific to SQLite
        failwithf "SQLite Exception: %s" ex.Message
    | ex ->
        // Handle other exceptions
        failwithf "Exception: %s" ex.Message

let scanToPeak1DArray (scans:(uint32[]*uint32[])[]) (td: TimsData) frameID =
    let mutable i = 0.
    scans
    |> Array.collect (fun (indices,intensities) ->
        let mz = td.IndexToMz(int64 frameID, (indices |> Array.map float))
        let voltage = 
            td.ScanNumToOneOverK0(int64 frameID, [|i|])
            |> Array.head
            |> fun x -> Array.init mz.Length (fun _ -> x)
        i <- i + 1.
        Array.zip3 mz (intensities|> Array.map float) voltage
    )
    |> Array.unzip3

let processFile (processParams:MzMLtoMzLiteParams) (outputDir:string) (instrumentOutput:string) =

    //let logger = Logging.createLogger (Path.GetFileNameWithoutExtension instrumentOutput)

    //logger.Trace (sprintf "Input file: %s" instrumentOutput)
    //logger.Trace (sprintf "Output directory: %s" outputDir)
    //logger.Trace (sprintf "Parameters: %A" processParams)

    
    //let tmp = File.ReadAllText instrumentOutput
    //File.WriteAllText(instrumentOutput, tmp.Replace("&quot;", ""))

    let inReaderMS = new TimsData(instrumentOutput)

    let inReaderMS' = new MzMLReaderMIRIM(@"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\Ara_60min_wTrap_Aurora_DDA_Slot1-4_113.mzML")
    let inRunID  = "sample=0"
    let inTrMS = inReaderMS'.BeginTransaction()
    let spectra = inReaderMS'.ReadMassSpectra(inRunID)

    inReaderMS'.ResetReader()

    let ms1PeakPicking = initPeakPicking (*logger*) processParams.MS1PeakPicking
    let ms2PeakPicking = initPeakPicking (*logger*) processParams.MS2PeakPicking

    let outDirPath =
        let fileName = Path.GetFileNameWithoutExtension instrumentOutput
        Path.Combine(outputDir, fileName)
            
    //logger.Trace $"Creating directory {outDirPath} for binned results of {instrumentOutput}"
        
    Directory.CreateDirectory outDirPath |> ignore

    //logger.Trace $"Reading spectra from {instrumentOutput}"



    //logger.Trace "Done reading spectra"
    //logger.Trace $"Reading model from {instrumentOutput}"
        
    //let model = inReaderMS.Model
    //inReaderMS.ResetReader()
        
    //logger.Trace "Done reading model"

    //logger.Trace $"Start writing binned mzlite files"
    //logger.Trace $"Total number of binned files: {spectrumMap.Count}"
    let connectionMap = new Dictionary<string, MzSQL.MzSQL*System.Data.SQLite.SQLiteTransaction>()
    spectra
    |> Seq.iteri (fun i (spectrum) ->
        if i % 1000 = 0 then
            printfn "Processing spectrum %i" (i+1)
        let mz,intensities,voltage =
            let scans =
                let scanCount = 
                    executeSQLiteQueryFloat(inReaderMS.Conn, $"SELECT NumScans FROM Frames WHERE Id={i+1}")
                    |> List.head
                inReaderMS.ReadScans(int64 (i+1), 0u, uint32 scanCount)
            scanToPeak1DArray scans inReaderMS (i+1)
        let binResult = createBinnedPeaks 0.002 mz intensities voltage
        binResult
        |> Seq.iter(fun (bin, peaks) ->
            printfn "%A" spectrum
            let outFile = Path.Combine(outputDir, $"binned_spectra_%.3f{bin}.mzlite")
            let outReader,outTr =
                if connectionMap.ContainsKey(outFile) then
                    connectionMap.[outFile]
                else
                    let outReader = new MzSQL(outFile)
                    let _ = outReader.Open()
                    let outTr = outReader.BeginTransaction()
                    //try
                    //    outReader.InsertModel model
                    //    //logger.Trace "Model inserted."
                    //with
                    //    | ex -> failwith $"Inserting model failed: {ex}"
                    connectionMap.Add(outFile, (outReader, outTr))
                    outReader, outTr
            let outRunID  = "sample=0"
            insertSpectrum processParams.Compress outReader outRunID ms1PeakPicking ms2PeakPicking (spectrum) peaks
        )
    )
    for x in connectionMap do
        let outReader = fst x.Value
        let outTr = snd x.Value
        outTr.Commit()
        outTr.Dispose()
        outReader.Dispose()

    //logger.Trace "Done."

let deserialized = 
    System.IO.File.ReadAllText(@"C:\Users\jonat\source\repos\ProteomIQon\src\ProteomIQon\defaultParams\TIMsMzMLtoMzLiteParams.json")
    |> Json.deserialize<Dto.MzMLtoMzLiteParams>
    |> PreprocessingParams.toDomain
#time
//processFile deserialized @"C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\outTest" "C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\Ara_60min_wTrap_Aurora_DDA_Slot1-4_113.d"



let reader = new TimsData("C:\Users\jonat\OneDrive\Doktor\TIMsDataWrapper\Ara_60min_wTrap_Aurora_DDA_Slot1-4_113.d")

[|1 .. 10000|]
|> Array.map (fun i ->
    if i % 1000 = 0 then
        printfn "Processing spectrum %i" (i)
    //let scanCount = 
    //    //executeSQLiteQueryFloat(reader.Conn, $"SELECT NumScans FROM Frames WHERE Id={i}")
    //    |> List.head
    reader.ReadScans(int64 (i), 0u, uint32 944)

)
1
//let inTrMS = reader.BeginTransaction()
//let spectra = reader.ReadMassSpectra("sample=0") |> Array.ofSeq

//spectra.[0].ID
//reader.ResetReader()
//let data = reader.getSpecificPeak1DArraySequentialWithMIRIM(spectra.[50].ID)
//data.Peaks.Length
//(data?Mirim :?> float array).Length

//let binResult = createBinnedPeaks false 0.002 data

//let mz,intensity = data.Peaks |> Core.MzIO.Peaks.unzipIMzliteArray
//let ionMobility = (data?Mirim :?> float array)
//let a = Array.zip3 mz intensity ionMobility
//open Plotly.NET
//1
//Chart.Scatter3d(a, StyleParam.Mode.Markers)
//|> Chart.withX_AxisStyle "m/z"
//|> Chart.withY_AxisStyle "Intensity"
//|> Chart.withZ_AxisStyle "Voltage"
//|> Chart.withSize (1200., 900.)
//|> Chart.Show

//Chart.Point(mz,intensity)
//|> Chart.withX_AxisStyle "m/z"
//|> Chart.withY_AxisStyle "Intensity"
//|> Chart.withSize (1200., 900.)
//|> Chart.Show