namespace ProteomIQon

open Domain
open Core
open System.IO
open MzIO
open MzIO.Model
open MzIO.Processing
open MzIO.IO
open MzIO.IO.MzML

module MzMLConverter =

    let initGetData (reader:IMzIODataReader)=
        fun (massSpec:MassSpectrum) ->
            reader.ReadSpectrumPeaks(massSpec.ID).Peaks
            |> Core.MzIO.Peaks.unzipIMzliteArray

    let convertFile (converterParams: MzMLConverterParams)(outputDir:string) (instrumentOutput:string) =

        let logger = Logging.createLogger (Path.GetFileNameWithoutExtension instrumentOutput)

        logger.Trace (sprintf "Input file: %s" instrumentOutput)
        logger.Trace (sprintf "Output directory: %s" outputDir)
        logger.Trace (sprintf "Parameters: %A" converterParams)

        logger.Trace "Init connection to input data base."
        // initialize Reader and Transaction
        let inReader = Core.MzIO.Reader.getReader instrumentOutput :?> MzSQL.MzSQL
        let cn = inReader.Open()
        let inTr = inReader.BeginTransaction()
        let inRunID  = Core.MzIO.Reader.getDefaultRunID inReader

        logger.Trace "Creating mzML file."
        // initialize Reader and Transaction
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension instrumentOutput) + ".mzML"
            Path.Combine [|outputDir;fileName|]

        let outReader = new MzMLWriter(outFilePath)

        /// All files created by this application will have a unified runID.
        let outRunID  = "sample=0"

        logger.Trace "Getting mass spectra."
        // Get all mass spectra
        let massSpectra = inReader.ReadMassSpectra(inRunID)
        logger.Trace "Filtering mass spectra according to retention time."
        let massSpectraF =
            match converterParams.StartRetentionTime, converterParams.EndRetentionTime with
            | Some s, Some e ->
                massSpectra
                |> Seq.filter (fun ms ->
                                let scanTime = MassSpectrum.getScanTime ms
                                scanTime > s && scanTime < e
                              )
            | Some s, None ->
                massSpectra
                |> Seq.filter (fun ms ->
                                let scanTime = MassSpectrum.getScanTime ms
                                scanTime > s
                              )
            | None, Some e ->
                massSpectra
                |> Seq.filter (fun ms ->
                                let scanTime = MassSpectrum.getScanTime ms
                                scanTime < e
                              )
            | None, None ->
                massSpectra

        logger.Trace (sprintf "Writing %i mass spectra to output mzML file." (Seq.length massSpectraF))
        outReader.insertMSSpectraBy (outReader.insertMSSpectrum) outRunID inReader converterParams.Compress (massSpectraF)

        inTr.Commit()
        inTr.Dispose()
        inReader.Dispose()