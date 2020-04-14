namespace ProteomIQon


open System.IO
open Argu
open System.Data.SQLite
open System
open ProteomIQon.Core
open Core.MzIO
open Dto
open FSharp.Stats
open BioFSharp.Mz
open FSharpAux.IO.SchemaReader
open FSharp.Plotly
open BioFSharp
open MzIO.Processing
open Microsoft
open Microsoft.ML
open Microsoft.ML.Data   

module QuantBasedAlignment = 

    let printParams argv = printfn "%A" argv

    /////
    //type AlignmentResult = 
    //    {
    //    StringSequence               : string;
    //    GlobalMod                    : int;
    //    Charge                       : int;
    //    N14QuantMz                   : float
    //    N15QuantMz                   : float
    //    N15Minus1QuantMz             : float
    //    PredictedScanTime            : float
    //    ProteinNames                 : string;
    //    }  
    
    //type AlignmentItem = {
    //    Sequence             : string
    //    GlobalMod            : bool
    //    Charge               : int
    //    Mz                   : float          
    //    MeanScanTime         : float
    //    Stabw                : float
    //    Intensity            : float
    //    N14Quant             : float
    //    N15Quant             : float
    //    N14QuantMz           : float
    //    N15QuantMz           : float
    //    N15Minus1QuantMz     : float
    //    Ratio                : float
    //    }

    //type AlignmentItemForLearning = {
    //    Sequence             : string
    //    GlobalMod            : bool
    //    Charge               : int
    //    Mz                   : float          
    //    MeanScanTime         : float
    //    CorrespondingScanTime: float
    //    Stabw                : float
    //    Intensity            : float
    //    N14Quant             : float
    //    N15Quant             : float
    //    N14QuantMz           : float
    //    N15QuantMz           : float
    //    N15Minus1QuantMz     : float
    //    Ratio                : float
    //    }

    //let alignmentItemForLearningOf correspondingScanTime (aItem:AlignmentItem) =
    //    {
    //    Sequence             = aItem.Sequence             
    //    GlobalMod            = aItem.GlobalMod            
    //    Charge               = aItem.Charge               
    //    Mz                   = aItem.Mz                            
    //    MeanScanTime         = aItem.MeanScanTime         
    //    CorrespondingScanTime= correspondingScanTime
    //    Stabw                = aItem.Stabw                
    //    Intensity            = aItem.Intensity            
    //    N14Quant             = aItem.N14Quant             
    //    N15Quant             = aItem.N15Quant             
    //    N14QuantMz           = aItem.N14QuantMz           
    //    N15QuantMz           = aItem.N15QuantMz           
    //    N15Minus1QuantMz     = aItem.N15Minus1QuantMz     
    //    Ratio                = aItem.Ratio                
    //    }


    //type Peptide = {
    //    Sequence             : string
    //    GlobalMod            : int
    //    Charge               : int
    //    ScanTime             : float
    //    StandardDeviation    : float
    //    }

    //type AlignmentFile<'a> = {
    //    FileName                   : string
    //    Peptides                   : AlignmentItem [] 
    //    MissingPeptides            : Peptide []
    //    }

    //type AlignmentParams = {
    //    Placeholder : bool 
    //    }

    /////
    //let getQuantifiedPeptides (quantFilePath:string) = 
    //    ///
    //    let peptides =
    //        Csv.CsvReader<QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(quantFilePath,'\t',false,1)
    //        |> Array.ofSeq
    //    peptides
    //    |> Array.filter (fun qp -> 
    //        // Filter for peptides where fit assures a good estimation of the ScanTime
    //        (qp.GlobalMod = 0 && qp.N14Params <> "") || (qp.GlobalMod = 1 && qp.N15Params <> "")
    //        )
    //    |> Array.map (fun qp ->
                
    //            if qp.GlobalMod = 0 then 
    //                {
    //                Sequence          = qp.StringSequence
    //                GlobalMod         = qp.GlobalMod
    //                Charge            = qp.Charge
    //                ScanTime          = qp.N14Params
    //                StandardDeviation = qp.N14Params
    //                }
    //            else
    //                {
    //                Sequence          = qp.StringSequence
    //                GlobalMod         = qp.GlobalMod
    //                Charge            = qp.Charge
    //                ScanTime          = qp.N14Params
    //                StandardDeviation = qp.N14Params
    //                }
    //        )        


    //let alignFiles (processParams:AlignmentParams) (*(filesToAlign:string [])*) (filePath:string) = 

    //    let logger = Logging.createLogger (Path.GetFileNameWithoutExtension instrumentOutput)

    //    logger.Trace (sprintf "Input directory: %s" filePath)
    //    logger.Trace (sprintf "Output directory: %s" outputDir)
    //    logger.Trace (sprintf "Parameters: %A" processParams)
        
    //    let quantFiles = 
    //        System.IO.Directory.GetFiles (filePath, "*.quant")
    //        |> Array.toList

        

    //    let ctx = new ML.MLContext()
        




































