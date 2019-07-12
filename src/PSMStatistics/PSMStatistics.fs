namespace ProteomIQon

open System.Data.SQLite
open BioFSharp
open BioFSharp.Mz.SearchDB
open Domain
open Core 
open Logary
open System.IO
open BioFSharp.Mz
open MzLite 

open Core.MzLite.Reader
open Core.MzLite.Peaks
open Core.MzLite
open FSharpAux.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Csv
open FSharpAux.IO.SchemaReader.Attribute

module PSMStatistics = 
    open System.IO

    type PercolatorIn = {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
        PSMId                        : string;
        Label                        : int;
        ScanNr                       : int;
        Charge2                      : int;
        Charge3                      : int;
        Charge4                      : int;
        Charge5                      : int;
        Charge6                      : int;
        Charge7                      : int;
        Charge8                      : int;
        PrecursorMZ                  : float;
        TheoMass                     : float;
        AbsDeltaMass                 : float;
        PeptideLength                : int;
        MissCleavages                : int;
        SequestScore                 : float;
        SequestNormDeltaBestToRest   : float;
        SequestNormDeltaNext         : float;
        AndroScore                   : float;
        AndroNormDeltaBestToRest     : float;
        AndroNormDeltaNext           : float;
        Peptide                      : string;
        Protein                      : string;
        }


    type PercolatorPSMOut = {
        // a combination of the spectrum ID in the rawFile, the ascending ms2 id and the chargeState in the search space seperated by '_'
        [<FieldAttribute(0)>]
        PSMId                        : string;
        [<FieldAttribute(1)>]
        PercolatorScore              : float
        [<FieldAttribute(2)>]
        QValue                       : float
        [<FieldAttribute(3)>]
        PosteriorErrorProbability    : float
        [<FieldAttribute(4)>]
        StringSequence               : string
        }
     
    let printParams argv = printfn "%A" argv
