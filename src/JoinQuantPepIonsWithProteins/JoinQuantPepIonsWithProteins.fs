namespace ProteomIQon

open System.IO
open Argu
open Deedle
open FSharpAux
open FSharpAux.IO
open FSharpAux.IO.SchemaReader
open FSharpAux.IO.SchemaReader.Attribute
open ProteomIQon
open Dto


module JoinQuantPepIonsWithProteins = 
  
    ///
    let joinQuantPepIonsWithProteins outDirectory (quantifiedPeptides:string) (proteins:string) = //(param: TableSortParams) =
        let fileName = Path.GetFileNameWithoutExtension quantifiedPeptides
        // initialize Reader and Transaction
        let outFilePath =
            let fileName = (Path.GetFileNameWithoutExtension quantifiedPeptides) + ".quantAndProt"
            Path.Combine [|outDirectory;fileName|]
        let proteins = 
            Csv.CsvReader<Dto.ProteinInferenceResult>(SchemaMode=Csv.Fill).ReadFile(proteins,'\t',false,1)
            |> Array.ofSeq
            |> Frame.ofRecords
            |> Frame.indexRowsString "ProteinGroup"
            |> Frame.expandRowsByColumn "PeptideSequence" 
                (fun p (v:string) -> 
                    v.Split(';') 
                    |> Seq.map (fun s-> 
                        {| 
                            ProteinGroup = p
                            StringSequence = s 
                        |}
                        )
                )
            |> Frame.sliceCols ["ProteinGroup";"PeptideSequence";"QValue"]
            |> Frame.mapColKeys (fun ck -> if ck = "QValue" then "ProteinGroup_QValue" else ck)
            // 
        let peptideIons =
            Csv.CsvReader<Dto.QuantificationResult>(SchemaMode=Csv.Fill).ReadFile(quantifiedPeptides,'\t',false,1)
            |> Array.ofSeq
        let filteredPeptideIons =
            peptideIons
            |> Array.filter (fun qp -> Dto.QuantificationResult.tryTargetGetScanTime qp |> Option.isSome)
            |> Frame.ofRecords
            |> Frame.indexRowsUsing (fun s -> 
                {| 
                    StringSequence = s.GetAs<string>("StringSequence")
                    PepSequenceID = s.GetAs<int>("PepSequenceID")
                    ModSequenceID = s.GetAs<int>("ModSequenceID")
                    Charge = s.GetAs<int>("Charge") 
                    GlobalMod = s.GetAs<int>("GlobalMod")
                |}
                )
        let peptidesWithProteins =
            Frame.align 
                (fun (k: ({| ProteinGroup: string;StringSequence: string; |})) -> k.StringSequence ) 
                (fun (k: ({| StringSequence: string;PepSequenceID: int;ModSequenceID: int;Charge: int;GlobalMod: int; |})) -> k.StringSequence) 
                proteins filteredPeptideIons
            |> Frame.mapRowKeys (fun (s,prot,pep) -> 
                {| 
                    FileName = fileName
                    ProteinGroup    = prot.ProteinGroup
                    StringSequence  = pep.StringSequence
                    PepSequenceID   = pep.PepSequenceID
                    ModSequenceID   = pep.ModSequenceID
                    Charge          = pep.Charge
                    GlobalMod       = pep.GlobalMod
                |}
                )
        peptidesWithProteins
        |> Frame.mapRows (fun k s ->
            { 
            FileName                                    = k.FileName
            ProteinGroup                                = k.ProteinGroup
            ProteinGroup_QValue                         = s.GetAs<float>("ProteinGroup_QValue",nan)
            StringSequence                              = k.StringSequence
            PepSequenceID                               = k.PepSequenceID
            ModSequenceID                               = k.ModSequenceID
            Charge                                      = k.Charge
            GlobalMod                                   = k.GlobalMod
            PrecursorMZ                                 = s.GetAs<float>("PrecursorMZ",nan)
            MeasuredMass                                = s.GetAs<float>("MeasuredMass",nan)
            TheoMass                                    = s.GetAs<float>("TheoMass",nan)
            AbsDeltaMass                                = s.GetAs<float>("AbsDeltaMass",nan)
            MeanPercolatorScore                         = s.GetAs<float>("MeanPercolatorScore",nan)
            QValue                                      = s.GetAs<float>("QValue",nan)
            PEPValue                                    = s.GetAs<float>("PEPValue",nan)
            ProteinNames                                = s.GetAs<string>("ProteinNames","")
            QuantMz_Light                               = s.GetAs<float>("QuantMz_Light",nan)
            Quant_Light                                 = s.GetAs<float>("Quant_Light",nan)
            MeasuredApex_Light                          = s.GetAs<float>("MeasuredApex_Light",nan)
            Seo_Light                                   = s.GetAs<float>("Seo_Light",nan)
            Params_Light                                = s.GetAs<float[]>("Params_Light",[||])
            Difference_SearchRT_FittedRT_Light          = s.GetAs<float>("Difference_SearchRT_FittedRT_Light",nan)
            KLDiv_Observed_Theoretical_Light            = s.GetAs<float>("KLDiv_Observed_Theoretical_Light",nan)
            KLDiv_CorrectedObserved_Theoretical_Light   = s.GetAs<float>("KLDiv_CorrectedObserved_Theoretical_Light",nan)
            QuantMz_Heavy                               = s.GetAs<float>("QuantMz_Heavy",nan)
            Quant_Heavy                                 = s.GetAs<float>("Quant_Heavy",nan)
            MeasuredApex_Heavy                          = s.GetAs<float>("MeasuredApex_Heavy",nan)
            Seo_Heavy                                   = s.GetAs<float>("Seo_Heavy",nan)
            Params_Heavy                                = s.GetAs<float[]>("Params_Heavy",[||])        
            Difference_SearchRT_FittedRT_Heavy          = s.GetAs<float>("Difference_SearchRT_FittedRT_Heavy",nan)
            KLDiv_Observed_Theoretical_Heavy            = s.GetAs<float>("KLDiv_Observed_Theoretical_Heavy",nan)
            KLDiv_CorrectedObserved_Theoretical_Heavy   = s.GetAs<float>("KLDiv_CorrectedObserved_Theoretical_Heavy",nan)
            Correlation_Light_Heavy                     = s.GetAs<float>("Correlation_Light_Heavy",nan)
            QuantificationSource                        = s.GetAs<QuantificationSource>("QuantificationSource")
            IsotopicPatternMz_Light                     = s.GetAs<float[]>("IsotopicPatternMz_Light",[||]) 
            IsotopicPatternIntensity_Observed_Light     = s.GetAs<float[]>("IsotopicPatternIntensity_Observed_Light",[||]) 
            IsotopicPatternIntensity_Corrected_Light    = s.GetAs<float[]>("IsotopicPatternIntensity_Corrected_Light",[||]) 
            RtTrace_Light                               = s.GetAs<float[]>("RtTrace_Light",[||]) 
            IntensityTrace_Observed_Light               = s.GetAs<float[]>("IntensityTrace_Observed_Light",[||]) 
            IntensityTrace_Corrected_Light              = s.GetAs<float[]>("IntensityTrace_Corrected_Light",[||]) 
            IsotopicPatternMz_Heavy                     = s.GetAs<float[]>("IsotopicPatternMz_Heavy",[||]) 
            IsotopicPatternIntensity_Observed_Heavy     = s.GetAs<float[]>("IsotopicPatternIntensity_Observed_Heavy",[||]) 
            IsotopicPatternIntensity_Corrected_Heavy    = s.GetAs<float[]>("IsotopicPatternIntensity_Corrected_Heavy",[||]) 
            RtTrace_Heavy                               = s.GetAs<float[]>("RtTrace_Heavy",[||]) 
            IntensityTrace_Observed_Heavy               = s.GetAs<float[]>("IntensityTrace_Observed_Heavy",[||]) 
            IntensityTrace_Corrected_Heavy              = s.GetAs<float[]>("IntensityTrace_Corrected_Heavy",[||]) 
            AlignmentScore                              = s.GetAs<float>("AlignmentScore", nan)
            AlignmentQValue                             = s.GetAs<float>("AlignmentQValue", nan)
            }
            )
        |> Series.values
        |> SeqIO'.csv "\t" true false
        |> SeqIO.Seq.writeOrAppend (outFilePath)

    
 
