namespace ProteomIQon


open Deedle
open BioFSharp
open FSharp.Stats
open Plotly.NET
open FSharpAux
open System.IO

module LabelEfficiencyCalculator = 
  
    type PeptideIon = 
        {
            ProteinGroup  : string  
            StringSequence: string
            PepSequenceID : int
            Charge        : int
            GlobalMod     : bool
        }
        static member create proteinGroup stringSequence pepSequenceID charge globalMod =
            {
                ProteinGroup   = proteinGroup
                StringSequence = stringSequence
                PepSequenceID  = pepSequenceID
                Charge         = charge
                GlobalMod      = globalMod
            }

    // create chemical formula for amino acid and add water to reflect hydrolysed state in mass spectrometer
    let toFormula bioseq =  
        bioseq
        |> BioSeq.toFormula
        // peptides are hydrolysed in the mass spectrometer, so we add H2O
        |> Formula.add Formula.Table.H2O

    let label n15LableEfficiency formula =
        let heavyN15 = Elements.Di  (Elements.createDi "N15" (Isotopes.Table.N15,n15LableEfficiency) (Isotopes.Table.N14,1.-n15LableEfficiency) )
        Formula.replaceElement formula Elements.Table.N heavyN15

    // Predicts an isotopic distribution of the given formula at the given charge, 
    // normalized by the sum of probabilities, using the MIDAs algorithm
    let generateIsotopicDistribution (charge:int) (f:Formula.Formula) =
        IsotopicDistribution.MIDA.ofFormula 
            IsotopicDistribution.MIDA.normalizeByProbSum
            0.01
            0.000000001
            charge
            f

    let floatArrayOf s = 
        if String.isNullOrEmpty s then 
            None
        else
            s
            |> String.split (';') 
            |> Array.map float
            |> Some

    let getHeavyPattern (heavyPatternMz: string) (heavyPatternI: string) = 
        let mz = 
            heavyPatternMz 
            |> floatArrayOf
        let intensities = 
            heavyPatternI 
            |> floatArrayOf
        match mz, intensities with
        | Some mz', Some intensities' ->
            Array.zip mz' intensities'
            |> Some
        | _ -> None

    type ExtractedIsoPattern = 
        {
            PeptideInfo: PeptideIon
            Pattern    : array<(float*float)>
        }
        static member create peptideInfo pattern =
           {
               PeptideInfo = peptideInfo
               Pattern     = pattern
           } 

    type SimulatedIsoPattern = 
        {
            PeptideSequence: string
            Charge         : int
            LableEfficiency: float
            SimPattern     : list<(float*float)>
        }
        static member create peptideSequence charge lableEfficiency simPattern =
            {
                PeptideSequence = peptideSequence
                Charge          = charge
                LableEfficiency = lableEfficiency
                SimPattern      = simPattern
            }

    let simulateFrom peptideSequence charge lableEfficiency =
        let simPattern =
            peptideSequence
            |> BioSeq.ofAminoAcidString
            |> toFormula 
            |> label lableEfficiency
            |> generateIsotopicDistribution charge 
        SimulatedIsoPattern.create
            peptideSequence
            charge
            lableEfficiency
            simPattern

    let plotIsotopicPattern color mzsAndintensities =
        let min,max =
            mzsAndintensities |> Seq.minBy fst |> fst,
            mzsAndintensities |> Seq.maxBy fst |> fst
        Seq.map (fun (x,y) -> 
            Chart.Line([x;x],[0.;y], ShowLegend = false, Opacity = 0.5)
            |> Chart.withLineStyle (Width = 7)
        ) mzsAndintensities
        |> Chart.combine
        |> Chart.withMarkerStyle(Size=0,Color = Color.fromHex (FSharpAux.Colors.toWebColor color))
        |> Chart.withXAxisStyle ("m/z", MinMax = (min - 1., max + 1.))
        |> Chart.withYAxisStyle "relative Intensity"

    let normBySum (a:seq<float*float>) =
        let s = Seq.sumBy snd a 
        Seq.map (fun (x,y) -> x,y / s) a

    /// Calculates the Kullback-Leibler divergence Dkl(p||q) from q (theory, model, description, or approximation of p) 
    /// to p (the "true" distribution of data, observations, or a ___ ___ precisely measured).
    let klDiv (p:seq<float>) (q:seq<float>) = 
        Seq.fold2 (fun acc p q -> (System.Math.Log(p/q)*p) + acc ) 0. p q

    let compareIsotopicDistributions (measured:ExtractedIsoPattern) (simulated:SimulatedIsoPattern)= 
        let patternSim = 
            measured.Pattern 
            |> Seq.map (fun (mz,intensities) -> 
                    mz,
                    simulated.SimPattern
                    |> Seq.filter (fun (mzSim,intensitiesSim) -> abs(mzSim-mz) < 0.05 )
                    |> Seq.sumBy snd
                )
            |> normBySum
        let klDiv = klDiv (patternSim |> Seq.map snd)  (measured.Pattern |> Seq.map snd)
        klDiv

    let compareIsotopicDistributionsPlot (measured:ExtractedIsoPattern) (simulated:SimulatedIsoPattern) =
        let patternSim = 
            measured.Pattern 
            |> Seq.map (fun (mz,intensities) -> 
                    mz,
                    simulated.SimPattern
                    |> Seq.filter (fun (mzSim,intensitiesSim) -> abs(mzSim-mz) < 0.05 )
                    |> Seq.sumBy snd
                )
            |> normBySum
        [
            plotIsotopicPattern FSharpAux.Colors.Table.Office.blue measured.Pattern
            |> Chart.withTraceInfo "Measured"
            plotIsotopicPattern FSharpAux.Colors.Table.Office.orange patternSim
            |> Chart.withTraceInfo "Simulated"
        ]
        |> Chart.combine

    let calcKL extractedIsoPattern peptideSequence charge lableEfficiency = 
        let sim = simulateFrom peptideSequence charge lableEfficiency
        let comp = compareIsotopicDistributions extractedIsoPattern sim
        comp
    
    let readQuantAndProt (path: string) =
        Frame.ReadCsv(path, true, separators = "\t")
        |> Frame.indexRowsUsing (fun s ->
            PeptideIon.create
                (s.GetAs<string>("ProteinGroup"))
                (s.GetAs<string>("StringSequence"))
                (s.GetAs<int>("PepSequenceID"))
                (s.GetAs<int>("Charge"))
                (s.GetAs<bool>("GlobalMod"))
        )

    let calculateLE (expectedLELower,expectedLEUpper) (charts: bool) (filename: string) (frame: Frame<PeptideIon,string>) =
        let plotDir = Path.Combine [|Path.GetDirectoryName filename; Path.GetFileNameWithoutExtension filename|]
        if charts then
            Directory.CreateDirectory plotDir
            |> ignore
        frame
        |> Frame.mapRows (fun k s ->
            if k.GlobalMod && k.StringSequence.Contains("[") |> not then
                let heavyPatternOption =
                    let mz =
                        s.GetAs<string>("IsotopicPatternMz_Heavy")
                    let intensity =
                        s.GetAs<string>("IsotopicPatternIntensity_Corrected_Heavy")
                    let heavyPattern =
                        getHeavyPattern mz intensity
                    heavyPattern
                match heavyPatternOption with
                | Some heavyPattern ->
                    let extractedPattern =
                        ExtractedIsoPattern.create
                            k
                            heavyPattern
                    let labelEfficiency = 
                        FSharp.Stats.Optimization.Brent.minimize 
                            (calcKL extractedPattern k.StringSequence k.Charge)
                            expectedLELower
                            expectedLEUpper
                    if charts && labelEfficiency.IsSome then
                        let replaceAsterisk = k.StringSequence |> String.replace "*" ""
                        compareIsotopicDistributionsPlot
                            extractedPattern
                            (simulateFrom k.StringSequence k.Charge labelEfficiency.Value)
                        |> Chart.withTitle $"{k.StringSequence}; Charge: {k.Charge}; LabelEfficiency: {labelEfficiency.Value}"
                        |> Chart.saveHtml (Path.Combine [|plotDir; $"{replaceAsterisk}_{k.Charge}.html"|])
                    match labelEfficiency with
                    | Some le ->
                        (le, (calcKL extractedPattern k.StringSequence k.Charge le))
                        |> Some
                    | None -> None
                | None -> None
            else
                None
        )
        
    let visualizeLE (path: string) (frame: Frame<PeptideIon,string>) =
        frame
        |> Frame.getCol "LabelEfficiency"
        |> Series.values
        |> fun (x: seq<float>) -> Chart.BoxPlot (X = x)
        |> Chart.withTraceInfo ("")
        |> Chart.withXAxisStyle ("Label Efficiency", MinMax = (0.9, 1.0))
        |> Chart.withTitle $"{path |> System.IO.Path.GetFileNameWithoutExtension}"
        |> Chart.saveHtml (path + ".html")

    let addLE (expectedLELower,expectedLEUpper) (charts: bool) (path: string) =
        let frame =
            path
            |> readQuantAndProt
            |> Frame.take 100
        let leKlColumn =
            frame
            |> calculateLE (expectedLELower,expectedLEUpper) (charts: bool) path
        let leColumn =
            leKlColumn
            |> Series.map (fun k s ->
                match s with
                | Some (le,kl) -> OptionalValue<float> le
                | None -> OptionalValue.Missing
            )
        let klColumn =
            leKlColumn
            |> Series.map (fun k s ->
                match s with
                | Some (le,kl) -> OptionalValue<float> kl
                | None -> OptionalValue.Missing
            )
        let finalFrame =
            frame
            |> Frame.addCol "LabelEfficiency" leColumn
            |> Frame.addCol "KullbackLeiblerDivergence" klColumn
        visualizeLE path finalFrame
        finalFrame
        |> fun fr -> fr.SaveCsv(path + "LE",false, separator = '\t')
