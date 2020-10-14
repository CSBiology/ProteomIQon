namespace ProteomIQon

module Ontologies =

    type SampleProcessing =
        |PeakPicking
        |TOPPPeakPickerWavelet
        |SavitzkyGolaySmoothing
        |WaveletTransformationSmoothing
        |LowIntensityDataPointRemoval
        |MedianBaselineReduction
        //|ParagonDigestion of string
        /// Combine this with the used enzyme
        |EnzymeDigestion
        |Trypsin
        |ArgC
        |Chymotrypsin
        |LysC
        |SDSPage
        |CoomassieStaining
        |ExcisedBand
        |SampleDesalting
        |DesaltingMicrocolumn
        |DesaltingCartridge
        |DesaltingBeads

        member this.toParam =
            match this with
            |PeakPicking ->
                "[MS, MS:1000035, peak picking, ]"
            |TOPPPeakPickerWavelet ->
                "[MS, MS:1002136, TOPP PeakPickerWavelet, ]"
            |SavitzkyGolaySmoothing ->
                "[MS, MS:1000782, Savitzky-Golay smoothing, ]"
            |WaveletTransformationSmoothing ->
                "[MS, MS:1001997, wavelet transformation smoothing, ]"
            |LowIntensityDataPointRemoval ->
                "[MS, MS:1000594, low intensity data point removal, ]"
            |MedianBaselineReduction ->
                "[MS, MS:1001996, median baseline reduction, ]"
            //|ParagonDigestion x ->
            //    sprintf "[MS, MS:1002436, Paragon: digestion, %s]" x
            |EnzymeDigestion ->
                "[SEP, sep:00142 , enzyme digestion, ]"
            |Trypsin ->
                "[MS, MS:1001251, Trypsin, ]"
            |ArgC ->
                "[MS, MS:1001303, Arg-C, ]"
            |Chymotrypsin ->
                "[MS, MS:1001306, Chymotrypsin, ]"
            |LysC ->
                "[MS, MS:1001309, Lys-C, ]"
            |SDSPage ->
                "[SEP, sep:00173, sodium dodecyl sulfate polyacrylamide gel electrophoresis, ]"
            |CoomassieStaining ->
                "[SEP, sep:00162, Coomassie staining, ]"
            |ExcisedBand ->
                "[SEP, sep:00132, excised band, ]"
            |SampleDesalting ->
                "[SEP, sep:00204, sample desalting, ]"
            |DesaltingMicrocolumn ->
                "[SEP, sep:00205, microcolumn, ]"
            |DesaltingCartridge ->
                "[SEP, sep:00206, cartridge, ]"
            |DesaltingBeads ->
                "[SEP, sep:00207, beads, ]"

    type Modification =
        |NoFixedModSearched
        |NoVariableModsSearched
        |Carbamidomethyl
        |Carboxymethyl
        |Oxidation
        |Acetyl
        |Amidated

        member this.toParam =
            match this with
            |NoFixedModSearched ->
                "[MS, MS:1002453, No fixed modifications searched, ]"
            |NoVariableModsSearched ->
                "[MS, MS:1002454, No variable modifications searched, ]"
            |Carbamidomethyl ->
                "[UNIMOD, UNIMOD:4, Carbamidomethyl, ]"
            |Carboxymethyl ->
                "[UNIMOD, UNIMOD:6, Carboxymethyl, ]"
            |Oxidation -> 
                "[UNIMOD, UNIMOD:35, Oxidation, ]"
            |Acetyl ->
                "[UNIMOD, UNIMOD:1, Acetyl, ]"
            |Amidated ->
                "[UNIMOD, UNIMOD:2, Amidated, ]"

    type ModificationPosition =
        |Anywhere
        |AnyNterm
        |AnyCterm
        |ProteinNterm
        |ProteinCterm

        member this.toParam =
            match this with
            |Anywhere ->
                "Anywhere"
            |AnyNterm ->
                "Any N-term"
            |AnyCterm ->
                "Any C-term"
            |ProteinNterm ->
                "Protein N-term"
            |ProteinCterm ->
                "Protein C-term"

    type SearchEngineScore =
        |Percolator
        |Mascot
        |SequestConsensusScore
        |SequestXCorr
        |SequestDeltaCn
        |Andromeda
        |XTandem
        |PeptideQValue
        |ProteinQValue

        member this.toParam =
            match this with
            |Percolator ->
                "[MS, MS:1001492, percolator:score, ]"
            |Mascot ->
                "[MS, MS:1001171, Mascot:score, ]"
            |SequestConsensusScore ->
                "[MS, MS:1001163, SEQUEST:consensus score, ]"
            |SequestXCorr ->
                "[MS, MS:1001155, SEQUEST:xcorr, ]"
            |SequestDeltaCn ->
                "[MS, MS:1001156, SEQUEST:deltacn, ]"
            |Andromeda ->
                "[MS, MS:1002338, Andromeda:score, ]"
            |XTandem ->
                "[MS, MS:1001331, X\\!Tandem:hyperscore, ]"
            |PeptideQValue ->
                "[MS, MS:1001868, distinct peptide-level q-value, ]"
            |ProteinQValue ->
                "[MS, MS:1001869, protein-level q-value, ]"

    type IDFormats =
        |Thermo
        |Waters
        |WIFF
        |BrukerAgilent
        |BrukerBAF
        |BrukerFID
        |MultiPeakListNativeID

        member this.toParam =
            match this with
            |Thermo ->
                "[MS, MS:1000768, Thermo nativeID format, ]"
            |Waters ->
                "[MS, MS:1000769, Waters nativeID format, ]"
            |WIFF ->
                "[MS, MS:1000770, WIFF nativeID format, ]"
            |BrukerAgilent ->
                "[MS, MS:1000771, Bruker/Agilent YEP nativeID format, ]"
            |BrukerBAF ->
                "[MS, MS:1000772, Bruker BAF nativeID format, ]"
            |BrukerFID ->
                "[MS, MS:1000773, Bruker FID nativeID format, ]"
            |MultiPeakListNativeID ->
                "[MS, MS:1000774, multiple peak list nativeID format, ]"

    type FileFormats =
        |MzML
        |BrukerAgilent
        |Thermo
        |WIFF

        member this.toParam =
            match this with
            |MzML ->
                "[MS, MS:1000584, mzML format, ]"
            |BrukerAgilent ->
                "[MS, MS:1000567, Bruker/Agilent YEP format, ]"
            |Thermo ->
                "[MS, MS:1000563, Thermo RAW format, ]"
            |WIFF ->
                "[MS, MS:1000562, ABI WIFF format, ]"
