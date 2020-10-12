namespace ProteomIQon

module ParamMappings =

    type SampleProcessing =
        |PeakPicking
        |TOPPPeakPickerWavelet
        |SavitzkyGolaySmoothing
        |WaveletTransformationSmoothing
        |LowIntensityDataPointRemoval
        |MedianBaselineReduction
        |ParagonDigestion of string
        |EnzymeDigestion
        |Trypsin
        |ArgC
        |Chymotrypsin
        |LysC
        |SDSPage
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
            |ParagonDigestion x ->
                sprintf "[MS, MS:1002436, Paragon: digestion, %s]" x
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
            |SampleDesalting ->
                "[SEP, sep:00204, sample desalting, ]"
            |DesaltingMicrocolumn ->
                "[SEP, sep:00205, microcolumn, ]"
            |DesaltingCartridge ->
                "[SEP, sep:00206, cartridge, ]"
            |DesaltingBeads ->
                "[SEP, sep:00207, beads, ]"
