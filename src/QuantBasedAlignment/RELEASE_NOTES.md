#### Unreleased
* Model the scan time shift between two runs with a running median over the k nearest shared ions, after three rounds of trimming residuals beyond max(3 MAD, 0.05 min); k is chosen by 5-fold cross-validation grouped by peptide (PepSequenceID)
* Rank the source runs of a target by their median absolute out-of-fold scan time error and transfer each ion from the best ranked source that has it
* Add columns to the .align file after the existing ones: the local error scale of the predicted scan time, the source run, the source model error and the predicted ion mobility with its error scale
* Predict the target ion mobility with the same model on the IonMobilityCentre of the .quant rows when both runs have ion mobility; the IonMobility column of the .align and .alignmetric files holds the ion mobility centre of the ion in its source run
* Give ions the model was trained on their out-of-fold prediction and a local error scale without their own error
* Write the out-of-fold predictions of one cross-validation fold to .alignmetric, without dynamic time warping
* Stop with an error message on input paths that do not exist or have no .quant file, and on input files from different directories with the same name
* Warn once per target that has no source with enough shared ions
* Seed the cross-validation folds (-s, default 1)
* Log the ions that are not used with their reason
* Overwrite existing output files instead of appending to them
* Replace the packages ProteomIQon.QuantBasedAlignment_win_x64 and ProteomIQon.QuantBasedAlignment_linux_x64 by the tool package ProteomIQon.QuantBasedAlignment, until now a library, which has no MKL dependency and runs on all platforms
* Take the scan time of a fitted peak at the apex of its fitted shape and its width as the standard deviation of that shape, also for peaks with an exponential tail
* Breaking: the .quant files need the column IonMobilityCentre, so files of earlier PSMBasedQuantification and PSMBasedQuantificationTIMs releases have to be quantified again

#### 0.0.4 - Wednesday, September 2, 2026
* Update to .NET 10
* Update System.Data.SQLite.Core to 1.0.119 for the memory access fix
* Update BioFSharp to 2.0.0, FSharpAux to 2.1.0, FSharp.Stats to 0.6.0 and Plotly.NET to 6.0.0-preview.2

#### 0.0.3 - Sunday, March 23, 2025
* Add Ion Mobility to output type
* Update to .NET 8

#### 0.0.2 - Monday, February 13, 2023
* Add Alignment Score fields to quantification results

#### 0.0.1 - Monday, November 8, 2021
* Initial nuget release
