#### Unreleased
* Model the scan time shift between two runs with a running median over the k nearest shared ions, after three rounds of trimming residuals beyond max(3 MAD, 0.05 min); k is chosen by 5-fold cross-validation grouped by peptide (PepSequenceID)
* Rank the source runs of a target by their median absolute out-of-fold scan time error and transfer each ion from the best ranked source that has it
* Add columns to the .align file after the existing ones: the local error scale of the predicted scan time, the source run, the source model error and the predicted ion mobility with its error scale
* Predict the target ion mobility with the same model when both runs have ion mobility
* Give ions the model was trained on their out-of-fold prediction and a local error scale without their own error
* Write the out-of-fold predictions of one cross-validation fold to .alignmetric, without dynamic time warping
* Stop with an error message on input paths that do not exist or have no .quant file, and on input files from different directories with the same name
* Warn once per target that has no source with enough shared ions
* Seed the cross-validation folds (-s, default 1)
* Log the ions that are not used with their reason
* Overwrite existing output files instead of appending to them
* No MKL dependency, so no platform specific packages
