#### Unreleased
* Breaking: the tool reads the ion table of LabeledProteinQuantification with -i, writes RatioLFQ.txt and a log file into the directory given with -o and takes a parameter file with -p. The flags -pc, -rac and -rfc are removed.
* Compute the ratios from the peptide ion table of LabeledProteinQuantification by comparing the runs pairwise over shared ions
* Write a light to heavy ratio and an intensity per protein group and run
* Add a parameter file with the reference channel, the minimum number of shared ions and the centring of runs
* Read quoted run names and protein groups without their quotes, as LabeledProteinQuantification quotes run names with a space and protein groups with a semicolon
* Stop with a message when the largest intensity of a channel is below 100, since the table then looks log transformed
* Update NLog to 4.7.11

#### 0.0.2 - Wednesday, September 2, 2026
* Update to .NET 10
* Update FSharp.Stats to 0.6.0

#### 0.0.1 - Sunday, March 23, 2025
* Initial release
