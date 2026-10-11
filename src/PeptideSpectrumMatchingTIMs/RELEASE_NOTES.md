#### Unreleased
* Build the decoys by pseudo-reversal instead of reversal: the C-terminal residue stays in place and the rest is reversed, and a peptide that this leaves unchanged is reversed completely. The decoys keep the cleavage site of the targets. Their fragments are kept in a second index of the same size, which doubles the memory of the index. Target scores do not change, but the Expectscore and the NormDelta columns of targets do, as they are computed over targets and decoys
* Update BioFSharp.Mz to 0.2.3
* Refuse a peptide data base only when its residues exceed the longest byte array, not when the length of its sequence strings with the modification tokens does

#### 0.0.1 - Monday, September 28, 2026
* Initial release: fragment index search of timsTOF runs in mzlite files, with the decoys served from the target index
* Sum the mobility scans of every MS2 spectrum before the search, with the m/z width set by --mobility-merge-ppm
* Build the fragment index only for the peptides inside the precursor windows of the runs, in slices (-s), and keep the preprocessed spectra between the slices (-k)
* Search the runs of a group (-c) against one shared index, one thread per run
* Write the hyperscore, the X!Tandem expectation value and the SEQUEST-like, Andromeda-like and X!Tandem-like scores in the .psm layout of PeptideSpectrumMatching
