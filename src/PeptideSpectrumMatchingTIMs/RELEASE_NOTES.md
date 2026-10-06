#### 0.0.1 - Monday, September 28, 2026
* Initial release: fragment index search of timsTOF runs in mzlite files, with the decoys served from the target index
* Sum the mobility scans of every MS2 spectrum before the search, with the m/z width set by --mobility-merge-ppm
* Build the fragment index only for the peptides inside the precursor windows of the runs, in slices (-s), and keep the preprocessed spectra between the slices (-k)
* Search the runs of a group (-c) against one shared index, one thread per run
* Write the hyperscore, the X!Tandem expectation value and the SEQUEST-like, Andromeda-like and X!Tandem-like scores in the .psm layout of PeptideSpectrumMatching
