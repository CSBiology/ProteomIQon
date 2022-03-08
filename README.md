[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6335068.svg)](https://doi.org/10.5281/zenodo.6335068)
[![Made with F#](https://img.shields.io/badge/Made%20with-FSharp-rgb(184,69,252).svg)](https://fsharp.org/)
![GitHub contributors](https://img.shields.io/github/contributors/CSBiology/ProteomIQon)

The ProteomIQon 
----------------------

The ProteomIQon is a collection of open source computational proteomics tools to build pipelines for the evaluation of MS derived proteomics data
written in F#. Each tool is described in detail on its corresponding [documentation page](http://csbiology.github.io/ProteomIQon).

<img src="https://github.com/CSBiology/ProteomIQon/blob/master/docs/img/PillarsOfCompProt.png" width="750" height="400" />

The Core Project
------------------

The ProteomIQon core is referenced by all tools. It contains mainly serializable data transfer objects such as tool results and tool parameters, as well as their mapping to domain specific types. This is also the place for any kind of code reusable across tools such as thin wrappers around data readers, logging or CLI formatting. 

Documentation
-------------

The documentation and tutorials for this library are automatically generated (using the F# Formatting) from *.fsx and *.md files in the docs folder. If you find a typo, please submit a pull request!

Contributing
------------

Please refer to the CSB [Contribution guidelines](.github/CONTRIBUTING.md)

Community/Social
----------------

Want to get in touch with us? We recently joined the twitter crowd:

[![Twitter Follow](https://img.shields.io/twitter/follow/cs_biology.svg?style=social)](https://twitter.com/cs_biology)

Citation
----------------

When using ProteomIQon in scientifc or commercial releases, please cite us using the following DOI: `10.5281/zenodo.6335068`

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6335068.svg)](https://doi.org/10.5281/zenodo.6335068)
