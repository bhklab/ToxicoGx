# *News*

## Future Relases

ToxicoGx >v1.0.0
* Continue to abstract functionality into CoreGx
* Add additional plotting functions such as grouped boxplots
* Extend coverage of unit tests to >90%
* Implement a faster version of drugPertubationSignature
* Add additional plotting functions
* Include scripts for differential expression analysis and GSEA of 
  toxico-genomic data (limma)

## Current Release

ToxicoGx v1.0.0
* Modified package to depend on updated CoreGx
* All molecularProfiles are now SummarizedExperiment instead of ExpressionSet
* Abstracted some additional functions to CoreGx

## Previous Releases

ToxicoGx v0.1.2
* Updated downloadTSet function to use published Zenodo DOIs to retrieve data
* Modified rankGeneDrugsPerturbation to fix a bad unit conversion which would return concentrations in the wrong unit

ToxicoGx v0.1.1
* Bug Fix: Regenerated TGGATESsmall (sample dataset) to fix make a result in the vignette consistent with previous releases.

ToxicoGx v0.1.0
* Rewrote plots using ggplot2 to improve aesthetics
  * Also can now extend plotting functions using standard ggplot2 syntax
* Improved package documentation  

## Intial Submission

ToxicoGx v0.0.1
* Minimal package submitted
