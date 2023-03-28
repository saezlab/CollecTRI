# CollecTRI: **Collec**tion of **T**ranscriptional **R**egulatory **I**nteractions <img src="man/figures/CollecTRI_logo.png" align="right" width="120" />

<!-- badges: start -->
<!-- badges: end -->

## Overview
The CollecTRI-derived regulons contain signed transcription factor (TF) - target gene 
interactions compiled from 12 different resources. This collection provides 
an increased coverage of transcription factors and was benchmarked against 
other known GRNs, showing a superior performance in identifying perturbed TFs 
based on gene expression data using the knockTF data sets.


Graphical abstract


## Resources included in CollecTRI
- ExTRI
- HTRI
- TRRUST
- TFActS
- IntAct
- SIGNOR
- CytReg
- GEREDB
- Pavlidis
- DoRothEA A
- NTNU curations


## Data availability 
CollecTRI is available through the
[OmniPath](https://omnipathdb.org/) or  
[DoRothEA](https://saezlab.github.io/dorothea/) packages.


## Scripts
For more information about CollecTRI, please check out the following scripts:

- [Construction of CollecTRI](https://github.com/saezlab/CollecTRI/tree/main/scripts/CollecTRI)
- [Benchmark](https://github.com/saezlab/CollecTRI/tree/main/scripts/benchmark)
  - [Systematic comparison](https://github.com/saezlab/CollecTRI/blob/main/scripts/benchmark/benchmark.ipynb)
  - [Evaluation of binding weights](https://github.com/saezlab/CollecTRI/blob/main/scripts/benchmark/benchmark_weights.ipynb)
  - [Statistical evaluation](https://github.com/saezlab/CollecTRI/blob/main/scripts/benchmark/statistics.R)
- [Case study](https://github.com/saezlab/CollecTRI/blob/main/scripts/Case_study/case_study.R)
- [Manuscript figures](https://github.com/saezlab/CollecTRI/blob/main/scripts/figures/figures_manuscript.R)



## License
The CollecTRI-derived regulons are freely available to the community. The original licenses of all 
resources included in CollecTRI can be found [here](https://github.com/saezlab/pypath/blob/master/pypath/resources/data/resources.json)


## Citation
> 
