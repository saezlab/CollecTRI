# CollecTRI: **Collec**tion of **T**ranscriptional **R**egulatory **I**nteractions <img src="man/figures/CollecTRI_logo.png" align="right" width="120" />

<!-- badges: start -->
<!-- badges: end -->

## Overview
The CollecTRI-derived regulons contain signed transcription factor (TF) - target gene 
interactions compiled from 12 different resources. This collection provides 
an increased coverage of transcription factors and was benchmarked against 
other known GRNs, showing a superior performance in identifying perturbed TFs 
based on gene expression data using the knockTF data sets.

<p align="center" width="100%">
<img src="man/figures/overview.png" align="center" width="550">
</p>

## Data availability 
The CollecTRI regulons are available in the [DoRothEA](https://saezlab.github.io/dorothea/) and [decoupler](https://saezlab.github.io/decoupleR/) packages through [OmniPath](https://omnipathdb.org/).
A tutroial on how to perform TF activity estimation using CollecTRI is available in [R](https://saezlab.github.io/decoupleR/articles/tf_bk.html) and [python](https://decoupler-py.readthedocs.io/en/latest/notebooks/dorothea.html). 

To load the CollecTRI regulons through R or python you can use the following lines:
```r
# processed regulons
decoupleR::get_collectri(organism='human', split_complexes=FALSE)

# raw regulons
OmnipathR::collectri(organism=9606L, genesymbols=TRUE, loops=TRUE)
```

```python
# processed regulons
import decoupler as dc
dc.get_collectri(organism='human', split_complexes=False)

# raw regulons
import omnipath as op
op.interactions.CollecTRI.get(genesymbols=True, organism=9606L, loops=True)
```

## Resources included in CollecTRI
ExTRI, HTRI, TRRUST, TFActS, IntAct, SIGNOR, CytReg, GEREDB, Pavlidis, DoRothEA A, NTNU curations

## Scripts
For more information about the CollecTRI-derived regulons, please check out the following scripts:

- [Construction of CollecTRI regulons](https://github.com/saezlab/CollecTRI/tree/main/scripts/CollecTRI)
- [Benchmark](https://github.com/saezlab/CollecTRI/tree/main/scripts/benchmark)
  - [Systematic comparison](https://github.com/saezlab/CollecTRI/blob/main/scripts/benchmark/02_benchmark.ipynb)
  - [Evaluation of binding weights](https://github.com/saezlab/CollecTRI/blob/main/scripts/benchmark/03_benchmark_weights.ipynb)
  - [Statistical evaluation](https://github.com/saezlab/CollecTRI/blob/main/scripts/benchmark/05_statistics.R)
- [Case study](https://github.com/saezlab/CollecTRI/blob/main/scripts/casestudy/case_study.R)
- [Manuscript figures](https://github.com/saezlab/CollecTRI/blob/main/scripts/figures/figures_manuscript.R)

If you are interested in the construction of the CollecTRI meta-resource check
out this [repository](https://github.com/Rbbt-Workflows/ExTRI)

## License
The CollecTRI-derived regulons are freely available for academic use. For commercial use please remove TF-gene interactions from SIGNOR.
The original licenses of all resources included in CollecTRI can be found [here](https://github.com/saezlab/pypath/blob/master/pypath/resources/data/resources.json)


## Citation
> Müller-Dott, S., Tsirvouli, E., Vazquez, M., Ramirez Flores, R. O., Badia-I-Mompel, P., Fallegger, R., Türei, D., Lægreid, A., & Saez-Rodriguez, J. (2023).
> Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities.
> Nucleic Acids Research. https://doi.org/10.1093/nar/gkad841
