# Junction coverage compatibility scores

[![R build status](https://github.com/csoneson/jcc/workflows/R-CMD-check/badge.svg)](https://github.com/csoneson/jcc/actions)
[![No Maintenance Intended](http://unmaintained.tech/badge.svg)](http://unmaintained.tech/)
<!--[![Codecov.io coverage status](https://codecov.io/github/csoneson/jcc/coverage.svg?branch=master)](https://codecov.io/github/csoneson/jcc)-->


This repository contains an R package aimed at facilitating the calculation of gene-wise junction coverage compatibility (JCC) scores, as described in 

- Soneson C, Love MI, Patro R, Hussain S, Malhotra D and Robinson MD: A junction coverage compatibility score to quantify the reliability of transcript abundance estimates and annotation catalogs. [Life Science Alliance 2019](http://www.life-science-alliance.org/content/2/1/e201800175.abstract).

The JCC scores are aimed at detecting genes for which estimated transcript abundances are unreliable, either because of problems in the transcript abundance estimation or because of missing or wrongly annotated reference transcripts. It does so by comparing the number of reads predicted to align across each junction, inferred from the transcript abundances and a fragment bias model, to the observed number of junction-spanning reads, obtained via alignment of the reads to the genome. A high JCC score for a gene indicates that the estimated abundances for the corresponding transcripts are unreliable and should be treated with care in downstream analyses. 

## Installation

The `jcc` package can be installed using `devtools`:

```
library(devtools)
install_github("csoneson/jcc")
```

or the `BiocManager` CRAN package

```
BiocManager::install("csoneson/jcc")
```

Note that `jcc` depends on the `alpine` package, which was removed from Bioconductor in release 3.18. Thus, in order to use `jcc`, please install Bioconductor release 3.17 or earlier. 
