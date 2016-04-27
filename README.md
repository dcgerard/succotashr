<!-- README.md is generated from README.Rmd. Please edit that file -->
SUCCOTASH: Surrogate and Confounder Correction Occuring Together with Adaptive Shrinkage
========================================================================================

[![Build Status](https://travis-ci.org/dcgerard/succotashr.svg?branch=master)](https://travis-ci.org/dcgerard/succotashr)

Description
-----------

Let

Y = XB + ZA + E,

for

-   Y an n by p matrix of gene expression data with n samples and p genes,
-   X an n by q matrix of of q covariates,
-   B a q by p matrix of unobserved coefficients for the observed covariates,
-   Z an n by k matrix of hidden confounders,
-   A an k by p matrix of hidden coefficients for the hidden confounders, and
-   E an n by p matrix of independent normal errors with column variances s1,...,sp.

Not accounting for the hidden covariates, Z, can reduce power and result in poor control of false discovery rate.

`succotashr` fits this model under a two-step empirical Bayesian approach. It places a non-parametric prior on B and jointly estimates B and ZA. The main function is `succotash`.

Installation
------------

To install, run the following code in `R`:

``` r
install.packages(c("devtools", "SQUAREM"))
devtools::install_github("stephens999/ashr")
devtools::install_github("dcgerard/succotashr")
```

The following packages are suggested:

``` r
install.packages("cate")
devtools::install_github("NKweiwang/flash")
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
```
