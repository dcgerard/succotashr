<!-- README.md is generated from README.Rmd. Please edit that file -->
Description
-----------

Let \[
Y = X\beta + Z\alpha + E,
\] for

-   \(Y \in \mathbb{R}^{n\times p}\), a matrix of gene expression data with \(n\) samples and \(p\) genes,
-   \(X \in \mathbb{R}^{n \times q}\), a matrix of of \(q\) covariates,
-   \(\beta \in \mathbb{R}^{q \times p}\), the unobserved matrix of coefficients for the observed covariates,
-   \(Z \in \mathbb{R}^{n \times k}\), a matrix of hidden confounders,
-   \(\alpha \in \mathbb{R}^{k \times p}\), the matrix of hidden coefficients for the hidden confounders, and
-   \(E \sim N_{n \times p}(0, \Sigma \otimes I_n)\), the error matrix following a matrix normal with identity row covariance and diagonal column covariance \(\Sigma = diag(\sigma_1^2,\ldots,\sigma_p^2)\).

Not accounting for the hidden covariates, \(Z\), can reduce power and result in poor control of false discovery rate.

`succotashr` fits this model under a two-step empirical Bayesian approach. It places a non-parametric prior on \(\beta\) and jointly estimates \(\beta\) and \(Z\alpha\). The main function is `succotash`.

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
