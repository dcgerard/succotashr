#'\code{succotashr}: An \code{R} package for running Surrogate and Confounder
#'Correction Occuring Together with Adaptive SHrinkage.
#'
#'This package contains functions for implementing the SUCCOTASH method in
#'\code{R}. This is method to account for hiddent confounders when performing
#'linear regression. The important functions in \code{succotashr} are
#'\code{\link{succotash}}, \code{\link{succotash_given_alpha}}, and
#'\code{\link{factor_mle}}.
#'
#'@section \code{succotashr} functions:
#'
#'  \code{\link{draw_beta}}: Draw from a mixture of normals.
#'
#'  \code{\link{factor_mle}}: Regularized maximum likelihood factor analysis
#'  with heteroscedastic columns.
#'
#'  \code{\link{f_val}}: Regularized normal log-likelihood.
#'
#'  \code{\link{lfdr_to_q}}: Transform local false discovery rates to q-values.
#'
#'  \code{\link{succotash}}: Surrogate and Confounder Correction Occuring
#'  Together with Adaptive SHrinkage.
#'
#'  \code{\link{succotash_em}}: An EM algorithm for maximizing the SUCCOTASH
#'  log-likelihood.
#'
#'  \code{\link{succotash_fixed}}: A fixed-point iteration of the EM algorithm.
#'
#'  \code{\link{succotash_given_alpha}}: Maximize the SUCCOTASH log-likelihood
#'  and return posterior summaries.
#'
#'  \code{\link{succotash_llike}}: The SUCCOTASH log-likelihood.
#'
#'  \code{\link{succotash_summaries}}: Provides posterior summaries in the
#'  SUCCOTASH model.
#'
#'  \code{\link{update_A}}: Update low rank mean in regularized maximum
#'  likelihood factor analysis.
#'
#'  \code{\link{update_sigma}}: Update variances in regularized maximum
#'  likelihood factor analysis.
#'
#'  \code{\link{trim}}: Truncates small numbers to 0.
#'
#'@docType package
#'@name succotashr
NULL
#> NULL
