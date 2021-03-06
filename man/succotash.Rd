% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/succotash.R
\name{succotash}
\alias{succotash}
\title{Surrogate and Confounder Correction Occuring Together with Adaptive
SHrinkage.}
\usage{
succotash(Y, X, k = NULL, sig_reg = 0.01, num_em_runs = 2,
  z_start_sd = 1, two_step = TRUE, fa_method = c("pca", "reg_mle",
  "quasi_mle", "homoPCA", "pca_shrinkvar", "mod_fa", "flash_hetero", "non_homo",
  "non_hetero", "non_shrinkvar"), lambda_type = c("zero_conc", "ones"),
  mix_type = c("normal", "uniform"), likelihood = c("normal", "t"),
  lambda0 = 10, tau_seq = NULL, em_pi_init = NULL,
  plot_new_ests = FALSE, em_itermax = 200, var_scale = TRUE,
  inflate_var = 1, optmethod = c("coord", "em"), use_ols_se = FALSE,
  z_init_type = c("null_mle", "random"), var_scale_init_type = c("null_mle",
  "one", "random"))
}
\arguments{
\item{Y}{An \code{n} by \code{p} matrix of response variables.}

\item{X}{An \code{n} by \code{q} matrix of covariates. Only the
variable in the last column is of interest.}

\item{k}{An integer. The number of hidden confounders. If
\code{NULL} and \code{sva} is installed, this will be
estimated, by the \code{num.sv} function in the \code{sva}
package available on Bioconductor.}

\item{sig_reg}{A numeric. If \code{fa_method} is \code{"reg_mle"},
then this is the value of the regularization parameter.}

\item{num_em_runs}{An integer. The number of times we should run
the EM algorithm.}

\item{z_start_sd}{A positive numeric. At the beginning of each EM
algorithm, \code{Z} is initiated with independent mean zero
normals with standard deviation \code{z_start_sd}.}

\item{two_step}{A logical. Should we run the two-step SUCCOTASH
procedure of inflating the variance (\code{TRUE}) or not
(\code{FALSE})? Defaults to \code{TRUE}.}

\item{fa_method}{Which factor analysis method should we use? The
regularized MLE implemented in \code{\link{factor_mle}}
(\code{"reg_mle"}), two methods fromthe package \code{cate}:
the quasi-MLE (\code{"quasi_mle"}) from
\href{http://projecteuclid.org/euclid.aos/1334581749}{Bai and
Li (2012)}, just naive PCA (\code{"pca"}), FLASH
(\code{"flash_hetero"}), homoscedastic PCA (\code{"homoPCA"}),
PCA followed by shrinking the variances using limma
(\code{"pca_shrinkvar"}), or moderated factor analysis
(\code{"mod_fa"}).  Three methods for no confounder adjustment
are available, \code{"non_homo"}, \code{"non_shrinkvar"}, and
\code{"non_hetero"}.}

\item{lambda_type}{See \code{\link{succotash_given_alpha}} for
options on the regularization parameter of the mixing
proportions.}

\item{mix_type}{Should the prior be a mixture of normals
\code{mix_type = 'normal'} or a mixture of uniforms
\code{mix_type = 'uniform'}?}

\item{likelihood}{Which likelihood should we use? Normal
(\code{"normal"}) or t (\code{"t"})?}

\item{lambda0}{If \code{lambda_type = "zero_conc"}, then
\code{lambda0} is the amount to penalize \code{pi0}.}

\item{tau_seq}{A vector of length \code{M} containing the standard
deviations (not variances) of the mixing distributions.}

\item{em_pi_init}{A vector of length \code{M} containing the
starting values of \eqn{\pi}. If \code{NULL}, then one of three
options are implemented in calculating \code{pi_init} based on
the value of \code{pi_init_type}. Only available in normal
mixtures for now.}

\item{plot_new_ests}{A logical. Should we plot the mixing
proportions at each iteration of the EM algorithm?}

\item{em_itermax}{A positive numeric. The maximum number of
iterations to run during the EM algorithm.}

\item{var_scale}{A logical. Should we update the scaling on the
variances (\code{TRUE}) or not (\code{FALSE}). Only works for
the normal mixtures case right now. Defaults to \code{TRUE}.}

\item{inflate_var}{A positive numeric. The multiplicative amount to
inflate the variance estimates by. There is no theoretical
justification for it to be anything but 1, but I have it in
here to play around with it.}

\item{optmethod}{Either coordinate ascent (\code{"coord"}) or an EM
algorithm (\code{"em"}). Coordinate ascent is currently only
implemented in the uniform mixtures case, for which it is the
default.}

\item{use_ols_se}{A logical. Should we use the standard formulas
for OLS of X on Y to get the estimates of the variances
(\code{TRUE}) or not (\code{FALSE})}

\item{z_init_type}{How should we initiate the confounders? At the
all-null MLE (\code{"null_mle"}) or from iid standard normals
(\code{"random"})?}

\item{var_scale_init_type}{If \code{var_scale = TRUE}, how should
we initiate the variance inflaiton parameter? From the all-null
MLE (\code{"null_mle"}), at no inflation (\code{"one"}), or
from a chi-squared distribution with one degree of freedom
(\code{"random"})?}
}
\value{
See \code{\link{succotash_given_alpha}} for details of
    output.

  \code{Y1_scaled} The OLS estimates.

  \code{sig_diag_scaled} The estimated standard errors of the
  estimated effects (calculated from the factor analysis step)
  times \code{scale_val}.

  \code{sig_diag} The estimates of the gene-wise variances (but not
  times \code{scale_val}).

  \code{pi0} A non-negative numeric. The marginal probability of
  zero.

  \code{alpha_scaled} The scaled version of the estimated
  coefficients of the hidden confounders.

  \code{Z} A vector of numerics. Estimated rotated confounder in
  second step of succotash.

  \code{pi_vals} A vector of numerics between 0 and 1. The mixing
  proportions.

  \code{tau_seq} A vector of non-negative numerics. The mixing
  standard deviations (not variances).

  \code{lfdr} A vector of numerics between 0 and 1. The local false
  discovery rate. I.e. the posterior probability of a coefficient
  being zero.

  \code{lfsr} A vector of numerics between 0 and 1. The local false
  sign rate. I.e. the posterior probability of making a sign error
  if one chose the most probable sign.

  \code{qvals} A vector of numerics between 0 and 1. The
  q-values. The average error rate if we reject all hypotheses that
  have smaller q-value.

  \code{betahat} A vector of numerics. The posterior mean of the
  coefficients.
}
\description{
This function implements the full SUCCOTASH method. First, it
rotates the response and explanatory variables into a part that we
use to estimate the confounding variables and the variances, and a
part that we use to estimate the coefficients of the observed
covariates. This function will implement a factor analysis for the
first part then run \code{\link{succotash_given_alpha}} for the
second part.
}
\details{
The assumed mode is \deqn{Y = X\beta + Z\alpha + E.} \eqn{Y} is a
\eqn{n} by \code{p} matrix of response varaibles. For example, each
row might be an array of log-transformed and quantile normalized
gene-expression data.  \eqn{X} is a \eqn{n} by \eqn{q} matrix of
observed covariates. It is assumed that all but the last column of
which contains nuisance parameters. For example, the first column
might be a vector of ones to include an intercept.  \eqn{\beta} is
a \eqn{q} by \eqn{p} matrix of corresponding coefficients.  \eqn{Z}
is a \eqn{n} by \eqn{k} matrix of confounder
variables. \eqn{\alpha} is the corresponding \eqn{k} by \eqn{p}
matrix of coefficients for the unobserved confounders. \eqn{E} is a
\eqn{n} by \eqn{p} matrix of error terms. \eqn{E} is assumed to be
matrix normal with identity row covariance and diagonal column
covariance \eqn{\Sigma}. That is, the columns are heteroscedastic
while the rows are homoscedastic independent.

This function will first rotate \eqn{Y} and \eqn{X} using the QR
decomposition. This separates the model into three parts. The first
part only contains nuisance parameters, the second part contains
the coefficients of interest, and the third part contains the
confounders. \code{succotash} applies a factor analysis to the
third part to estimate the confounding factors, then runs an EM
algorithm on the second part to estimate the coefficients of
interest.

Many forms of factor analyses are avaiable. The default is PCA with
the column-wise residual mean-squares as the estimates of the
column-wise variances.
}
\seealso{
\code{\link{succotash_given_alpha}}, \code{\link{factor_mle}},
  \code{\link{succotash_summaries}}.
}

