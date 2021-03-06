% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/succotash.R
\name{succotash_llike}
\alias{succotash_llike}
\title{The SUCCOTASH log-likelihood.}
\usage{
succotash_llike(pi_Z, lambda, alpha, Y, tau_seq, sig_diag,
  plot_new_ests = FALSE, var_scale = TRUE)
}
\arguments{
\item{pi_Z}{A vector. The first \code{M} values are the current
values of \eqn{\pi}. The last \code{k} values are the current
values of \eqn{Z}.}

\item{lambda}{A vector. This is a length \code{M} vector with the
regularization parameters for the mixing proportions. If
\code{NULL} then refer to \code{lambda_type}.}

\item{alpha}{A matrix. This is of dimension \code{p} by \code{k}
and are the coefficients to the confounding variables.}

\item{Y}{A matrix of dimension \code{p} by \code{1}. These are the
observed regression coefficients of the observed variables.}

\item{tau_seq}{A vector of length \code{M} containing the standard
deviations (not variances) of the mixing distributions.}

\item{sig_diag}{A vector of length \code{p} containing the
variances of the observations.}

\item{plot_new_ests}{A logical. Should we plot the new estimates of
pi?}

\item{var_scale}{A logical. Should we update the scaling on the
variances (\code{TRUE}) or not (\code{FALSE}).}
}
\value{
\code{llike_new} A numeric. The value of the SUCCOTASH
  log-likelihood.
}
\description{
\code{succotash_llike} returns the SUCCOTASH log-likelihood. Returning the
regularized log-likelihood is currently not implemented.
}

