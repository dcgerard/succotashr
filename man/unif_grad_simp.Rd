% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniform_coordinate.R
\name{unif_grad_simp}
\alias{unif_grad_simp}
\title{This is me being very careful about calculating the gradient of the
likelihood function wrt Z.}
\usage{
unif_grad_simp(pi_Z, lambda, alpha, Y, sig_diag, left_seq = NULL,
  right_seq = NULL, a_seq = NULL, b_seq = NULL, var_scale = TRUE,
  likelihood = c("normal", "t"), df = NULL)
}
\arguments{
\item{pi_Z}{A vector. The first \code{M} values are the current
values of \eqn{\pi}. The last \code{k} values are the current
values of \eqn{Z}.}

\item{lambda}{A vector. This is a length \code{M} vector with the
regularization parameters for the mixing proportions.}

\item{alpha}{A matrix. This is of dimension \code{p} by \code{k} and
are the coefficients to the confounding variables.}

\item{Y}{A matrix of dimension \code{p} by \code{1}. These are the
observed regression coefficients of the observed variables.}

\item{sig_diag}{A vector of length \code{p} containing the variances
of the observations.}

\item{left_seq}{The left endpoints of the uniforms.}

\item{right_seq}{The right endpoints of the uniforms}

\item{a_seq}{A vector of negative numerics containing the left
endpoints of the mixing uniforms.}

\item{b_seq}{A vector of positiv numerics containing the right
endpoints of the mixing uniforms.}

\item{var_scale}{A logical. Should we update the scaling on the
variances (\code{TRUE}) or not (\code{FALSE}).}

\item{likelihood}{Can be \code{"normal"} or \code{"t"}.}

\item{df}{A positive numeric. The degrees of freedom if the
likelihood is t.}
}
\value{
The gradient for Z.
}
\description{
This is me being very careful about calculating the gradient of the
likelihood function wrt Z.
}

