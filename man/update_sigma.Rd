% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factor_mle.R
\name{update_sigma}
\alias{update_sigma}
\title{Update variances in regularized maximum likelihood factor analysis.}
\usage{
update_sigma(Y, A, lambda = 0.01)
}
\arguments{
\item{Y}{A matrix. This is the \code{n} by \code{p} data matrix.}

\item{A}{A matrix. This the the \code{n} by {p} low rank mean
matrix.}

\item{lambda}{A numeric. This is the tuning parameter. The MLE is
unbounded for lambda = 0.}
}
\value{
\code{sig_diga} The update for for the variances.
}
\description{
The is the update of the variances of the columns in the
regularized factor analysis.
}

