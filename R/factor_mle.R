#' Truncates small numbers to 0.
#'
#' Given an array, matrix, or vector, \code{trim} will truncate all elements
#' smaller than \code{epsilon} (in absolute value) to zero.
#'
#' All elements in \code{X} that are smaller than \code{epsilon} (in absolute
#' value) will be set to zero then returned.
#'
#' @param x An array, a matrix, or a vector.
#' @param epsilon A numeric.
#'
#' @export
#'
#' @examples
#' X <- c(0, 1, 10^-7, -1, -10^-7)
#' X
#' trim(X)
trim <- function(x, epsilon = 10 ^ -6) {
    x[abs(x) < epsilon] <- 0
    return(x)
}

#' Update variances in regularized maximum likelihood factor analysis.
#'
#' The is the update of the variances of the columns in the regularized factor
#' analysis.
#'
#'
#' @param Y A matrix. This is the \code{n} by \code{p} data matrix.
#' @param A A matrix. This the the \code{n} by {p} low rank mean matrix.
#' @param lambda A numeric. This is the tuning parameter. The MLE is unbounded
#'   for lambda = 0.
#'
#' @return \code{sig_diga} The update for for the variances.
update_sigma <- function(Y, A, lambda = 0.01) {
    eps <- Y - A
    sig_diag <- apply(eps ^ 2, 2, mean) + lambda / nrow(Y)
    return(sig_diag)
}

#' Update low rank mean in regularized maximum likelihood factor analysis.
#'
#' This is the update of the low rank mean in the regularized factor analysis.
#' The low rank mean is not regularized, so this is the exact same update as in
#' \href{http://projecteuclid.org/euclid.aoas/1356629055}{Sun et al (2012)}.
#'
#' The update is just a truncated SVD of a scaled data matrix, where the scaling
#' depends on the current estimates of the variances of the columns.
#'
#' @param Y A matrix. This is the \code{n} by \code{p} data matrix.
#' @param sig_diag A vector of length \code{p}. These are the estimates of the
#'   variances of the columns of \code{Y}.
#' @param k An integer. The rank of the mean matrix.
#'
#' @return \code{A_new} The update of the mean matrix.
update_A <- function(Y, sig_diag, k) {
    svd_a <- svd((1 / sqrt(sig_diag)) * t(Y), nu = k, nv = k)
    A_new <- t(sqrt(sig_diag) * tcrossprod(t(svd_a$d[1:k] * t(svd_a$u)), svd_a$v))
    return(A_new)
}

#' Regularized normal log-likelihood.
#'
#' \code{f_val} will return the regularized normal log-likelihood where thre is
#' a low rank mean matrix and only a diagonal covariance matrix along the
#' columns.
#'
#' @param Y A matrix. This is the \code{n} by \code{p} data matrix.
#' @param sig_diag A vector of length \code{p}. These are the estimates of the
#'   variances of the columns of \code{Y}.
#' @param A A matrix. This the the \code{n} by {p} low rank mean matrix.
#' @param lambda A numeric. The regularization parameter.
#'
#' @return \code{llike} A numeric. The regularized log-likelihood.
f_val <- function(Y, sig_diag, A, lambda = 0.01) {
    n <- nrow(Y)
    p <- ncol(Y)
    llike <- -n / 2 * sum(log(sig_diag)) -
      sum(((1 / sqrt(sig_diag)) * t(Y - A)) ^ 2) / 2 - log(2 * pi) * n * p / 2 - lambda * sum(1 / sig_diag) / 2
    return(llike)
}

#' Regularized maximum likelihood factor analysis with heteroscedastic columns.
#'
#' \code{factor_mle} implements regularized maximum likelihood estimation on a
#' data matrix where the mean is low rank (with the rank known) and the columns
#' are heteroscedastic but independent.
#'
#' This function calculates the regularized MLE under a normal model with a
#' low-rank mean, a diagonal column covariance matrix, and an identity row
#' covariance matrix. The regularization is on the covariance matrix.
#'
#' The unregularized version of this factor analysis can be found in
#' \href{http://projecteuclid.org/euclid.aoas/1356629055}{Sun et al (2012)}.
#' However, the likelihood is unbounded and there are many "MLE's" that have
#' that unbounded likelihood. This apparently was unnoticed in
#' \href{http://projecteuclid.org/euclid.aoas/1356629055}{Sun et al (2012)}.
#' \code{sig_reg} should never be set to 0.
#'
#' @param Y A matrix. This is an \code{n} by \code{p} matrix, where the the
#'   columns are heteroscedastic.
#' @param k A numeric. The rank of the mean matrix.
#' @param itermax An integer. The maximum number of block-coordinate ascent
#'   iterations to perform when calculating the MLE.
#' @param tol A numeric. When the difference from one of the ratio of two
#'   successive log-likelihoods is less than \code{tol}, the algorithm will
#'   stop.
#' @param print_diff A logical. Should we print to the screen the updates?
#' @param sig_reg The regularization parameter. Never set to 0.
#'
#' @export
#'
#' @return \code{A} The low rank mean estimate.
#'
#'   \code{sig_diag} A vector of the variance estimates of the columns.
factor_mle <- function(Y, k, itermax = 100, tol = 10 ^ -6, print_diff = FALSE, sig_reg = 0.01) {
    cat("Working on factor_mle().\n")
    p <- ncol(Y)
    ## set inital value
    sig_diag <- rep(1, p)
    A <- update_A(Y, sig_diag, k)
    f_old <- f_val(Y, sig_diag, A, sig_reg)
    err <- tol + 1
    iter <- 0
    while (iter < itermax & err > tol) {
        sig_diag <- update_sigma(Y, A, lambda = sig_reg)
        A <- update_A(Y, sig_diag, k)
        f_new <- f_val(Y, sig_diag, A, sig_reg)
        err <- abs(f_new / f_old - 1)

        if ( (f_new > f_old) | (iter < 2)) {
            sig_diag_old <- sig_diag
            A_old <- A
            f_old <- f_new
            iter <- iter + 1
        } else {
            ## if the function value goes down then it's no longer doing
            ## anything
            cat("Running into numerical issues in factor_mle.r\n")
            return(list(A = A_old, sig_diag = sig_diag_old))
        }

        if (print_diff) {
            cat("          Diff:", err, "\n")
            cat("Function Value:", f_new, "\n")
            cat("     Iteration:", iter, "\n\n")
        }
    }
    cat("Did", iter, "iterations in factor_mle().\n")

    return(list(A = A, sig_diag = sig_diag))
}
