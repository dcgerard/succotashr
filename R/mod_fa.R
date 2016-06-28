#######################
## Filename: mod_fa.R
## Created by: David Gerard
## Created on: 03/08/2016
## Synopsis: Factor analysis with a prior on column variances but nothing else.
#######################

#' Update the factors given the loadings and the data matrix.
#'
#'@param L A matrix of numerics. n by k. The current loadings.
#'@param Y A matrix of numerics. The data matrix. n by p.
#'
update_f <- function(L, Y) {
    Fnew <- solve(t(L) %*% L) %*% t(L) %*% Y
    return(Fnew)
}

#' Update the loadings given the factors, the data matrix, and the
#' rate of the variances.
#'
#' @param L A matrix of numerics. n by k. The current loadings.
#' @param F A matrix of numerics. k by p. The factors.
#' @param Y A matrix of numerics. n by p. The data.
#' @param beta A positive numeric. The rate of the variances.
#'
#'
#'
update_l <- function(L, F, Y, beta) {
    optim_out <- stats::optim(par = L, fn = fn_l, gr = gr_l, F = F,
                              Y = Y, beta = beta, method = "L-BFGS-B")
    return(L = optim_out$par)
}


#' Function to minimize wrt the loadings.
#'
#' @inheritParams update_l
#'
#'
#'
#'
fn_l <- function(L, F, Y, beta) {
    L <- matrix(L, nrow = nrow(Y))
    resid_fit <- Y -  L %*% F
    sse_vec <- colSums(resid_fit ^ 2)
    return(sum(log(sse_vec + 2 * beta)))
}

#' Gradient of function to minimize wrt loadings.
#'
#'
#' @inheritParams update_l
#'
#'
#'
gr_l <- function(L, F, Y, beta) {
    L <- matrix(L, nrow = nrow(Y))
    fit <- L %*% F
    resid_fit <- Y - fit
    sse_vec <- colSums(resid_fit ^ 2)
    sum_tot <- 0
    for (index in 1:ncol(Y)) {
        sum_tot <- sum_tot + (fit[, index] %*% t(F[, index]) - Y[, index] %*% t(F[, index])) /
            (sse_vec[index] + 2 * beta)
    }
    sum_tot <- sum_tot * 2
    return(sum_tot)
}

#' Updates the hyperparameters for the gamma prior.
#'
#' @param alpha A positive numeric. The shape parameter of the variances.
#' @inheritParams update_l
#'
update_ab <- function(alpha, beta, Y, L, F) {
    optim_out <- stats::optim(par = c(alpha, beta), fn = fn_ab,
                              gr = gr_ab, method = "L-BFGS-B", Y = Y,
                              L = L, F = F,
                              control = list(fnscale = -1))
    return(optim_out$par, )
}

#' Function to optimize wrt alpha and beta
#'
#' @param ab A 2 vector of positive numerics. \code{ab[1]} is the current value
#'   of alpha, \code{ab[2]} is the current value of beta.
#' @inheritParams update_ab
#'
#'
#'
#'
fn_ab <- function(ab, Y, L, F) {
    alpha <- ab[1]
    beta <- ab[2]
    n <- nrow(Y)
    p <- ncol(Y)
    resid_fit <- Y -  L %*% F
    sse_vec <- colSums(resid_fit ^ 2) ## the "omega"s from the write-up.

    fval <- p * alpha * log(beta) + p * lgamma(n / 2 + alpha) - p * lgamma(alpha) -
        (n / 2 + alpha) * sum(log(sse_vec / 2 + beta))
    return(fval)
}

#' Gradient of function to optimize wrt alpha and beta.
#'
#' @inheritParams fn_ab
#'
#'
#'
#'
gr_ab <- function(ab, Y, L, F) {
    alpha <- ab[1]
    beta <- ab[2]
    n <- nrow(Y)
    p <- ncol(Y)
    resid_fit <- Y -  L %*% F
    sse_vec <- colSums(resid_fit ^ 2) ## the "omega"s from the write-up.

    da <- p * log(beta) + p * digamma(n / 2 + alpha) - p * digamma(alpha) -
        sum(log(sse_vec / 2 + beta))
    db <- p * alpha / beta - (n / 2 + alpha) * sum(1 / (sse_vec / 2 + beta))
    grad_val <- c(da, db)
    return(grad_val)
}

#' Moderated factor analysis. Optimize over loadings, factors, and
#' hyperparameters of variances.
#'
#' @param Y A matrix of numerics. n by p.
#' @param k A positive integer. The rank of the mean.
#' @param tol A positive numeric. The stopping criterion.
#' @param itermax A positive integer. The maximum number of iterations
#'     to run through in the optimization.
#'
#' @export
#'
#'
mod_fa <- function(Y, k, tol = 10 ^ -6, itermax =  100) {

    ## svd for initial values of L and F
    svdY <- svd(Y)
    L <- (svdY$u %*% diag(svdY$d))[, 1:k]
    F <- t(svdY$v[, 1:k])

    resid_vals <- Y - L %*% F
    mse_vec <- colSums(resid_vals ^ 2) / (nrow(Y) - 1)

    ## method of moments for starting alpha and beta based on log-gamma distributed variances.
    emp_var <- stats::var(log(mse_vec))
    emp_mean <- mean(log(mse_vec))
    alpha <- stats::uniroot(f = minus_trigamma, emp_var = emp_var, interval = c(0, 100))$root
    beta <- exp(digamma(alpha) - emp_mean)

    current_fit <- L %*% F
    err <- tol + 1
    iter_index <- 1
    while (err > tol & iter_index < itermax) {
        old_fit <- current_fit
        L <- update_l(L = L, F = F, Y = Y, beta = beta)
        F <- update_f(L = L, Y = Y)
        current_fit <- L %*% F
        err <- sum(abs(current_fit - old_fit))
        cat("Err:", err, "\n")
        iter_index <- iter_index + 1
    ## ab <- update_ab(alpha = alpha, beta = beta, Y = Y, L = L, F = F)
    ## alpha <- ab[1]
    ## beta <- ab[2]
    }

    resid_vals <- Y - L %*% F
    sse_vec <- colSums(resid_vals ^ 2)
    sigma2est <- (sse_vec + 2 * beta) / (nrow(Y) + 2 * alpha)
    return(list(F = F, sigma2est = sigma2est))
}

#' Trigamma minus a constant.
#'
#' @param x A positive numeric.
#' @param emp_var A numeric.
minus_trigamma <- function(x, emp_var) {
    return(trigamma(x) - emp_var)
}


#' Use PCA to estimate confounders, then shrink variances using Smyth's method
#' implemented in limma.
#'
#' Factors are estimated using PCA. Variances are shrunk using an empirical
#' Bayesian approach assuming a hierarchical model as described in Smyth (2004).
#'
#' @param Y A matrix of numerics. n by p.
#' @param k A positive integer. The known rank of the mean.
#' @param df A string. Should we subtract off the rank (\code{"rank_based"}) or
#'   just 1 (\code{"minus_one"})?
#'
#' @export
#'
#' @references Smyth, G. K. (2004).
#'   \href{http://www.statsci.org/smyth/pubs/ebayes.pdf}{Linear models and
#'   empirical Bayes methods for assessing differential expression in microarray
#'   experiments}. Statistical Applications in Genetics and Molecular Biology
#'   Volume 3, Issue 1, Article 3.
pca_shrinkvar <- function(Y, k, df = "rank_based") {
    svdY <- svd(Y)
    L <- (svdY$u %*% diag(svdY$d))[, 1:k]
    F <- t(svdY$v[, 1:k])

    resid_vals <- Y - L %*% F
    if (df == "rank_based") {
        df <- (nrow(Y) - k)
    } else if (df == "minus_one") {
        df <- (nrow(Y) - 1)
    }
    mse_vec <- colSums(resid_vals ^ 2) / df
    ## Need to think harder about what I should divide by.
    ## In particular, shouldn't it depend on k? With larger k needing smaller df?
    sv_out <- limma::squeezeVar(var = mse_vec, df = df)
    sig_diag <- sv_out$var.post
    df <- sv_out$df.prior
    return(list(F = F, sigma2est = sig_diag, df = df))
}



#' Basic PCA.
#'
#' Most of this code is from the package \code{cate}. I corrected some
#' problems. Specifically, I allow \code{r = 0} and I included a few
#' needed \code{drop = FALSE} terms. I also divide by \code{nrow(Y) -
#' r} rather than by \code{nrow(Y)}.
#'
#'
#' @param Y A matrix of numerics. The data.
#' @param r the rank.
#'
#' @author David Gerard
pca_naive <- function (Y, r) {
    if (r == 0) {
        Gamma <- NULL
        Z <- NULL
        Sigma <- apply(Y, 2, function(x) mean(x ^ 2))
    } else {
        svd_Y <- svd(Y)
        Gamma <- svd_Y$v[, 1:r, drop = FALSE] %*% diag(svd_Y$d[1:r], r, r) /
            sqrt(nrow(Y))
        Z <- sqrt(nrow(Y)) * svd_Y$u[, 1:r, drop = FALSE]
        Sigma <- apply(Y - Z %*% t(Gamma), 2, function(x) sum(x ^ 2)) / (nrow(Y) - r)
    }
    return(list(Gamma = Gamma, Z = Z, Sigma = Sigma))
}
