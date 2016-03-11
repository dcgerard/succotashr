#' A fixed-point iteration of the EM algorithm.
#'
#' This is a fixed-point iteration for the SUCCOTASH EM algorithm. This updates
#' the estimte of the prior and the estimate of the hidden covariates.
#'
#' @param pi_Z A vector. The first \code{M} values are the current values of
#'   \eqn{\pi}. The last \code{k} values are the current values of \eqn{Z}.
#' @inheritParams succotash_em
#' @inheritParams succotash_given_alpha
#'
#' @return \code{pi_new} A vector of length \code{M}. The update for the mixing
#'   components.
#'
#'   \code{Z_new} A vector of length \code{k}. The update for the confounder
#'   covariates.
succotash_fixed <- function(pi_Z, lambda, alpha, Y, tau_seq, sig_diag, plot_new_ests = FALSE) {
    M <- length(tau_seq)
    p <- nrow(Y)
    k <- length(pi_Z) - M
    pi_old <- pi_Z[1:M]
    if (k != 0) {
        Z_old <- matrix(pi_Z[(M + 1):(M + k)], nrow = k)
    }

    if (!is.null(alpha)) {
        mean_mat <- matrix(rep(alpha %*% Z_old, M), ncol = M, nrow = p)
    } else {
        mean_mat <- matrix(0, ncol = M, nrow = p)
    }
    var_mat <- outer(sig_diag, tau_seq ^ 2, "+")
    top_vals <- t(pi_old * t(dnorm(matrix(rep(Y, M), ncol = M, nrow = p), mean = mean_mat,
                                   sd = sqrt(var_mat))))
    T <- t(1 / rowSums(top_vals) * top_vals)

    Theta_diag <- colSums((T / t(var_mat)) / 2)

    T_sum <- rowSums(T)

    if (!is.null(alpha)) {
        Z_new <- solve(t(alpha) %*% (Theta_diag * alpha)) %*% t(Theta_diag * alpha) %*% Y
    } else if (k != 0) {
        Z_new <- rep(NA, length = k)
    } else {
        Z_new <- NULL
    }
    pi_new <- (T_sum + lambda - 1) / (p - M + sum(lambda))

    if(plot_new_ests) {
        plot(tau_seq, pi_new, type = "h", xlab = expression(tau[k]), ylab = expression(pi[k]))
    }

    return(c(pi_new, Z_new))
}

#' The SUCCOTASH log-likelihood.
#'
#' \code{succotash_llike} returns the SUCCOTASH log-likelihood. Returning the
#' regularized log-likelihood is currently not implemented.
#'
#' @param pi_Z A vector. The first \code{M} values are the current values of
#'   \eqn{\pi}. The last \code{k} values are the current values of \eqn{Z}.
#' @param lambda A vector. This is a length \code{M} vector with the
#'   regularization parameters for the mixing proportions.
#' @param alpha A matrix. This is of dimension \code{p} by \code{k} and are the
#'   coefficients to the confounding variables.
#' @param Y A matrix of dimension \code{p} by \code{1}. These are the observed
#'   regression coefficients of the observed variables.
#' @param tau_seq A vector of length \code{M} containing the standard deviations
#'   (not variances) of the mixing distributions.
#' @param sig_diag A vector of length \code{p} containing the variances of the
#'   observations.
#' @param plot_new_ests A logical. Needed to be here so can use SQUAREM package
#'   without an error because being used in \code{succotash_fixed}.
#'
#' @export
#'
#' @return \code{llike_new} A numeric. The value of the SUCCOTASH
#'   log-likelihood.
succotash_llike <- function(pi_Z, lambda, alpha, Y, tau_seq, sig_diag, plot_new_ests = FALSE) {
    M <- length(tau_seq)
    p <- nrow(Y)
    k <- length(pi_Z) - M
    pi_current <- pi_Z[1:M]
    if (k != 0) {
        Z_current <- matrix(pi_Z[(M + 1):(M + k)], nrow = k)
    }

    if (!is.null(alpha)) {
        mean_mat <- matrix(rep(alpha %*% Z_current, M), ncol = M, nrow = p)
    } else {
        mean_mat <- matrix(0, ncol = M, nrow = p)
    }
    var_mat <- outer(sig_diag, tau_seq ^ 2, "+")
    top_vals <- t(pi_current * t(dnorm(matrix(rep(Y, M), ncol = M, nrow = p), mean = mean_mat,
                                       sd = sqrt(var_mat))))

    llike_new <- sum(log(rowSums(top_vals)))
    return(llike_new)
}

#' An EM algorithm for maximizing the SUCCOTASH log-likelihood.
#'
#' Let \eqn{Y} (\eqn{p} by \eqn{1}) be multivariate normal with mean
#' \eqn{\beta + \alpha Z} and diagonal covariance \eqn{\Sigma}, where
#' \eqn{\alpha} and \eqn{\Sigma} are both known. If \eqn{\beta} is
#' assumed to be a mixture of normals with known variances and unknown
#' mixing proportions \eqn{\pi} (\eqn{p} by \eqn{1}), then this
#' function will maximize the marginal likelihood over \eqn{Z} and
#' \eqn{\pi}.
#'
#' This function uses the \code{SQUAREM} package with fixed point
#' iteration \code{succotash_fixed} to run the EM. There can be a lot
#' of local modes, so this function should be run at many starting
#' locations.
#'
#' @param pi_init A vector of length \code{M} containing the starting
#'     values of \eqn{\pi}. If \code{NULL}, then one of three options
#'     are implemented in calculating \code{pi_init} based on the
#'     value of \code{pi_init_type}.
#' @param Z_init A \code{k} by \code{1} matrix. These are the initial
#'     values of the unobserved covariates. If its value is
#'     \code{NULL}, then each element of \code{Z_init} will be drawn
#'     from a mean zero normal with standard deviation
#'     \code{z_start_sd}.
#' @param itermax An integer. The maximum number of fixed-point
#'     iterations to run the EM algorithm.
#' @param tol A numeric. The stopping criterion is the absolute
#'     difference of the ratio of subsequent iterations'
#'     log-likelihoods from 1.
#' @param z_start_sd A positive numeric. If \code{Z_init} is
#'     \code{NULL}, then the starting values for \eqn{Z} are drawn
#'     from a mean zero normal with standard devation
#'     \code{z_start_sd}.
#' @param print_note Should we print that we're doing an EM?
#' @param pi_init_type How should we choose the initial values of
#'     \eqn{\pi}.  Possible values of \code{"random"},
#'     \code{"uniform"}, and \code{"zero_conc"}. If \code{"random"}
#'     then the initial values of \eqn{\pi} are drawn uniformly over
#'     the probability simplex. If \code{"uniform"}, then each element
#'     of \code{pi_init} is given mass \code{1 / M}. If
#'     \code{"zero_conc"} then the last \code{M - 1} elements of
#'     \code{pi_init} are given mass \code{1 / p} and
#'     \code{pi_init[1]} is given mass \code{1 - sum(pi_init[2:M])}.
#' @inheritParams succotash_given_alpha
#'
#' @export
#'
#' @return \code{pi_vals} A vector of length \code{M}. The estimates
#'     of the mixing proportions.
#'
#'   \code{Z} A matrix of dimension \code{k} by \code{1}. The
#'   estimates of the confounder covariates.
#'
#'   \code{llike} A numeric. The final value of the SUCCOTASH
#'   log-likelihood.
#'
#'   \code{tau_seq} A vector of length \code{M}. The variances of the
#'   mixing distribution.
succotash_em <- function(Y, alpha, sig_diag, tau_seq = NULL, pi_init = NULL, lambda = NULL,
                         Z_init = NULL, itermax = 1500, tol = 10 ^ -6, z_start_sd = 1,
                         print_note = FALSE, pi_init_type = "random", lambda_type = "zero_conc",
                         lambda0 = 10, plot_new_ests = FALSE) {
    if (print_note) {
        cat("Working on EM.\n")
    }
    p <- nrow(Y)
    if (!is.null(alpha)) {
        k <- ncol(alpha)
    } else {
        k <- 0
    }

    if (is.null(tau_seq)) {
        ## default grid to be same as in ASH
        tau_min <- min(sig_diag) / 10
        tau_max <- 2 * sqrt(max(Y ^ 2 - sig_diag))
        if (tau_max < 0) {
            tau_max <- 8 * tau_min
        }
        tau_current <- tau_min
        tau_seq <- c(0, tau_current)
        mult_fact <- sqrt(2)
        while (tau_current <= tau_max) {
            tau_current <- tau_current * mult_fact
            tau_seq <- c(tau_seq, tau_current)
        }
    }

    M <- length(tau_seq)  # number of classes in mixture, = K+1 in paper

    if (is.null(pi_init) | length(pi_init) != M) {
        if (pi_init_type == "random") {
            ## random start points
            pi_init <- rep(NA, length = M)
            temp <- abs(rnorm(M))
            pi_init[1:M] <- temp / sum(temp)
        } else if (pi_init_type == "uniform") {
            ## uniform mass
            pi_init[1:M] <- 1 / M
        } else if (pi_init_type == "zero_conc") {
            ## most mass at 0
            pi_init[2:M] <- min(1 / p, 1 / M)
            pi_init[1] <- 1 - sum(pi_init[2:M])
        } else {
            warning("pi_init_type is bad")
        }
    }

    if (is.null(lambda) | length(lambda) != M) {
        if (lambda_type == "zero_conc") {
            ## use same values as in ASH
            lambda <- rep(NA, length = M)
            lambda[1] <- lambda0
            lambda[2:M] <- 1
        } else if (lambda_type == "ones") {
            lambda <- rep(1, length = M)
        } else {
            warning("lambda not specified correctly")
        }
    }

    if (is.null(Z_init)) {
        Z_init <- matrix(rnorm(k, sd = z_start_sd), nrow = k)
    } else if (length(Z_init) != k) {
        Z_init <- matrix(rnorm(k, sd = z_start_sd), nrow = k)
    }

    pi_Z <- c(pi_init, Z_init)

    sq_out <- SQUAREM::fpiter(par = pi_Z, lambda = lambda,
                              alpha = alpha, Y = Y, tau_seq = tau_seq,
                              sig_diag = sig_diag,
                              plot_new_ests = plot_new_ests,
                              fixptfn = succotash_fixed,
                              objfn = succotash_llike,
                              control = list(maxiter = itermax, tol = tol))
                              ## from the SQUAREM package

    llike <- succotash_llike(pi_Z = sq_out$par, lambda = lambda, alpha = alpha, Y = Y,
                             tau_seq = tau_seq, sig_diag = sig_diag)

    pi_vals <- sq_out$par[1:M]
    if (k != 0) {
        Z <- matrix(sq_out$par[(M + 1):(M + k)], nrow = k)
    } else {
        Z <- NULL
    }

    return(list(pi_vals = pi_vals, Z = Z, llike = llike, tau_seq = tau_seq))
}

#' Maximize the SUCCOTASH log-likelihood and return posterior
#' summaries.
#'
#' This function runs \code{\link{succotash_em}} repetitively, keeping
#' the highest local mode. It then returns posterior summaries.
#'
#' Let \eqn{Y} (\eqn{p} by \eqn{1}) be multivariate normal with mean
#' \eqn{\beta + \alpha Z} and diagonal covariance \eqn{\Sigma}, where
#' \eqn{\alpha} and \eqn{\Sigma} are both known. If \eqn{\beta} is
#' assumed to be a mixture of normals with known variances and unknown
#' mixing proportions \eqn{\pi} (\eqn{p} by \eqn{1}), then this
#' function will maximize the likelihood over \eqn{Z} and
#' \eqn{\pi}. It does this by running the EM algorithm implemented in
#' \code{\link{succotash_em}} many times at different starting points.
#'
#' The defaults are to run the first EM algorithm using the
#' \code{"zero_conc"} option for \code{pi_init_type}, then use the
#' \code{"random"} option for every other EM run.
#'
#' @param Y A matrix of dimension \code{p} by \code{1}. These are the
#'     observed regression coefficients of the observed variables.
#' @param alpha A matrix. This is of dimension \code{p} by \code{k}
#'     and are the coefficients to the confounding variables.
#' @param sig_diag A vector of length \code{p} containing the
#'     variances of the observations.
#' @param num_em_runs How many times should we run the EM algorithm?
#' @param print_steps A logical. Should we write the updates after
#'     each EM algorithm?
#' @param tau_seq A vector of length \code{M} containing the standard
#'     deviations (not variances) of the mixing distributions.
#' @param em_pi_init A vector of length \code{M} containing the
#'     starting values of \eqn{\pi}. If \code{NULL}, then one of three
#'     options are implemented in calculating \code{pi_init} based on
#'     the value of \code{pi_init_type}.
#' @param lambda A vector. This is a length \code{M} vector with the
#'     regularization parameters for the mixing proportions. If
#'     \code{NULL} then refer to \code{lambda_type}.
#' @param em_Z_init A \code{k} by \code{1} matrix. These are the
#'     initial values of the unobserved covariates. If its value is
#'     \code{NULL}, then each element of \code{Z_init} will be drawn
#'     from a mean zero normal with standard deviation
#'     \code{z_start_sd}.
#' @param em_itermax An integer. The maximum number of fixed-point
#'     iterations to run the EM algorithm.
#' @param em_tol A numeric. The stopping criterion is the absolute
#'     difference of the ratio of subsequent iterations'
#'     log-likelihoods from 1.
#' @param em_z_start_sd A positive numeric. If \code{Z_init} is
#'     \code{NULL}, then the starting values for \eqn{Z} are drawn
#'     from a mean zero normal with standard devation
#'     \code{z_start_sd}.
#' @param lambda_type If \code{lambda} is \code{NULL}, then how should
#'     we choose the regularization parameters. Two options are
#'     available. If \code{lambda_type} is \code{"zero_conc"}, then
#'     \code{lambda[1] = 10} and \code{lambda[2:M] = 1}. If
#'     \code{lambda_type} is \code{"ones"} then \code{lambda = 1}.
#' @param lambda0 If \code{lambda_type = "zero_conc"}, then
#'     \code{lambda0} is the amount to penalize \code{pi0}.
#' @param plot_new_ests A logical. Should we plot the new estimates of
#'     pi?
#'
#' @return  \code{Z} A matrix  of dimension \code{k} by  \code{1}. The
#'     estimates of the confounder covariates.
#'
#'   \code{pi_vals} A vector of length \code{M}. The estimates of the
#'   mixing proportions.
#'
#'   \code{tau_seq} A vector of length \code{M}. The variances of the
#'   mixing distribution.
#'
#'   \code{lfdr} (local false discovery rate) A vector of length
#'   \code{p}. The posterior probability that \eqn{\beta_j = 0}.
#'
#'   \code{lfsr} (local false sign rate) A vector of length
#'   \code{p}. The posterior probability of making a sign error.
#'
#'   \code{qvals} A vector of length \code{p}. The q-values.
#'
#'   \code{betahat} A vector of length \code{p}. The posterior
#'   estimates of \eqn{\beta}.
#'
#' @export
#'
#' @seealso \code{\link{succotash_em}},
#'     \code{\link{succotash_summaries}}.
succotash_given_alpha <- function(Y, alpha, sig_diag, num_em_runs = 2, print_steps = FALSE,
                                  tau_seq = NULL, em_pi_init = NULL, lambda = NULL,
                                  em_Z_init = NULL, em_itermax = 1500, em_tol = 10 ^ -6,
                                  em_z_start_sd = 1,
                                  lambda_type = "zero_conc", lambda0 = 10,
                                  plot_new_ests = FALSE) {

    ## @param em_pi_init_type How should we choose the initial values of
    ## \eqn{\pi}.  Possible values of \code{"random"},
    ## \code{"uniform"}, and \code{"zero_conc"}. If \code{"random"}
    ## then the initial values of \eqn{\pi} are drawn uniformly over
    ## the probability simplex. If \code{"uniform"}, then each element
    ## of \code{pi_init} is given mass \code{1 / M}. If
    ## \code{"zero_conc"} then the last \code{M - 1} elements of
    ## \code{pi_init} are given mass \code{1 / p} and
    ## \code{pi_init[1]} is given mass \code{1 - sum(pi_init[2:M])}.

    em_out <- succotash_em(Y = Y, alpha = alpha, sig_diag = sig_diag,
                           tau_seq = tau_seq, pi_init = em_pi_init,
                           lambda = lambda, Z_init = em_Z_init,
                           itermax = em_itermax, tol = em_tol,
                           z_start_sd = em_z_start_sd,
                           pi_init_type = "zero_conc",
                           lambda_type = lambda_type,
                           lambda0 = lambda0,
                           plot_new_ests = plot_new_ests)

    if(num_em_runs > 1) {
        for (index in 2:num_em_runs) {
            em_new <- succotash_em(Y = Y, alpha = alpha,
                                   sig_diag = sig_diag, tau_seq = tau_seq,
                                   pi_init = em_pi_init, lambda = lambda,
                                   Z_init = em_Z_init,
                                   itermax = em_itermax, tol = em_tol,
                                   z_start_sd = em_z_start_sd,
                                   pi_init_type = "random",
                                   lambda_type = lambda_type,
                                   lambda0 = lambda0,
                                   plot_new_ests = plot_new_ests)
            pi_diff <- sum(abs(em_new$pi_vals - em_out$pi_vals))
            z_diff <- sum(abs(em_new$Z - em_out$Z))
            if (em_out$llike < em_new$llike) {
                em_out <- em_new
            }
            if (print_steps) {
                cat("   Rep:", index, "\n")
                cat(" lbest:", em_out$llike, "\n")
                cat(" llike:", em_new$llike, "\n")
                cat("pidiff:", pi_diff, "\n")
                cat(" zdiff:", z_diff, "\n\n")
            }
        }
    }
    
    sum_out <- succotash_summaries(Y = Y, Z = em_out$Z, pi_vals = em_out$pi_vals,
                                   alpha = alpha, sig_diag = sig_diag, tau_seq = em_out$tau_seq)
    q_vals <- lfdr_to_q(lfdr = sum_out$lfdr)
    
    return(list(Z = em_out$Z, pi_vals = em_out$pi_vals, tau_seq = em_out$tau_seq,
                lfdr = sum_out$lfdr, lfsr = sum_out$lfsr, qvals = q_vals,
                betahat = sum_out$beta_hat))
}


#' Surrogate and Confounder Correction Occuring Together with Adaptive
#' SHrinkage.
#'
#' This function implements the full SUCCOTASH method. First, it rotates the
#' response and explanatory variables into a part that we use to estimate the
#' confounding variables and the variances, and a part that we use to estimate
#' the coefficients of the observed covariates. This function will implement a
#' factor analysis for the first part then run
#' \code{\link{succotash_given_alpha}} for the second part.
#'
#' The assumed mode is \deqn{Y = X\beta + Z\alpha + E.} \eqn{Y} is a \eqn{n} by
#' \code{p} matrix of response varaibles. For example, each row might be an
#' array of log-transformed and quantile normalized gene-expression data.
#' \eqn{X} is a \eqn{n} by \eqn{q} matrix of observed covariates. It is assumed
#' that all but the last column of which contains nuisance parameters. For
#' example, the first column might be a vector of ones to include an intercept.
#' \eqn{\beta} is a \eqn{q} by \eqn{p} matrix of corresponding coefficients.
#' \eqn{Z} is a \eqn{n} by \eqn{k} matrix of confounder variables. \eqn{\alpha}
#' is the corresponding \eqn{k} by \eqn{p} matrix of coefficients for the
#' unobserved confounders. \eqn{E} is a \eqn{n} by \eqn{p} matrix of error
#' terms. \eqn{E} is assumed to be matrix normal with identity row covariance
#' and diagonal column covariance \eqn{\Sigma}. That is, the columns are
#' heteroscedastic while the rows are homoscedastic independent.
#'
#' This function will first rotate \eqn{Y} and \eqn{X} using the QR
#' decomposition. This separates the model into three parts. The first part only
#' contains nuisance parameters, the second part contains the coefficients of
#' interest, and the third part contains the confounders. \code{succotash}
#' applies a factor analysis to the third part to estimate the confounding
#' factors, then runs an EM algorithm on the second part to estimate the
#' coefficients of interest.
#'
#' The possible forms of factor analysis are a regularized maximum likelihood
#' estimator, or a quasi-mle implemented in the package \code{cate}.
#'
#' @param Y An \code{n} by \code{p} matrix of response variables.
#' @param X An \code{n} by \code{q} matrix of covariates. Only the
#'     variable in the last column is of interest.
#' @param k An integer. The number of hidden confounders. This can be
#'     estimated, for example by the \code{num.sv} function in the
#'     \code{sva} package available on Bioconductor.
#' @param sig_reg A numeric. If \code{fa_method} is \code{"reg_mle"},
#'     then this is the value of the regularization parameter.
#' @param num_em_runs An integer. The number of times we should run
#'     the EM algorithm.
#' @param z_start_sd A positive numeric. At the beginning of each EM
#'     algorithm, \code{Z} is initiated with independent mean zero
#'     normals with standard deviation \code{z_start_sd}.
#' @param fa_method Which factor analysis method should we use? The
#'     regularized MLE implemented in \code{\link{factor_mle}}
#'     (\code{"reg_mle"}), two methods fromthe package \code{cate}:
#'     the quasi-MLE (\code{"quasi_mle"}) from
#'     \href{http://projecteuclid.org/euclid.aos/1334581749}{Bai and
#'     Li (2012)}, just naive PCA (\code{"pca"}), FLASH
#'     (\code{"flash"}), homoscedastic PCA (\code{"homoPCA"}), PCA
#'     followed by shrinking the variances using limma
#'     (\code{"pca_shrinkvar"}), or moderated factor analysis
#'     (\code{"mod_fa"}).
#' @param lambda_type See \code{\link{succotash_given_alpha}} for
#'     options on the regularization parameter of the mixing
#'     proportions.
#' @param mix_type Should the prior be a mixture of normals
#'     \code{mix_type = 'normal'} or a mixture of uniforms
#'     \code{mix_type = 'uniform'}?
#' @param lambda0 If \code{lambda_type = "zero_conc"}, then
#'     \code{lambda0} is the amount to penalize \code{pi0}.
#' @param tau_seq A vector of length \code{M} containing the standard
#'     deviations (not variances) of the mixing distributions.
#' @param em_pi_init A vector of length \code{M} containing the
#'     starting values of \eqn{\pi}. If \code{NULL}, then one of three
#'     options are implemented in calculating \code{pi_init} based on
#'     the value of \code{pi_init_type}. Only available in normal
#'     mixtures for now.
#' @param likelihood Which likelihood should we use? Normal
#'     (\code{"normal"}) or t (\code{"t"})?
#'
#' @return See \code{\link{succotash_given_alpha}} for details of output.
#'
#' @export
#'
#' @seealso \code{\link{succotash_given_alpha}}, \code{\link{factor_mle}},
#'   \code{\link{succotash_summaries}}.
succotash <- function(Y, X, k, sig_reg = 0.01, num_em_runs = 2,
                      z_start_sd = 1,
                      fa_method = c("reg_mle", "quasi_mle", "pca", "flash",
                                    "homoPCA", "pca_shrinkvar", "mod_fa"),
                      lambda_type = "zero_conc", mix_type = 'normal',
                      likelihood = c("normal", "t"), lambda0 = 10,
                      tau_seq = NULL, em_pi_init = NULL) {
    ncol_x <- ncol(X)

    fa_method <- match.arg(fa_method, c("reg_mle", "quasi_mle", "pca", "flash",
                                        "homoPCA", "pca_shrinkvar", "mod_fa"))

    likelihood <- match.arg(likelihood, c("normal", "t"))

    qr_x <- qr(X)
    ## multiply by sign so that it matches with beta_hat_ols
    Q <- qr.Q(qr_x, complete = TRUE) * sign(qr.R(qr_x)[ncol_x, ncol_x])
    Y_tilde <- crossprod(Q, Y)[ncol_x:nrow(Y), ]  # discard first q-1 rows.


    n <- nrow(Y_tilde)


    ## Factor Analysis Methods -----------------------------------------------
    if (fa_method == "quasi_mle" & requireNamespace("cate", quietly = TRUE)) {
        ## then use the maximum quasi-likelihood method of wang et al
        qml_out <- cate::fa.em(Y = Y_tilde[2:n, ], r = k)
        alpha <- qml_out$Gamma
        sig_diag <- qml_out$Sigma
        nu <- n - 1
    } else if (fa_method == "reg_mle") {
        mle_out <- factor_mle(Y_tilde[2:n, ], k = k, sig_reg = sig_reg, itermax = 1000)
        alpha <- svd(mle_out$A, nv = k, nu = 0)$v
        sig_diag <- mle_out$sig_diag
        nu <- n - 1
    } else if (fa_method == "pca" & requireNamespace("cate", quietly = TRUE)) {
        pca_out <- cate::fa.pc(Y = Y_tilde[2:n, ], r = k)
        alpha <- pca_out$Gamma
        sig_diag <- pca_out$Sigma
        nu <- n - 1
    } else if (fa_method == "flash" & requireNamespace("flash", quietly = TRUE)) {
        Y_current <- Y_tilde[2:n, ]
        L_mat <- matrix(NA, nrow = nrow(Y_current), ncol = k)
        F_mat <- matrix(NA, nrow = ncol(Y_current), ncol = k)
        var_vec <- rep(NA, length = k)
        for(conf_index in 1:k) {
            flash_out <- flash::flash(Y_current)
            Y_current <- Y_current - flash_out$l %*% t(flash_out$f)
            L_mat[, conf_index] <- flash_out$l
            F_mat[, conf_index] <- flash_out$f
            var_vec[conf_index] <- flash_out$sigmae2
        }
        sig_diag <- rep(mean(var_vec), length = nrow(F_mat))
        alpha <- F_mat
        nu <- n - 1
    } else if (fa_method == "homoPCA") {
        Y_current <- Y_tilde[2:n, ]
        svdY <- svd(Y_current)
        alpha <- svdY$v[, 1:k] %*% diag(svdY$d[1:k], k, k) / sqrt(nrow(Y_current))
        Ztemp <- sqrt(nrow(Y_current)) * svdY$u[, 1:k]
        sig_diag <- rep(sum((Y_current - Ztemp %*% t(alpha)) ^ 2) /
                        (max(dim(Y_current)) * (min(dim(Y_current)) - k)), ncol(Y_current))
        nu <- n - 1
    } else if (fa_method == "pca_shrinkvar" & requireNamespace("limma", quietly = TRUE)) {
        Y_current <- Y_tilde[2:n, ]
        pca_shrinkvar_out <- pca_shrinkvar(Y_current, k, df = "minus_one")
        alpha <- t(pca_shrinkvar_out$F)
        sig_diag <- pca_shrinkvar_out$sigma2est
        nu <- pca_shrinkvar_out$df
    } else if (fa_method == "mod_fa") {
        Y_current <- Y_tilde[2:n, ]
        mod_fa_out <- mod_fa(Y_current, k)
        alpha <- t(mod_fa_out$F)
        sig_diag <- mod_fa_out$sigma2est
        nu <- n - 1
    }

    ## absorb fnorm(X) into Y_tilde[1,], alpha, and sig_diag -----------------
    fnorm_x <- abs(qr.R(qr_x)[ncol_x, ncol_x])  ## since dealt with sign earlier
    Y1_scaled <- matrix(Y_tilde[1, ] / fnorm_x, ncol = 1)
    alpha_scaled <- alpha / fnorm_x
    sig_diag_scaled <- sig_diag / (fnorm_x ^ 2)

    ## Fit succotash ---------------------------------------------------------
    if (likelihood == "normal") {
        if (mix_type == 'normal') {
            suc_out <- succotash_given_alpha(Y = Y1_scaled, alpha = alpha_scaled,
                                             sig_diag = sig_diag_scaled,
                                             num_em_runs = num_em_runs, em_z_start_sd = z_start_sd,
                                             lambda_type = lambda_type, lambda0 = lambda0,
                                             tau_seq = tau_seq, em_pi_init = em_pi_init)
        } else if (mix_type == 'uniform') {
            ## right now only runs one em
            ## does not return lfsr
            suc_out <-
                uniform_succ_given_alpha(Y = Y1_scaled, alpha = alpha_scaled,
                                         sig_diag = sig_diag_scaled, num_em_runs = num_em_runs,
                                         em_z_start_sd = z_start_sd, lambda_type = lambda_type)
        }
    } else if (likelihood == "t") {
        suc_out <-
            t_uniform_succ_given_alpha(Y = Y1_scaled, alpha = alpha_scaled, nu = nu,
                                       sig_diag = sig_diag_scaled, num_em_runs = num_em_runs,
                                       em_z_start_sd = z_start_sd, lambda_type = lambda_type)
    }

    suc_out$Y1_scaled <- Y1_scaled  ## ols estimates
    suc_out$alpha_scaled <- alpha_scaled
    suc_out$sig_diag_scaled <- sig_diag_scaled
    return(suc_out)
}

#' Provides posterior summaries in the SUCCOTASH model.
#'
#' \code{succotash_summaries} will return useful posterior summaries
#' used in estimation and testing/FDR control.
#'
#' The posterior distribution is just a mixture of normals. This
#' function will the posterior means as a point estimate. It will also
#' return local false sign and discovery rates. These are useful for
#' controlling for multiple testing.
#'
#' @param Z A matrix of dimension \code{k} by \code{1}. The
#'     (estimated) confounding covariates.
#' @param pi_vals A vector of length \code{M}. The (estimated) mixing
#'     proportions.
#' @inheritParams succotash_given_alpha
#'
#' @return \code{lfdr} (local false discovery rate) A vector of length
#'     \code{p}.  The posterior probability that \eqn{\beta_j = 0}.
#'
#'   \code{lfsr} (local false sign rate) A vector of length
#'   \code{p}. The posterior probability of making a sign error.
#'
#'   \code{betahat} A vector of length \code{p}. The posterior
#'   estimates of \eqn{\beta}.
#'
#' @export
#'
succotash_summaries <- function(Y, Z, pi_vals, alpha, sig_diag, tau_seq) {
    M <- length(tau_seq)  # number of classes in mixture, = K+1 in paper
    p <- nrow(Y)

    if (!is.null(alpha)) {
        mean_mat <- matrix(rep(alpha %*% Z, M), ncol = M, nrow = p)
    } else {
        mean_mat <- matrix(0, ncol = M, nrow = p)
    }
    var_mat <- outer(sig_diag, tau_seq ^ 2, "+")
    top_vals <- t(pi_vals * t(dnorm(matrix(rep(Y, M), ncol = M, nrow = p), mean = mean_mat,
                                    sd = sqrt(var_mat))))
    T <- t(1 / rowSums(top_vals) * top_vals)

    prob.zero <- T[1, ]

    ## get posterior means and variances
    cov_plus <- outer(tau_seq ^ 2, sig_diag, "+")  ## k by p
    if (!is.null(alpha)) {
        post.mean <- (tau_seq ^ 2 / cov_plus) * matrix(Y - alpha %*% Z, ncol = p, nrow = M,
                                                       byrow = TRUE)
    } else {
        post.mean <- (tau_seq ^ 2 / cov_plus) * matrix(Y, ncol = p, nrow = M, byrow = TRUE)
    }
    post.var <- outer(tau_seq ^ 2, sig_diag, "*") / cov_plus

    post.less <- pnorm(q = matrix(0, ncol = p, nrow = M), mean = post.mean, sd = sqrt(post.var))
    post.less[1, ] <- 0

    prob.less <- colSums(T * post.less)

    lfsr <- apply(cbind(prob.less + prob.zero, 1 - prob.less), 1, min)


    beta_hat <- colSums(T * post.mean)

    return(list(lfdr = prob.zero, lfsr = lfsr, beta_hat = beta_hat))
}

#' Transform local false discovery rates to q-values.
#'
#' This function will simply average over lfdr's larger than each
#' component's lfdr.
#'
#' These q-values provide average error rates over subsets of
#' observations. See
#' \href{http://projecteuclid.org/euclid.aos/1074290335}{Storey
#' (2003)} for details.
#'
#' @param lfdr A vector. The local false discovery rates.
#'
#' @return \code{q_vals} A vector of the same length as
#'     \code{lfdr}. Contains the q-values for each observation.
#'
#' @seealso \code{\link{succotash_summaries}},
#'     \code{\link{succotash_given_alpha}}.
lfdr_to_q <- function(lfdr) {
    order_lfdr <- order(lfdr)
    q_vals <- rep(NA, length = length(lfdr))
    q_vals[order_lfdr] <- (cumsum(lfdr[order_lfdr]) / 1:length(lfdr))
    return(q_vals)
}


#' Draw from a mixture of normals.
#'
#' Draw from a mean zero mixture of normals given mixing proportions
#' and mixing standard deviations.
#'
#' Given mixing proportions \code{pi_vals}, and mixing standard
#' deviations (not variances) \code{tau_seq}, this function will draw
#' \code{p} samples from the mean zero mixture of normals.
#'
#' @param pi_vals A vector length \code{M}. The mixing proportions.
#' @param tau_seq A vector of length \code{M}. The mixing standard
#'     deviations.
#' @param p An integer. The number of samples to draw.
#'
#' @export
#'
#' @return beta A vector of length \code{p} that are drawn from the
#'     mean zero mixture of normals.
draw_beta <- function(pi_vals, tau_seq, p) {
    M <- length(pi_vals)
    beta <- rep(NA, length = p)
    which.mix <- sample(1:M, size = p, replace = TRUE, prob = pi_vals)
    for (index in 1:M) {
        current.ind <- which.mix == index
        n_m <- sum(current.ind)
        if (n_m > 0) {
            beta[current.ind] <- rnorm(n = n_m, mean = 0, sd = tau_seq[index])
        }
    }
    return(beta)
}
