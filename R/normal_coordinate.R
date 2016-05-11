#' Coordinate ascent for normal mixtures with normal likelihood.
#'
#' @inheritParams succotash_fixed
#' @inheritParams succotash_given_alpha
#'
normal_coord <- function(pi_Z, lambda, alpha, Y, tau_seq, sig_diag,
                         plot_new_ests = FALSE, var_scale = TRUE,
                         itermax = 100, tol = 10 ^ -4) {
    M <- length(tau_seq)
    p <- nrow(Y)
    k <- length(pi_Z) - M - var_scale ## var_scale is 0 if FALSE, 1 if TRUE
    pi_new <- pi_Z[1:M]
    if (k != 0) {
        Z_new <- matrix(pi_Z[(M + 1):(M + k)], nrow = k)
    }

    if (var_scale) {
        scale_val <- pi_Z[M + k + 1] ## the amount to scale the variance by
    } else {
        scale_val <- 1
    }

    assertthat::are_equal(sum(pi_new), 1)
    assertthat::are_equal(length(lambda), M)
    assertthat::assert_that(all(tau_seq >= 0))

    llike_new <- succotash_llike(pi_Z = pi_Z, lambda = lambda,
                                 alpha = alpha, Y = Y,
                                 tau_seq = tau_seq,
                                 sig_diag = sig_diag,
                                 var_scale = var_scale)
    llike_vec <- llike_new

    err <- tol + 1
    iter_index <- 1
    while (err > tol & iter_index < itermax) {
        llike_old <- llike_new
        pi_old    <- pi_new
        Z_old     <- Z_new
        pi_Z_old  <- pi_Z

        ## Update Z ---------------------------------------------------------
        optim_out <- stats::optim(par = Z_new, fn = normal_only_z,
                                  gr = normal_llike_grad,
                                  scale_val = scale_val,
                                  pi_vals = pi_new,
                                  alpha = alpha, Y = Y,
                                  tau_seq = tau_seq,
                                  sig_diag = sig_diag,
                                  method = "BFGS",
                                  control = list(fnscale = -1, maxit = 50, reltol = 10 ^ -4))
        Z_new <- optim_out$par

        if (var_scale) {
            pi_Z <- c(pi_new, optim_out$par, scale_val)
        } else {
            pi_Z <- c(pi_new, optim_out$par)
        }

        cat("llike after Z:", optim_out$value, "\n")
        llike_new <- optim_out$value
        llike_vec <- c(llike_vec, llike_new)

        ## Update pi with ashr --------------------------------------------
        betahat <- Y - alpha %*% Z
        sebetahat <- sqrt(sig_diag * scale_val)
        g <- ashr::normalmix(pi = pi_new, mean = rep(0, length = M), sd = tau_seq)

        if (requireNamespace(package = "Rmosek", quietly = TRUE)) {
            optmethod <- "mixIP"
        } else {
            optmethod <- "mixEM"
        }
        control.default <- list(K = 1, method = 3, square = TRUE,
                               step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1,
                               tol = 1e-07, maxiter = 500, trace = FALSE)
        ash_out <- ashr:::estimate_mixprop(betahat = c(betahat), sebetahat = sebetahat,
                                           g = g, prior = lambda, null.comp = 1,
                                           optmethod = optmethod, control = control.default)

        pi_new <- ash_out$g$pi

        assertthat::are_equal(sum(pi_new), 1)
        if (var_scale) {
            pi_Z <- c(pi_new, Z_new, scale_val)
        } else {
            pi_Z <- c(pi_new, Z_new)
        }

        llike_new <- succotash_llike(pi_Z = pi_Z, lambda = lambda,
                                     alpha = alpha, Y = Y,
                                     tau_seq = tau_seq,
                                     sig_diag = sig_diag,
                                     var_scale = var_scale)
        cat("llike after ash:", llike_new, "\n")


        ## Update scale_val with Brent's method -----------------------------------------
        if (var_scale) {
            oout <- stats::optim(par = scale_val, fn = normal_only_z,
                                 Z = Z_new, pi_vals = pi_new,
                                 lambda = lambda, alpha = alpha,
                                 Y = Y, tau_seq = tau_seq,
                                 sig_diag = sig_diag,
                                 method = "Brent", lower = 0,
                                 upper = 10,
                                 control = list(fnscale = -1))
            pi_Z <- c(pi_Z[1:(length(pi_Z) - 1)], oout$par)
            scale_val <- oout$par
        }

        err <- sum(abs(pi_Z - pi_Z_old))
    }

    return(list(pi_Z = pi_Z, llike_vec = llike_vec))
}

#' Wrapper for \code{\link{succotash_llike}}, but useful when calling
#' optim.
#'
#' @inheritParams succotash_llike
#' @param Z A k by 1 matrix of numerics.
#' @param pi_vals An M vector of numerics that sum to one.
#' @param scale_val A positive numeric. The variance inflation
#'     parameter.
#'
normal_only_z <- function(Z, pi_vals, scale_val, lambda, alpha, Y,
                          tau_seq, sig_diag) {
    pi_Z <- c(pi_vals, Z, scale_val)
    llike_new <- succotash_llike(pi_Z = pi_Z, lambda = lambda,
                                 alpha = alpha, Y = Y,
                                 tau_seq = tau_seq,
                                 sig_diag = sig_diag,
                                 var_scale = TRUE)
    return(llike_new)
}

#' Gradient of the log-likelihood wrt Z.
#'
#'
normal_llike_grad <- function(Z, Y, alpha, sig_diag, tau_seq, scale_val, pi_vals) {

    p <- nrow(Y)
    k <- nrow(Z)
    M <- length(pi_vals)

    resid_vec <- Y - alpha %*% Z

    quant_matrix <- matrix(rep(c(resid_vec), times = M), nrow = p, byrow = FALSE)
    var_mat <- outer(sig_diag, tau_seq ^ 2, FUN = "+")

    dnorm_mat <- stats::dnorm(x = quant_matrix, mean = 0, sd = sqrt(var_mat))
    like_mat <- dnorm_mat %*% diag(pi_vals)
    resid_mat <- quant_matrix / var_mat

    dlike_mat <- like_mat * resid_mat

    top_vals <- rowSums(dlike_mat)
    bottom_vals <- rowSums(like_mat)

    weight_vec <- top_vals / bottom_vals

    alpha_weighted <- diag(weight_vec) %*% alpha

    grad_final <- colSums(alpha)

    return(grad_final)
}

#' Me being super careful in calculating the normal mixtures likelihood.
normal_llike_simp <- function(Z, Y, alpha, sig_diag, tau_seq, scale_val, pi_vals) {

    p <- nrow(Y)
    k <- nrow(Z)
    M <- length(pi_vals)

    resid_vec <- Y - alpha %*% Z

    quant_matrix <- matrix(rep(c(resid_vec), times = M), nrow = p, byrow = FALSE)
    var_mat <- outer(sig_diag, tau_seq ^ 2, FUN = "+")

    dnorm_mat <- stats::dnorm(x = quant_matrix, mean = 0, sd = sqrt(var_mat))
    like_mat <- dnorm_mat %*% diag(pi_vals)

    llike <- sum(log(rowSums(like_mat)))

    return(llike)
}
