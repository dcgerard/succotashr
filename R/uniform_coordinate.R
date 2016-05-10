#' Fixed point for uniform succotash when doing coordinate descent.
#'
#'
#' @inheritParams succotash_unif_fixed
#'
#'
#'
#'
#'
#'
#'
fit_succotash_unif_coord <- function(pi_Z, lambda, alpha, Y, a_seq, b_seq, sig_diag,
                                       print_ziter = FALSE, newt_itermax = 100, tol = 10^-4,
                                       var_scale = TRUE) {

    M <- length(a_seq) + length(b_seq) + 1
    p <- nrow(Y)
    if (var_scale) {
        k         <- length(pi_Z) - M - 1
        scale_val <- pi_Z[length(pi_Z)]
    } else {
        k         <- length(pi_Z) - M
        scale_val <- 1
    }
    pi_old <- pi_Z[1:M]
    if (k != 0) {
        Z_old <- matrix(pi_Z[(M + 1):(M + k)], nrow = k)
    }

    assertthat::assert_that(length(lambda) == M)
    assertthat::are_equal(sum(pi_old), 1)
    assertthat::assert_that(all(a_seq < 0))
    assertthat::assert_that(all(b_seq > 0))


    Z_new <- Z_old
    pi_new <- pi_old

    llike_new <- succotash_llike_unif(pi_Z = pi_Z, lambda = lambda,
                                      alpha = alpha, Y = Y, a_seq = a_seq, b_seq = b_seq,
                                      sig_diag = sig_diag, var_scale = var_scale)

    llike_vec <- llike_new

    err <- tol + 1
    iter_index <- 1
    while (err > tol & iter_index < newt_itermax) {
        llike_old <- llike_new
        pi_Z_old  <- pi_Z
        Z_old     <- Z_new
        pi_old    <- pi_new
        ## update Z with newton step --------------------------------------------------
        optim_out <- optim(par = Z_new, fn = only_Z, gr = only_Z_grad, scale_val = scale_val,
                           pi_vals = pi_new,
                           lambda = lambda, alpha = alpha, Y = Y, a_seq = a_seq, b_seq = b_seq,
                           sig_diag = sig_diag, method = "BFGS",
                           control = list(fnscale = -1, maxit = 50, reltol = 10 ^ -4))
        Z_new <- optim_out$par

        ## optim_out1 <- optim(Z_new, fn = only_Z, scale_val = scale_val, pi_vals = pi_new,
        ##                    lambda = lambda, alpha = alpha, Y = Y, a_seq = a_seq,
        ##                    method = "BFGS",
        ##                    b_seq = b_seq, sig_diag = sig_diag, control = list(fnscale = -1))
        if (var_scale) {
           pi_Z <- c(pi_new, optim_out$par, scale_val)
        } else {
           pi_Z <- c(pi_new, optim_out$par)
        }
        ## Z_new <- optim_out$par


        ## cat("llike after Z:", optim_out$value, "\n")
        llike_new <- optim_out$value
        llike_vec <- c(llike_vec, llike_new)

        assertthat::are_equal(pi_Z_old[1:M], pi_Z[1:M])
        assertthat::are_equal(sum(pi_new), 1)

        ## update pi with ashr -------------------------------------------------------
        betahat <- Y - alpha %*% Z_new
        if (var_scale) {
            sebetahat <- sqrt(sig_diag * scale_val)
        } else {
            sebetahat <- sqrt(sig_diag)
        }


        g = ashr::unimix(pi_new, c(a_seq, rep(0, length(b_seq) + 1)),
                         c(rep(0, length(a_seq) + 1), b_seq))

        if (requireNamespace(package = "Rmosek", quietly = TRUE)) {
            optmethod = "mixIP"
        } else {
            optmethod = "mixEM"
        }
        control.default = list(K = 1, method = 3, square = TRUE,
                               step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1,
                               tol = 1e-07, maxiter = 500, trace = FALSE)
        ash_out <- ashr:::estimate_mixprop(betahat = c(betahat), sebetahat = sebetahat,
                                           g = g, prior = lambda, null.comp = length(a_seq + 1),
                                           optmethod = optmethod, control = control.default)

        ## ashr subtracts off the largest value of the loglikelihood
        ## to increase numerical stability during the optimization
        ## step, so the log-likelihood will be artificially smaller
        ## than what I calculate.

        pi_new <- ash_out$g$pi

        assertthat::are_equal(sum(pi_new), 1)
        if (var_scale) {
            pi_Z <- c(pi_new, Z_new, scale_val)
        } else {
            pi_Z <- c(pi_new, Z_new)
        }

        ## llike_new <- succotash_llike_unif(pi_Z = pi_Z, lambda = lambda,
        ##                                   alpha = alpha, Y = Y, a_seq = a_seq, b_seq = b_seq,
        ##                                   sig_diag = sig_diag, var_scale = var_scale)

        ## cat("llike after ash:", llike_new, "\n")
        ## cat(" llike ash says:", ash_out$loglik, "\n")


        if (var_scale) {
            oout <- optim(par = scale_val, fn = only_scale,
                          pi_Z_minus_scale = pi_Z[-length(pi_Z)],
                          lambda = lambda, alpha = alpha, Y = Y, a_seq = a_seq,
                          b_seq = b_seq, sig_diag = sig_diag, method = "Brent", lower = 0,
                          upper = 10, control = list(fnscale = -1))
            pi_Z <- c(pi_Z[1:(length(pi_Z) - 1)], oout$par)
            scale_val <- oout$par
        }


        err <- sum(abs(pi_Z - pi_Z_old))
        ## err <- abs(llike_new / llike_old - 1)
        ## cat("err:", err, "\n")
    }

    return(list(pi_Z = pi_Z, llike_vec = llike_vec))
}

#' Wrapper for \code{\link{succotash_llike_unif}} but useful for optim where only update scale_val.
only_scale <- function(scale_val, pi_Z_minus_scale, lambda, alpha, Y, a_seq, b_seq, sig_diag) {
    pi_Z <- c(pi_Z_minus_scale, scale_val)
    llike <- succotash_llike_unif(pi_Z = pi_Z, lambda = lambda,
                                  alpha = alpha, Y = Y, a_seq = a_seq,
                                  b_seq = b_seq, sig_diag = sig_diag,
                                  var_scale = TRUE)
    return(llike)
}

#' Wrapper for \code{\link{succotash_llike_unif}} but useful for optim where only update Z.
#'
#' @param Z A k by 1 matrix of numerics. The confounders.
#' @param pi_vals A M vector of numerics. The mixing probs.
#' @param scale_val A positive numeric. The variance inflation parameter.
#' @inheritParams succotash_llike_unif
#'
only_Z <- function(Z, pi_vals, scale_val, lambda, alpha, Y, a_seq, b_seq, sig_diag) {
    pi_Z <- c(pi_vals, Z, scale_val)
    llike <- succotash_llike_unif(pi_Z = pi_Z, lambda = lambda,
                                  alpha = alpha, Y = Y, a_seq = a_seq,
                                  b_seq = b_seq, sig_diag = sig_diag,
                                  var_scale = TRUE)
    return(llike)
}

#' Wrapper for \code{\link{unif_grad_simp}}, but useful for optim where only update Z.
#'
#' @inheritParams only_Z
#'
only_Z_grad <- function(Z, pi_vals, scale_val, lambda, alpha, Y, a_seq, b_seq, sig_diag) {
  pi_Z <- c(pi_vals, Z, scale_val)
  grad_final <- unif_grad_simp(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,  sig_diag = sig_diag,
                               a_seq = a_seq, b_seq = b_seq, var_scale = TRUE)
  return(grad_final)
}

#' My first attempt at calculating the gradient of the likelihood
#' function using uniform mixtures wrt Z.
#'
#' This doesn't seem to work too well.
#'
#' @inheritParams succotash_llike_unif
#'
unif_grad_llike <- function(pi_Z, lambda, alpha, Y, a_seq, b_seq, sig_diag,
                                 var_scale = TRUE) {
    M <- length(a_seq) + length(b_seq) + 1
    ## p <- nrow(Y)
    if (var_scale) {
        k <- length(pi_Z) - M - 1
        scale_val <- pi_Z[length(pi_Z)]
    } else {
        k <- length(pi_Z) - M
        scale_val <- 1
    }
    pi_current <- pi_Z[1:M]
    if (k != 0) {
        Z_current <- matrix(pi_Z[(M + 1):(M + k)], nrow = k)
    }

    assertthat::are_equal(length(lambda), M)
    assertthat::assert_that(all(lambda >= 1))
    assertthat::assert_that(scale_val > 0)
    assertthat::are_equal(sum(pi_current), 1, tol = 10 ^ -4)

    sig_diag <- sig_diag * scale_val

    az <- alpha %*% Z_current

    left_means <- diag(1 / sqrt(sig_diag)) %*% outer(c( (Y - az)), a_seq, "-")
    left_means_zero <- diag(1 / sqrt(sig_diag)) %*% outer(c( (Y - az)), rep(0, length(a_seq)), "-")
    right_means <- diag(1 / sqrt(sig_diag)) %*% outer(c( (Y - az)), b_seq, "-")
    right_means_zero <- diag(1 / sqrt(sig_diag)) %*% outer(c( (Y - az)), rep(0, length(b_seq)), "-")
    zero_means <- dnorm(Y, mean = az, sd = sqrt(sig_diag))

    left_ispos <- left_means > 0
    right_ispos <- right_means > 0

    pnorm_left_diff <- matrix(NA, ncol = ncol(left_means), nrow = nrow(left_means))
    pnorm_right_diff <- matrix(NA, ncol = ncol(right_means), nrow = nrow(right_means))

    pnorm_left_diff[!left_ispos] <-
        pnorm(left_means[!left_ispos]) - pnorm(left_means_zero[!left_ispos])
    pnorm_right_diff[!right_ispos] <-
        pnorm(right_means_zero[!right_ispos]) - pnorm(right_means[!right_ispos])

    pnorm_left_diff[left_ispos] <-
        pnorm(-1 * left_means_zero[left_ispos]) - pnorm(-1 * left_means[left_ispos])
    pnorm_right_diff[right_ispos] <-
        pnorm(-1 * right_means[right_ispos]) - pnorm(-1 * right_means_zero[right_ispos])

    pnorm_left_diff  <- t(t(pnorm_left_diff) / a_seq)
    pnorm_right_diff <- t(t(pnorm_right_diff) / b_seq)

    fkj <- abs(cbind(pnorm_left_diff,
                     zero_means,
                     pnorm_right_diff))

    denom_grad <- rowSums(t(t(fkj) * pi_current))
    ## fkj %*% diag(pi_current)

    mid_part <- zero_means * (Y - az) / sig_diag


    dnorm_left_diff <- dnorm(left_means_zero) - dnorm(left_means)
    dnorm_right_diff <- dnorm(right_means) - dnorm(right_means_zero)

    ##diag(1 / sig_diag) %*% dnorm_left_diff %*% diag(1 / a_seq)
    left_part <- t(t(1 / sqrt(sig_diag) * dnorm_left_diff) / a_seq)
    right_part <- t(t(1 / sqrt(sig_diag) * dnorm_right_diff) / b_seq)

    top_mult <- rowSums(t(t(cbind(left_part, mid_part, right_part)) * pi_current))

    ## t(t(cbind(left_part, mid_part, right_part)) * pi_current) -
    ##      cbind(left_part, mid_part, right_part) %*% diag(pi_current)

    mult_vals <- top_mult / denom_grad

    gradient_val <- colSums(mult_vals * alpha)
    ## diag(mult_vals) %*% alpha

    ## cat(gradient_val, "\n\n")

    if (var_scale) {
        augmented_grad <- c(rep(0, length = length(pi_current)), gradient_val, 0)
    } else {
        augmented_grad <- c(rep(0, length = length(pi_current)), gradient_val)
    }

    return(augmented_grad)
}


#' Me being as careful as possible when calculating the succotash log-likelihood.
#'
#' @param Y a p by 1 matrix of numerics.
#' @param Z a k by 1 matrix of numerics.
#' @param alpha A p by k matrix of numerics.
#' @param sig_diag A p-vector of numerics.
#' @param scale_val A positivie numeric.
#' @param left_seq The left endpoints of the uniforms.
#' @param right_seq The right endpoints of the uniforms
#'
#'
#'
llike_unif_simp <- function(Y, Z, pi_vals, alpha, sig_diag, left_seq, right_seq,
                            scale_val = 1) {
    sig_diag <- sig_diag * scale_val

    p <- nrow(Y)
    k <- nrow(Z)
    M <- length(left_seq)
    assertthat::are_equal(M, length(right_seq))

    null_spot <- which(abs(left_seq) < 10^-12 & abs(right_seq) < 10^-12)


    resid_vec <- Y - alpha %*% Z

    mean_mat_left  <- outer(c(resid_vec), left_seq, "-")
    mean_mat_right <- outer(c(resid_vec), right_seq, "-")

    fkj_mat <- matrix(NA, nrow = p, ncol = M)

    var_mat <- matrix(rep(sqrt(sig_diag), M - length(null_spot)), nrow = p, byrow = FALSE)
    pnorm_diff <- pnorm(mean_mat_left[, -null_spot], mean = 0, sd = var_mat) -
        pnorm(mean_mat_right[, -null_spot], mean = 0, sd = var_mat)

    pnorm_ratio <- pnorm_diff %*% diag(1 / (right_seq[-null_spot] - left_seq[-null_spot]))


    fkj_mat[, -null_spot] <- pnorm_ratio

    assertthat::are_equal(mean_mat_left[, null_spot], mean_mat_right[, null_spot])

    fkj_mat[, null_spot] <- dnorm(mean_mat_left[, null_spot], mean = 0, sd = sqrt(sig_diag))

    mat_lik <- fkj_mat %*% diag(pi_vals)

    llike <- sum(log(rowSums(mat_lik)))
    return(llike)
}

#' This is me being very careful about calculating the gradient of the
#' likelihood function wrt Z.
#'
#' @inheritParams succotash_llike_unif
#' @inheritParams llike_unif_simp
#'
#' @return The gradient for Z.
#'
unif_grad_simp <- function(pi_Z, lambda, alpha, Y,  sig_diag, left_seq = NULL, right_seq = NULL,
                           a_seq = NULL, b_seq = NULL, var_scale = TRUE) {

    assertthat::assert_that(!((is.null(left_seq) | is.null(right_seq)) &
                              (is.null(a_seq) | is.null(b_seq))))

    if (is.null(left_seq) | is.null(right_seq)) {
        left_seq <- c(a_seq, rep(0, length = length(b_seq) + 1))
        right_seq <- c(rep(0, length = length(a_seq) + 1), b_seq)
    }

    M <- length(left_seq)
    p <- nrow(Y)
    if (var_scale) {
        k <- length(pi_Z) - M - 1
        scale_val <- pi_Z[length(pi_Z)]
    } else {
        k <- length(pi_Z) - M
        scale_val <- 1
    }
    pi_current <- pi_Z[1:M]
    if (k != 0) {
        Z_current <- matrix(pi_Z[(M + 1):(M + k)], nrow = k)
    }

    assertthat::are_equal(length(left_seq), length(right_seq))
    assertthat::are_equal(length(lambda), M)
    assertthat::assert_that(all(lambda >= 1))
    assertthat::assert_that(scale_val > 0)
    assertthat::are_equal(sum(pi_current), 1, tol = 10 ^ -4)

    sig_diag <- sig_diag * scale_val

    az <- alpha %*% Z_current

    resid_vec <- Y - az

    null_spot <- which(abs(left_seq) < 10^-12 & abs(right_seq) < 10^-12)

    left_means  <- matrix(rep(left_seq, p), nrow = p,byrow = TRUE)
    right_means <- matrix(rep(right_seq, p), nrow = p,byrow = TRUE)

    quants <- matrix(rep(resid_vec, M), nrow = p, byrow = FALSE)
    sd_mat <- matrix(rep(sig_diag, M), nrow = p, byrow = FALSE)

    dnorm_mat_left  <- dnorm(x = quants, mean = left_means, sd = sd_mat)
    dnorm_mat_right <- dnorm(x = quants, mean = right_means, sd = sd_mat)

    dnorm_diff_mat <- dnorm_mat_right[, -null_spot] - dnorm_mat_left[, -null_spot]

    dnorm_diff_mat_scaled <- matrix(NA, nrow = p, ncol = M)
    dnorm_diff_mat_scaled[, -null_spot] <-
        dnorm_diff_mat %*% diag(1 / (right_seq[-null_spot] - left_seq[-null_spot]))

    assertthat::are_equal(dnorm_mat_left[, null_spot], dnorm_mat_right[, null_spot])

    dnorm_diff_mat_scaled[, null_spot] <- dnorm_mat_left[, null_spot] * (resid_vec / sig_diag)

    dnorm_diff_mat_pluspi <- dnorm_diff_mat_scaled %*% diag(pi_current)

    top_vec <- rowSums(dnorm_diff_mat_pluspi)


    ## Now do likelihood for bottom part ----------------------------------------------------
    mean_mat_left  <- outer(c(resid_vec), left_seq, "-")
    mean_mat_right <- outer(c(resid_vec), right_seq, "-")

    fkj_mat <- matrix(NA, nrow = p, ncol = M)

    var_mat <- matrix(rep(sqrt(sig_diag), M - length(null_spot)), nrow = p, byrow = FALSE)
    pnorm_diff <- pnorm(mean_mat_left[, -null_spot], mean = 0, sd = var_mat) -
        pnorm(mean_mat_right[, -null_spot], mean = 0, sd = var_mat)

    pnorm_ratio <- pnorm_diff %*% diag(1 / (right_seq[-null_spot] - left_seq[-null_spot]))


    fkj_mat[, -null_spot] <- pnorm_ratio

    assertthat::are_equal(mean_mat_left[, null_spot], mean_mat_right[, null_spot])

    fkj_mat[, null_spot] <- dnorm(mean_mat_left[, null_spot], mean = 0, sd = sqrt(sig_diag))

    mat_lik <- fkj_mat %*% diag(pi_current)

    bottom_vec <- rowSums(mat_lik)

    ## Now combine -------------------------------------------------------------------------

    weights_vec <- top_vec / bottom_vec

    weighted_alpha <- diag(weights_vec) %*% alpha

    grad_final <- colSums(weighted_alpha)

    if (var_scale) {
        augmented_grad <- c(rep(0, length = M), grad_final, 0)
    } else {
        augmented_grad <- c(rep(0, length = M), grad_final)
    }

    return(grad_final)
}
