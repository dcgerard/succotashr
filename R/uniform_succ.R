
#'A fixed point iteration in the mixture of uniforms EM.
#'
#'This is a fixed-point iteration for the SUCCOTASH EM algorithm. This
#' updates the estimate of the prior and the estimate of the hidden
#' covariates.
#'
#'@param pi_Z A vector. The first \code{M} values are the current
#'     values of \eqn{\pi}. The last \code{k} values are the current
#'     values of \eqn{Z}.
#'@param lambda A vector. This is a length \code{M} vector with the
#'     regularization parameters for the mixing proportions.
#'@param alpha A matrix. This is of dimension \code{p} by \code{k} and
#'     are the coefficients to the confounding variables.
#'@param a_seq A vector of negative numerics containing the left
#'     endpoints of the mixing uniforms.
#'@param b_seq A vector of positiv numerics containing the right
#'     endpoints of the mixing uniforms.
#'@param print_ziter A logical. Should we we print each iteration of
#'     the Z optimization?
#'@param newt_itermax A positive integer. The maximum number of Newton
#'     steps to perform in updating \eqn{Z}.
#'@param tol A positive numeric. The stopping criterion for Newton's
#'     method in updating Z.
#' @inheritParams succotash_given_alpha
#'
#'@return \code{pi_new} A vector of length \code{M}. The update for
#'     the mixing components.
#'
#'  \code{Z_new} A vector of length \code{k}. The update for the
#'  confounder covariates.
#'
#'@seealso \code{\link{uniform_succ_given_alpha}}
#'     \code{\link{succotash_llike_unif}}
#'
#' @export
succotash_unif_fixed <- function(pi_Z, lambda, alpha, Y, a_seq, b_seq, sig_diag,
                                 print_ziter = FALSE, newt_itermax = 4, tol = 10 ^ -4,
                                 var_scale = TRUE) {
    M <- length(a_seq) + length(b_seq) + 1
    ## p <- nrow(Y)
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
    assertthat::are_equal(sum(pi_old), 1, tol = 10 ^ -4)
    assertthat::assert_that(all(a_seq < 0))
    assertthat::assert_that(all(b_seq > 0))


    update_pi_out <- unif_update_pi(pi_old = pi_old, Y = Y, sig_diag = sig_diag, alpha = alpha,
                                    a_seq = a_seq, b_seq = b_seq, Z_old = Z_old, lambda = lambda)
    pi_new    <- update_pi_out$pi_new
    llike_old <- update_pi_out$llike_old
    Tkj       <- update_pi_out$Tkj

    scale_val_new <- scale_val
    if (var_scale) {
        scale_val_tol   <- 10 ^ -3
        scale_val_err   <- scale_val_tol + 1
        scale_val_index <- 1
        while (scale_val_index <= 10 & scale_val_err > scale_val_tol) {
            scale_val_old <- scale_val_new
            update_z_out  <- unif_update_z(Y = Y, sig_diag = sig_diag,
                                           alpha = alpha,
                                           lambda = lambda,
                                           a_seq = a_seq,
                                           b_seq = b_seq,
                                           pi_new = pi_new, Tkj = Tkj,
                                           Z_old = Z_old,
                                           llike_old = llike_old,
                                           tol = tol,
                                           newt_itermax = newt_itermax,
                                           print_ziter = print_ziter,
                                           scale_val = scale_val_new)
            Z_new     <- update_z_out$Z_new
            llike_new <- update_z_out$llike_new
            ## cat(llike_new, "\n")
            oout <- stats::optim(par = scale_val, fn = fun_unif_scale,
                                 Y = Y, az = alpha %*% Z_new,
                                 sig_diag = sig_diag, Tkj = Tkj,
                                 a_seq = a_seq, b_seq = b_seq,
                                 method = "Brent", lower = 0,
                                 upper = 10,
                                 control = list(fnscale = -1, maxit = 10))
            scale_val_new <- oout$par

            scale_val_err   <- abs(scale_val_new / scale_val_old - 1)
            scale_val_index <- scale_val_index + 1
        }
    } else {
        update_z_out <- unif_update_z(Y = Y, sig_diag = sig_diag,
                                      alpha = alpha, lambda = lambda,
                                      a_seq = a_seq, b_seq = b_seq,
                                      pi_new = pi_new, Tkj = Tkj,
                                      Z_old = Z_old,
                                      llike_old = llike_old,
                                      tol = tol,
                                      newt_itermax = newt_itermax,
                                      print_ziter = print_ziter,
                                      scale_val = scale_val_new)
        Z_new     <- update_z_out$Z_new
        llike_new <- update_z_out$llike_new
        ## cat(llike_new, "\n")
    }


    if (var_scale) {
        pi_Z <- c(pi_new, Z_new, scale_val_new)
    } else {
        pi_Z <- c(pi_new, Z_new)
    }

    return(pi_Z = pi_Z)
}

#' E step in EM algorithm.
#'
#'
#' @inheritParams succotash_unif_fixed
#' @param Z_old A vector of numerics. The old confounders.
#' @param pi_old A vector of numerics. The old pi values.
unif_update_pi <- function(pi_old, Y, sig_diag, alpha, a_seq, b_seq,
                           Z_old, lambda) {
    az <- alpha %*% Z_old
    left_seq  <- c(a_seq, rep(0, length = length(b_seq)))
    right_seq <- c(rep(0, length = length(a_seq)), b_seq)

    p <- nrow(Y)
    M <- length(pi_old)

    left_means <- diag(1 / sqrt(sig_diag)) %*%
        outer(c(Y - az), left_seq, "-")
    right_means <- diag(1 / sqrt(sig_diag)) %*%
        outer(c(Y - az), right_seq, "-")

    ## pnorm(left_means) - pnorm(right_means) ## this is equal to pdiff

    left_bigger <- abs(left_means) > abs(right_means)
    which_switch <- (left_bigger * sign(left_means) + (!left_bigger) * sign(right_means)) == 1

    pdiff <- matrix(NA, nrow = nrow(left_means), ncol = ncol(left_means))
    pdiff[!which_switch] <- stats::pnorm(left_means[!which_switch]) -
        stats::pnorm(right_means[!which_switch])
    pdiff[which_switch] <- stats::pnorm(-1 * right_means[which_switch]) -
        stats::pnorm(-1 * left_means[which_switch])

    ## calculate new pi values
    fkj    <- pdiff %*% diag(1 / abs(right_seq - left_seq))
    f0j    <- stats::dnorm(Y, az, sqrt(sig_diag))
    fkj    <- cbind(fkj[, 1:length(a_seq)], f0j, fkj[, (length(a_seq) + 1):ncol(fkj)])
    fkj_pi <- fkj %*% diag(pi_old)
    Tkj    <- diag(1 / rowSums(fkj_pi)) %*% fkj_pi

    llike_old <- sum(log(rowSums(fkj_pi)))

    pi_new <- (colSums(Tkj) + lambda - 1) / (p - M + sum(lambda))

    return(list(pi_new = pi_new, llike_old = llike_old, Tkj = Tkj))
}


#' Runs a few newton steps to update Z given scale_val.
#'
#' @inheritParams succotash_unif_fixed
#' @param pi_new A vector of numerics. The current estimates of of the
#'     mixing proportions.
#' @param Tkj A matrix of numerics. The T matrix.
#' @param Z_old A vector of numerics. The old confounders.
#' @param llike_old The old log-likelihood
#' @param scale_val A positive numeric. The variance inflation parameter.
unif_update_z <- function(Y, sig_diag, alpha, lambda, a_seq, b_seq, pi_new, Tkj, Z_old, llike_old,
                          tol = 10 ^ -3, newt_itermax = 10, print_ziter = FALSE,
                          scale_val = 1) {
    z_diff <- tol + 1
    newt_iter <- 1
    Z_new <- Z_old
    mult_val <- 1
    llike_new <- llike_old
    while (z_diff > tol & newt_iter < newt_itermax) {
        did_I_move <- TRUE

        llike_old <- llike_new
        Z_old <- Z_new

        az <- alpha %*% Z_old

        left_means <- diag(1 / sqrt(sig_diag)) %*%
            outer(c(Y - az), a_seq, "-")
        right_means <- diag(1 / sqrt(sig_diag)) %*%
            outer(c(Y - az), b_seq, "-")
        zero_means_left <- diag(1 / sqrt(sig_diag)) %*%
            outer(c(Y - az), rep(0, length = length(a_seq)), "-")
        zero_means_right <- diag(1 / sqrt(sig_diag)) %*%
            outer(c(Y - az), rep(0, length = length(b_seq)), "-")

        ## negative means are more accurate, so switch
        left_ispos <- left_means > 0
        right_ispos <- right_means > 0

        pnorm_diff_left <- matrix(NA, nrow = nrow(left_means), ncol = ncol(left_means))
        pnorm_diff_right <- matrix(NA, nrow = nrow(right_means), ncol = ncol(right_means))

        pnorm_diff_left[left_ispos] <-
            stats::pnorm(-1 * zero_means_left[left_ispos]) - stats::pnorm(-1 * left_means[left_ispos])
        pnorm_diff_right[right_ispos] <-
            stats::pnorm(-1 * right_means[right_ispos]) - stats::pnorm(-1 * zero_means_right[right_ispos])

        pnorm_diff_left[!left_ispos] <-
            stats::pnorm(left_means[!left_ispos]) - stats::pnorm(zero_means_left[!left_ispos])
        pnorm_diff_right[!right_ispos] <-
            stats::pnorm(zero_means_right[!right_ispos]) - stats::pnorm(right_means[!right_ispos])

        ## calculate gradient
        dnorm_diff_left <- diag(1 / sqrt(sig_diag)) %*%
            (stats::dnorm(zero_means_left) - stats::dnorm(left_means))
        dnorm_diff_right <- diag(1 / sqrt(sig_diag)) %*%
            (stats::dnorm(right_means) - stats::dnorm(zero_means_right))
        dpratios_left <- dnorm_diff_left / pnorm_diff_left
        dpratios_right <- dnorm_diff_right / pnorm_diff_right

        ## ad-hoc adjustment for numerical instability of truncated normal. May need to remove -------------------
        dpratios_left[is.na(dpratios_left)] <- 0
        dpratios_left[dpratios_left == -Inf] <- min(dpratios_left[dpratios_left != -Inf], na.rm = TRUE)
        dpratios_left[dpratios_left == Inf] <- max(dpratios_left[dpratios_left != Inf], na.rm = TRUE)

        dpratios_right[is.na(dpratios_right)] <- 0
        dpratios_right[dpratios_right == -Inf] <- min(dpratios_right[dpratios_right != -Inf], na.rm = TRUE)
        dpratios_right[dpratios_right == Inf] <- min(dpratios_right[dpratios_right != Inf], na.rm = TRUE)
        ## -------------------------------------------------------------------------------------------------------

        zero_part <- (Y - az) / sig_diag
        alpha_weights <- rowSums(Tkj * cbind(dpratios_left, zero_part, dpratios_right))


        gradient_val <- colSums(diag(alpha_weights) %*% alpha)

        ##calculate Hessian
        top_left <- left_means * stats::dnorm(left_means) -
            zero_means_left *  stats::dnorm(zero_means_left)
        top_right <- zero_means_right * stats::dnorm(zero_means_right) -
            right_means * stats::dnorm(right_means)



        sum1 <- top_left / pnorm_diff_left - dpratios_left ^ 2
        sum2 <- top_right / pnorm_diff_right - dpratios_right ^ 2
        sum0 <- -1 / sig_diag

        ## ad-hoc adjustment for numerical instability of truncated normal. May need to remove -------------------
        sum1[is.na(sum1)] <- 0
        sum2[is.na(sum2)] <- 0

        sum1[sum1 == Inf] <- max(sum1[sum1 != Inf])
        sum2[sum2 == Inf] <- max(sum2[sum2 != Inf])

        sum1[sum1 == -Inf] <- min(sum1[sum1 != -Inf])
        sum2[sum2 == -Inf] <- min(sum2[sum2 != -Inf])
        ## -------------------------------------------------------------------------------------------------------

        diag_weights <- rowSums(Tkj * cbind(sum1, sum0, sum2))

        hessian_val <- t(alpha) %*% diag(diag_weights) %*% alpha

        Z_new <- Z_old -  mult_val * solve(hessian_val) %*% gradient_val

        if (any(is.na(Z_new))) {
            stop("Z_new contains NA's")
        }

        z_diff <- sum(abs(Z_new - Z_old))

        llike_new <- succotash_llike_unif(pi_Z = c(pi_new, Z_new, scale_val), lambda = lambda,
                                          alpha = alpha, Y = Y, a_seq = a_seq,
                                          b_seq = b_seq, sig_diag = sig_diag, var_scale = TRUE)

        if (llike_new < llike_old) {
            ## halve the step size and reset to old Z value
            mult_val <- mult_val / 2
            Z_new <- Z_old
            llike_new <- llike_old
            did_I_move <- FALSE
        }

        if (did_I_move == FALSE & newt_iter == 1) {
            break
        }

        if (print_ziter) {
            cat("Iteration =", newt_iter, "\n")
            cat("Did I move?", did_I_move, "\n")
            cat("Z Diff = ", z_diff, "\n")
            cat(llike_new)
        }
        newt_iter <- newt_iter + 1
    }
    return(list(Z_new = Z_new, llike_new = llike_new))
}


#' Calculate the difference in cdfs.
#'
#' @param resid_vec A vector of numerics such that should be distributed mean zero.
#' @param sig_diag A p-vector of numerics. The estimated variances.
#' @param a_seq A vector of numerics. The non-zero left ends of the
#'     uniforms.
#' @param b_seq A vector of numerics. The non-zero right ends of the
#'     uniforms.
#' @param like_type A character. Available options are "normal" and "t".
#'
#'
#'
get_pdiffs <- function(resid_vec, sig_diag, a_seq, b_seq, like_type = "normal") {
    left_seq  <- c(a_seq, rep(0, length = length(b_seq)))
    right_seq <- c(rep(0, length = length(a_seq)), b_seq)

    left_means <- diag(1 / sqrt(sig_diag)) %*%
        outer(resid_vec, left_seq, "-")
    right_means <- diag(1 / sqrt(sig_diag)) %*%
        outer(resid_vec, right_seq, "-")

    ## pnorm(left_means) - pnorm(right_means) ## this is equal to pdiff

    left_bigger <- abs(left_means) > abs(right_means)
    which_switch <- (left_bigger * sign(left_means) + (!left_bigger) * sign(right_means)) == 1

    pdiff <- matrix(NA, nrow = nrow(left_means), ncol = ncol(left_means))
    pdiff[!which_switch] <- stats::pnorm(left_means[!which_switch]) -
        stats::pnorm(right_means[!which_switch])
    pdiff[which_switch] <- stats::pnorm(-1 * right_means[which_switch]) -
        stats::pnorm(-1 * left_means[which_switch])
    return(pdiff)
}



#' Criterion function for the scale parameter as an intermediate step
#' to the EM algorithm when using uniform mixtures.
#'
#'
#' @param scale_val A positive numeric. The variance inflation
#'     parameter.
#' @param Tkj An m by p matrix of numerics. The "T-values".
#' @param Y A p by 1 matrix of numerics. The data.
#' @param az A p by 1 matrix. This is \eqn{alpha Z}.
#' @param sig_diag A p-vector of numerics. The estimated variances.
#' @param a_seq A vector of numerics. The non-zero left ends of the
#'     uniforms.
#' @param b_seq A vector of numerics. The non-zero right ends of the
#'     uniforms.
#'
#'
fun_unif_scale <- function(scale_val, Y, az, sig_diag, Tkj, a_seq, b_seq) {
    sig_diag  <- scale_val * sig_diag
    pdiff_mat <- get_pdiffs(resid_vec = c(Y - az), sig_diag = sig_diag,
                            a_seq = a_seq, b_seq = b_seq)
    center_part <- stats::dnorm(Y, mean = az, sd = sqrt(sig_diag))

    log_mat <- log(cbind(pdiff_mat[, 1:length(a_seq)], center_part,
                         pdiff_mat[, (length(a_seq) + 1):ncol(pdiff_mat)]))
    fval <- sum(Tkj * log_mat)
    return(fval)
}

#'Calculates the loglikelihood of the SUCCOTASH model under uniform
#' mixtures.
#'
#'@param pi_Z A vector. The first \code{M} values are the current
#'     values of \eqn{\pi}. The last \code{k} values are the current
#'     values of \eqn{Z}.
#'@param lambda A vector. This is a length \code{M} vector with the
#'     regularization parameters for the mixing proportions.
#'@param alpha A matrix. This is of dimension \code{p} by \code{k} and
#'     are the coefficients to the confounding variables.
#'@param Y A matrix of dimension \code{p} by \code{1}. These are the
#'     observed regression coefficients of the observed variables.
#'@param a_seq A vector of negative numerics containing the left
#'     endpoints of the mixing uniforms.
#'@param b_seq A vector of positiv numerics containing the right
#'     endpoints of the mixing uniforms.
#'@param sig_diag A vector of length \code{p} containing the variances
#'     of the observations.
#'@param var_scale A logical. Should we update the scaling on the
#'     variances (\code{TRUE}) or not (\code{FALSE}).
#' @param likelihood Can be \code{"normal"} or \code{"t"}.
#' @param df A positive numeric. The degrees of freedom if the
#'     likelihood is t.
#'
#'
#'@seealso \code{\link{uniform_succ_given_alpha}}
#'     \code{\link{succotash_unif_fixed}}.
#'
#' @export
succotash_llike_unif <- function(pi_Z, lambda, alpha, Y, a_seq, b_seq, sig_diag,
                                 var_scale = TRUE, likelihood = c("normal", "t"),
                                 df = NULL) {
    M <- length(a_seq) + length(b_seq) + 1

    likelihood = match.arg(likelihood, c("normal", "t"))
    if (likelihood == "normal") {
        df <- 1000
    } else if (likelihood == "t" & is.null(df)) {
        stop("t likelihood specified but df is null")
    }

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
    zero_means <- dt_wrap(x = Y, df = df, mean = az, sd = sqrt(sig_diag))

    left_ispos <- left_means > 0
    right_ispos <- right_means > 0

    pnorm_left_diff <- matrix(NA, ncol = ncol(left_means), nrow = nrow(left_means))
    pnorm_right_diff <- matrix(NA, ncol = ncol(right_means), nrow = nrow(right_means))

    pnorm_left_diff[!left_ispos] <- pt_wrap(x = left_means[!left_ispos], df = df) -
        pt_wrap(x = left_means_zero[!left_ispos], df = df)
    pnorm_right_diff[!right_ispos] <- pt_wrap(x = right_means_zero[!right_ispos], df = df) -
        pt_wrap(x = right_means[!right_ispos], df = df)

    pnorm_left_diff[left_ispos] <- pt_wrap(x = -1 * left_means_zero[left_ispos], df = df) -
        pt_wrap(x = -1 * left_means[left_ispos], df = df)
    pnorm_right_diff[right_ispos] <- pt_wrap(x = -1 * right_means[right_ispos], df = df) -
        pt_wrap(x = -1 * right_means_zero[right_ispos], df = df)

    pnorm_left_diff <- pnorm_left_diff %*% diag(1 / a_seq)
    pnorm_right_diff <- pnorm_right_diff %*% diag(1 / b_seq)

    overall_fmat <- abs(cbind(pnorm_left_diff, zero_means, pnorm_right_diff))

    ## which_lambda <- lambda == 1
    ## if (sum(which_lambda) > 0) {
    ##     prior_weight <- sum(log(pi_current[which_lambda]) *
    ##                         lambda[which_lambda])
    ## } else {
    ##     prior_weight <- 0
    ## }

    llike_new <- sum(log(rowSums(overall_fmat %*% diag(pi_current)))) ## + prior_weight

    return(llike_new)
}


#' Second step of SUCCOTASH with uniform mixture.
#'
#' @param Y A p by 1 matrix of numerics. The data.
#' @param alpha A p by k matrix of numerics. The confounder
#'     coefficients.
#' @param num_em_runs An integer. The number of em iterations to
#'     perform, starting at random locations.
#' @param a_seq A vector of negative numerics in increasing order. The
#'     negative end points in an [a, 0] grid.
#' @param b_seq A vector of positive numerics in increasing order. The
#'     positive end points in a [0, b] grid.
#' @param lambda A vector of numerics greater than or equal to 1, of
#'     length \code{length(a) + length(b) + 1}.
#' @param em_itermax A positive integer. The maximum number of
#'     iterations to perform on the em step.
#' @param em_tol A positive numeric. The stopping criterion for the EM
#'     algorithm.
#' @param pi_init A vector of non-negative numerics that sum of 1 of
#'     length \code{length(a) + length(b) + 1}. The initial values of
#'     the mixture probs.
#' @param Z_init A vector of length k of numerics. Starting values of
#'     Z.
#' @param em_z_start_sd A positive numeric. Z is initialized by iid
#'     normals with this standard deviation and mean 0.
#' @param pi_init_type Either "random", "uniform", or "zero_conc". How
#'     should we choose the initial mixture probabilities if pi_init
#'     is NULL? "random" will draw draw pi uniformly from the
#'     simplex. "uniform" will give each value equal mass. "zero_conc"
#'     will give more mass to 0 than any other probability.
#' @param lambda_type How should we regularize? 'unif' gives no
#'     regularization.  'zero_conf' gives regularization at zero
#'     alone.
#' @param print_progress A logical. Should we plot the progress?
#' @param print_ziter A logical. Should we print the progress of the
#'     Newton iterations for updating Z?
#' @param true_Z The true Z values. Used for testing.
#' @param sig_diag A vector of the variances of \code{Y}.
#' @param var_scale A logical. Should we update the scaling on the
#'     variances (\code{TRUE}) or not (\code{FALSE}).
#' @param df The degrees of freedom of the the t-likelihood if
#'     \code{likleihood = "t"}.
#' @inheritParams succotash
#'
#' @seealso \code{\link{succotash_llike_unif}}
#'     \code{\link{succotash_unif_fixed}}
#'
#' @export
uniform_succ_given_alpha <-
  function(Y, alpha, sig_diag, num_em_runs = 2,
           a_seq = NULL, b_seq = NULL, lambda = NULL,
           em_itermax = 200, em_tol = 10 ^ -6, pi_init = NULL, Z_init = NULL,
           em_z_start_sd = 1, pi_init_type = "random",
           lambda_type = "zero_conc", print_progress = TRUE, print_ziter = FALSE,
           true_Z = NULL, var_scale = TRUE, optmethod = c("coord", "em"),
           likelihood = c("normal", "t"), df = NULL) {
    p <- nrow(Y)
    ## k <- ncol(alpha)

    likelihood = match.arg(likelihood, c("normal", "t"))
    if (likelihood == "normal") {
        df <- 1000
    } else if (likelihood == "t" & is.null(df)) {
        stop("t likelihood specified but df is null")
    } else if (likelihood != "t" & likelihood != "normal") {
        stop("needs to be either a t likelihood or a normal likelihood")
    }

    optmethod <- match.arg(optmethod,  c("coord", "em"))

    ## set up grid
    if (is.null(a_seq)) {
      a_max <- 2 * sqrt(max(Y ^ 2 - sig_diag))
      a_min <- sqrt(min(sig_diag)) / 10
      if (a_max < 0) {
        a_max <- 8 * a_min
      }
      a_current <- a_min
      a_seq <- a_min
      mult_fact <- sqrt(2)
      while (a_current <= a_max) {
        a_current <- a_current * mult_fact
        a_seq <- c(a_seq, a_current)
      }
      a_seq <- sort(-1 * a_seq)
    }
    if (is.null(b_seq)) {
      b_max <- 2 * sqrt(max(Y ^ 2 - sig_diag))
      b_min <- sqrt(min(sig_diag)) / 10
      if (b_max < 0) {
        b_max <- 8 * b_min
      }
      b_current <- b_min
      b_seq <- b_min
      mult_fact <- sqrt(2)
      while (b_current <= b_max) {
        b_current <- b_current * mult_fact
        b_seq <- c(b_seq, b_current)
      }
    }

    M <- length(a_seq) + length(b_seq) + 1 ## plus 1 for pointmass at 0



    if (is.null(lambda)) {
      if (lambda_type == "unif") {
        lambda <- rep(1, M)
      } else if (lambda_type == "zero_conc") {
        lambda <- c(rep(1, length = length(a_seq)), 10, rep(1, length = length(b_seq)))
      }
    }


    ## zero conc for first EM algorithm.
    em_out <- uniform_succ_em(pi_init = pi_init, Z_init = Z_init,
                              a_seq = a_seq, b_seq = b_seq,
                              lambda = lambda, alpha = alpha, Y = Y,
                              sig_diag = sig_diag,
                              print_ziter = print_ziter,
                              print_progress = print_progress,
                              em_itermax = em_itermax,
                              em_tol = em_tol,
                              pi_init_type = "zero_conc",
                              em_z_start_sd = em_z_start_sd,
                              true_Z = true_Z, var_scale = var_scale,
                              optmethod = optmethod,
                              likelihood = likelihood,
                              df = df)

    ## Random init for other EM algorithms.
    if (num_em_runs > 1) {
        for (em_index in 2:num_em_runs) {
            em_new <- uniform_succ_em(pi_init = pi_init, Z_init = Z_init,
                                      a_seq = a_seq, b_seq = b_seq,
                                      lambda = lambda, alpha = alpha,
                                      Y = Y, sig_diag = sig_diag,
                                      print_ziter = print_ziter,
                                      print_progress = print_progress,
                                      em_itermax = em_itermax,
                                      em_tol = em_tol,
                                      pi_init_type = "random",
                                      em_z_start_sd = em_z_start_sd,
                                      true_Z = true_Z, var_scale = var_scale,
                                      optmethod = optmethod,
                                      likelihood = likelihood,
                                      df = df)
            if (em_new$llike > em_out$llike) {
                em_out <- em_new
            }
        }
    }
    pi_current    <- em_out$pi_new
    Z_current     <- em_out$Z_new
    ## llike_current <- em_out$llike
    scale_val     <- em_out$scale_val

    mix_fit <- ashr::unimix(pi = pi_current,
                            a = c(a_seq, rep(0, length(b_seq) + 1)),
                            b = c(rep(0, length(a_seq) + 1), b_seq))

    az <- alpha %*% matrix(Z_current, ncol = 1)

    if (length(df) == 1) {
        betahat <- ashr::postmean(m = mix_fit, betahat = c(Y - az),
                                  sebetahat = sqrt(sig_diag * scale_val),
                                  v = rep(df, p))
        probs <- ashr::comppostprob(m = mix_fit, x = c(Y - az),
                                    s = sqrt(sig_diag * scale_val),
                                    v = rep(df, p))
    } else if (length(df) == p) {
        betahat <- ashr::postmean(m = mix_fit, betahat = c(Y - az),
                                  sebetahat = sqrt(sig_diag * scale_val),
                                  v = df)
        probs <- ashr::comppostprob(m = mix_fit, x = c(Y - az),
                                    s = sqrt(sig_diag * scale_val),
                                    v = df)
    } else {
        stop("df not the correct length")
    }

    lfdr <- probs[length(a_seq) + 1, ]
    qvals <- ashr::qval.from.lfdr(lfdr)

    pi0 <- pi_current[length(a_seq) + 1]


    NegativeProb <- rep(0, length = p)
    if (length(df) == 1) {
        NegativeProb <- ashr::cdf_post(mix_fit, 0,
                                       betahat = c(Y - az),
                                       sebetahat = sqrt(sig_diag * scale_val),
                                       v = rep(df, p)) - lfdr
    } else if (length(df) == p) {
        NegativeProb <- ashr::cdf_post(mix_fit, 0,
                                       betahat = c(Y - az),
                                       sebetahat = sqrt(sig_diag * scale_val),
                                       v = df) - lfdr
    } else {
        stop("df not the correct length")
    }

    lfsr <- ifelse(NegativeProb > 0.5 * (1 - lfdr), 1 - NegativeProb,
                   NegativeProb + lfdr)

    return(list(Z = Z_current, pi_vals = pi_current,
                scale_val = scale_val, a_seq = a_seq, b_seq = b_seq,
                lfdr = lfdr, lfsr = lfsr, betahat = betahat,
                qvals = qvals, pi0 = pi0))
  }

#' EM algorithm for second step of SUCCOTASH
#'
#'
#'
#' @inheritParams uniform_succ_given_alpha
#'
#'
uniform_succ_em <- function(Y, alpha, sig_diag, a_seq, b_seq,
                            pi_init = NULL, Z_init = NULL, lambda = NULL,
                            print_ziter = FALSE, print_progress = FALSE, em_z_start_sd = 1,
                            em_itermax = 200,
                            em_tol = 10 ^ -3,
                            pi_init_type = c("random", "uniform", "zero_conc"), true_Z = NULL,
                            var_scale = TRUE, optmethod = c("coord", "em"),
                            likelihood = c("normal", "t"), df = NULL) {

    optmethod    <- match.arg(optmethod, c("coord", "em"))
    pi_init_type <- match.arg(pi_init_type, c("random", "uniform", "zero_conc"))

    likelihood = match.arg(likelihood, c("normal", "t"))
    if (likelihood == "normal") {
        df <- 1000
    } else if (likelihood == "t" & is.null(df)) {
        stop("t likelihood specified but df is null")
    } else if (likelihood != "t" & likelihood != "normal") {
        stop("needs to be either a t likelihood or a normal likelihood")
    }

    M <- length(a_seq) + length(b_seq) + 1

    assertthat::are_equal(nrow(Y), nrow(alpha))

    p <- nrow(Y)
    k <- ncol(alpha)

    if (is.null(lambda)) {
        lambda <- rep(1, M)
    }

    assertthat::are_equal(length(lambda), M)

    if (is.null(pi_init) | length(pi_init) != M) {
      if (pi_init_type == "random") {
        ## random start points
        pi_init <- rep(NA, length = M)
        temp <- abs(stats::rnorm(M))
        pi_init[1:M] <- temp / sum(temp)
      } else if (pi_init_type == "uniform") {
        ## uniform mass
        pi_init[1:M] <- 1 / M
      } else if (pi_init_type == "zero_conc") {
        ## most mass at 0
        pi_init[1:M] <- min(1 / p, 1 / M)
        pi_init[length(a_seq) + 1] <- 1 - sum(pi_init[2:M])
      } else {
          stop("pi_init_type is bad")
      }
    }

    if (is.null(Z_init)) {
        Z_init <- matrix(stats::rnorm(k, sd = em_z_start_sd), nrow = k)
    } else if (length(Z_init) != k) {
        message("Z_init not correct dimensions.\n Using random starting location")
        Z_init <- matrix(stats::rnorm(k, sd = em_z_start_sd), nrow = k)
    }

    if (var_scale) {
        pi_Z <- c(pi_init, Z_init, 1)
    } else {
        pi_Z <- c(pi_init, Z_init)
    }

    pi_new <- pi_Z[1:M]
    ##plot(c(a_seq, 0, b_seq), pi_new, type = "h", ylab = expression(pi), xlab = "a or b",
    ##     ylim = c(0,1))
    llike_current <- succotash_llike_unif(pi_Z = pi_Z,
                                          lambda = lambda,
                                          alpha = alpha, Y = Y,
                                          a_seq = a_seq,
                                          b_seq = b_seq,
                                          sig_diag = sig_diag,
                                          var_scale = var_scale,
                                          likelihood = likelihood,
                                          df = df)
    ##mtext(side = 3, paste("llike =", round(llike_current)))

    if (optmethod == "coord") {
        coord_out <- fit_succotash_unif_coord(pi_Z = pi_Z, lambda = lambda,
                                             alpha = alpha, Y = Y, a_seq = a_seq,
                                             b_seq = b_seq, sig_diag = sig_diag,
                                             print_ziter = print_ziter,
                                             newt_itermax = em_itermax, tol = em_tol,
                                             var_scale = var_scale, likelihood = likelihood,
                                             df = df)
        pi_Z_new <- coord_out$pi_Z
    } else if (optmethod == "em") {
        if (likelihood == "t") {
            stop("EM too janky for t-likelihood,\nuse optmethod = \"coord\" instead.")
        }
        message("EM selected when coordinate ascent is recommended for uniform mixtures,\nbut here goes nothing!")
        sq_out <- SQUAREM::fpiter(par = pi_Z, lambda = lambda, alpha = alpha,
                                  Y = Y, a_seq = a_seq, b_seq = b_seq,
                                  sig_diag = sig_diag, var_scale = var_scale,
                                  fixptfn = succotash_unif_fixed,
                                  control = list(maxiter = em_itermax, tol = em_tol))
        pi_Z_new <- sq_out$par
    } else {
        stop("no optmethod provided")
    }

    if (var_scale) {
        pi_new    <- pi_Z_new[1:M]
        Z_new     <- pi_Z_new[(M + 1):(length(pi_Z) - 1)]
        scale_val <- pi_Z_new[length(pi_Z)]
    } else {
        pi_new    <- pi_Z_new[1:M]
        Z_new     <- pi_Z_new[(M + 1):length(pi_Z)]
        scale_val <- 1
    }
    llike_current <- succotash_llike_unif(pi_Z = c(pi_new, Z_new, scale_val), lambda = lambda,
                                          alpha = alpha, Y = Y, a_seq = a_seq,
                                          b_seq = b_seq, sig_diag = sig_diag,
                                          var_scale = TRUE, likelihood = likelihood,
                                          df = df)
    ## var_scale = TRUE here b/c scale_val is 1 is var_scale = FALSE for realz

    ## DEFUNCT CODE -----------------------------------------------------------------
    ## em_index <- 1
    ## ldiff <- em_tol + 1
    ## zdiff <- 1
    ## Z_new <- Z_init
    ## while(em_index < em_itermax & ldiff > em_tol ) {
    ##   llike_old <- llike_current
    ##   Z_old <- Z_new
    ##   succ_fixed_out <- succotash_unif_fixed(pi_Z, lambda, alpha, Y, a_seq, b_seq,
    ##                                          sig_diag, print_ziter = print_ziter)
    ##   pi_Z <- succ_fixed_out$pi_Z
    ##   llike_current <- succ_fixed_out$llike_new
    ##   eval_hess <- succ_fixed_out$eval_hess
    ##   pi_new <- pi_Z[1:M]
    ##   Z_new <- pi_Z[(M+1):length(pi_Z)]
    ##   ldiff <- abs(llike_current / llike_old - 1)
    ##   zdiff <- sum(abs(Z_old - Z_new))

    ##   if(print_progress) {
    ##     cat(" Iter =", em_index, "\n")
    ##     cat("ldiff =", ldiff, "\n")
    ##     cat("zdiff =", zdiff, "\n\n")
    ##     par(mgp = c(1.5, 0.5, 0), mar = c(3, 3, 2, 0.5))
    ##     if(!is.null(true_Z)) {
    ##       par(mfrow = c(2,1))
    ##       lm_Z <- lm(Z_new ~ true_Z)
    ##       plot(true_Z, Z_new, xlab = "True Z", ylab = "Estimated Z")
    ##       abline(lm_Z, lty = 2, col = 2)
    ##       abline(0, 1)
    ##       legend("bottomright", c("Y = X", "Regression Line"),
    ##              lty = c(1,2), col = c(1,2))
    ##     }
    ##     plot(c(a_seq, 0, b_seq), pi_new, type = "h", ylab = expression(hat(pi)[k]),
    ##          xlab = "a or b", ylim = c(0,max(c(0.5, pi_new))))
    ##     mtext(side = 3, paste0("Iteration = ", em_index, ", llike = ",
    ##                            round(llike_current), ", Min/Max Eigenvalue = ",
    ##                            format(min(eval_hess), digits = 3), "/",
    ##                            format(max(eval_hess), digits = 3)))
    ##     par(mfrow = c(1,1))
    ##   }
    ##   em_index <- em_index + 1
    ## }
    ## ---------------------------------------------------------------------------

    return(list(pi_new = pi_new, Z_new = Z_new, llike = llike_current,
                scale_val = scale_val))
}
