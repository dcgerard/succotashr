
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
#'@param Y A matrix of dimension \code{p} by \code{1}. These are the
#'     observed regression coefficients of the observed variables.
#'@param a_seq A vector of negative numerics containing the left
#'     endpoints of the mixing uniforms.
#'@param b_seq A vector of positiv numerics containing the right
#'     endpoints of the mixing uniforms.
#'@param sig_diag A vector of length \code{p} containing the variances
#'     of the observations.
#'@param print_ziter A logical. Should we we print each iteration of
#'     the Z optimization?
#'@param newt_itermax A positive integer. The maximum number of Newton
#'     steps to perform in updating \eqn{Z}.
#'@param tol A positive numeric. The stopping criterion for Newton's
#'     method in updating Z.
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
                                 print_ziter = TRUE, newt_itermax = 100, tol = 10^-4) {
  M <- length(a_seq) + length(b_seq) + 1
  p <- nrow(Y)
  k <- length(pi_Z) - M
  pi_old <- pi_Z[1:M]
  if (k != 0) {
    Z_old <- matrix(pi_Z[(M + 1):(M + k)], nrow = k)
  }

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
    pnorm(-1 * zero_means_left[left_ispos]) - pnorm(-1 * left_means[left_ispos])
  pnorm_diff_right[right_ispos] <-
    pnorm(-1 * right_means[right_ispos]) - pnorm(-1 * zero_means_right[right_ispos])

  pnorm_diff_left[!left_ispos] <-
    pnorm(left_means[!left_ispos]) - pnorm(zero_means_left[!left_ispos])
  pnorm_diff_right[!right_ispos] <-
    pnorm(zero_means_right[!right_ispos]) - pnorm(right_means[!right_ispos])

  ## alternate way to calculate the pnorm_diff's
  #obs_mat_right <- matrix(rep(b_seq, p), nrow = p, byrow = TRUE)
  #mean_mat_right <- matrix(rep(Y-az, length(b_seq)), nrow = p)
  #obs_mat_left <- matrix(rep(a_seq, p), nrow = p, byrow = TRUE)
  #mean_mat_left <- matrix(rep(Y-az, length(a_seq)), nrow = p)
  #sig_left <- matrix(rep(sqrt(sig_diag), length(a_seq)), nrow = p)
  #sig_right <- matrix(rep(sqrt(sig_diag), length(b_seq)), nrow = p)
  #pnorm_diff_r <- pnorm(q = obs_mat_right, mean = mean_mat_right,
  #                      sd = sig_right) -  pnorm(q = 0, mean = mean_mat_right, sd = sig_left)


  ## calculate new pi values
  fkj_left <- pnorm_diff_left %*% diag(1 / abs(a_seq))
  fkj_right <- pnorm_diff_right %*% diag(1 / b_seq)
  f0j <- dnorm(Y, az, sqrt(sig_diag))
  fkj <- cbind(fkj_left, f0j, fkj_right)
  fkj_pi <- fkj %*% diag(pi_old)
  Tkj <- diag(1 / rowSums(fkj_pi)) %*% fkj_pi

  llike_old <- sum(log(rowSums(fkj_pi)))

  pi_new <- (colSums(Tkj) + lambda - 1) / (p - M + sum(lambda))

  z_diff <- tol + 1
  newt_iter <- 1
  Z_new <- Z_old
  mult_val <- 1
  llike_new <- llike_old
  while(z_diff > tol & newt_iter < newt_itermax) {
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
      pnorm(-1 * zero_means_left[left_ispos]) - pnorm(-1 * left_means[left_ispos])
    pnorm_diff_right[right_ispos] <-
      pnorm(-1 * right_means[right_ispos]) - pnorm(-1 * zero_means_right[right_ispos])

    pnorm_diff_left[!left_ispos] <-
      pnorm(left_means[!left_ispos]) - pnorm(zero_means_left[!left_ispos])
    pnorm_diff_right[!right_ispos] <-
      pnorm(zero_means_right[!right_ispos]) - pnorm(right_means[!right_ispos])

    ## calculate gradient
    dnorm_diff_left <- diag(1 / sqrt(sig_diag)) %*%
      (dnorm(zero_means_left) - dnorm(left_means))
    dnorm_diff_right <- diag(1 / sqrt(sig_diag)) %*%
      (dnorm(right_means) - dnorm(zero_means_right))
    dpratios_left <- dnorm_diff_left / pnorm_diff_left
    dpratios_right <- dnorm_diff_right / pnorm_diff_right
    zero_part <- (Y - az) / sig_diag
    alpha_weights <-
      rowSums(Tkj * cbind(dpratios_left, zero_part, dpratios_right))

    gradient_val <- colSums(diag(alpha_weights) %*% alpha)

    ##calculate Hessian
    top_left <- left_means * dnorm(left_means) -
      zero_means_left *  dnorm(zero_means_left)
    top_right <- zero_means_right * dnorm(zero_means_right) -
      right_means * dnorm(right_means)

    sum1 <- top_left / pnorm_diff_left - dpratios_left^2
    sum2 <- top_right / pnorm_diff_right - dpratios_right^2
    sum0 <- -1 / sig_diag

    diag_weights <- rowSums(Tkj * cbind(sum1, sum0, sum2))

    hessian_val <- t(alpha) %*% diag(diag_weights) %*% alpha

    Z_new <- Z_old -  mult_val * solve(hessian_val) %*% gradient_val

    z_diff <- sum(abs(Z_new - Z_old))

    llike_new <- succotash_llike_unif(c(pi_new, Z_new), lambda, alpha, Y, a_seq, b_seq, sig_diag)

    if(llike_new < llike_old) {
      ## halve the step size and reset to old Z value
      mult_val <- mult_val / 2
      Z_new <- Z_old
      llike_new <- llike_old
      did_I_move <- FALSE
    }

    if(print_ziter) {
      cat("Iteration =", newt_iter, "\n")
      cat("Did I move?", did_I_move, "\n")
      cat("Z Diff = ", z_diff, "\n")
      cat("llike = ", llike_new, '\n\n')
      newt_iter <- newt_iter + 1
    }
  }

  ## check to see if Hessian is negative definite at final Z val.
  eval_hess <- eigen(hessian_val, only.values = TRUE)$values

  pi_Z <- c(pi_new, Z_new)

  return(list(pi_Z = c(pi_new, Z_new), llike_new = llike_new, eval_hess = eval_hess))
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
#'
#'@seealso \code{\link{uniform_succ_given_alpha}}
#'     \code{\link{succotash_unif_fixed}}.
#'
#' @export
succotash_llike_unif <- function(pi_Z, lambda, alpha, Y, a_seq, b_seq, sig_diag) {
  M <- length(a_seq) + length(b_seq) + 1
  p <- nrow(Y)
  k <- length(pi_Z) - M
  pi_current <- pi_Z[1:M]
  if (k != 0) {
    Z_current <- matrix(pi_Z[(M + 1):(M + k)], nrow = k)
  }

  az <- alpha %*% Z_current

  left_means <- diag(1 / sqrt(sig_diag)) %*% outer(c((Y - az)), a_seq, "-")
  left_means_zero <- diag(1 / sqrt(sig_diag)) %*% outer(c((Y - az)), rep(0, length(a_seq)), "-")
  right_means <- diag(1 / sqrt(sig_diag)) %*% outer(c((Y - az)), b_seq, "-")
  right_means_zero <- diag(1 / sqrt(sig_diag)) %*% outer(c((Y - az)), rep(0, length(b_seq)), "-")
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

  pnorm_left_diff <- pnorm_left_diff %*% diag(1 / a_seq)
  pnorm_right_diff <- pnorm_right_diff %*% diag(1 / b_seq)

  overall_fmat <- abs(cbind(pnorm_left_diff, zero_means, pnorm_right_diff))

  llike_new <- sum(log(rowSums(overall_fmat %*% diag(pi_current))))

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
#' @param use_SQUAREM A logical. Should we use SQUAREM to run the EM?
#'     Not quite implemented yet.
#' @param sig_diag A vector of the variances of \code{Y}.
#'
#' @seealso \code{\link{succotash_llike_unif}}
#'     \code{\link{succotash_unif_fixed}}
#'
#' @exportn
uniform_succ_given_alpha <-
  function(Y, alpha, sig_diag, num_em_runs = 2,
           a_seq = NULL, b_seq = NULL, lambda = NULL,
           em_itermax = 200, em_tol = 10 ^ -6, pi_init = NULL, Z_init = NULL,
           em_z_start_sd = 1, pi_init_type = "random",
           lambda_type = "zero_conc", print_progress = TRUE, print_ziter = FALSE,
           true_Z = NULL, use_SQUAREM = FALSE) {
    p <- nrow(Y)
    k <- ncol(alpha)

    ## set up grid
    if(is.null(a_seq))
    {
      a_max <- 2 * sqrt(max(Y^2 - sig_diag))
      a_min <- sqrt(min(sig_diag)) / 10
      if(a_max < 0) {
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
    if(is.null(b_seq))
    {
      b_max <- 2 * sqrt(max(Y^2 - sig_diag))
      b_min <- sqrt(min(sig_diag)) / 10
      if(b_max < 0) {
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
      if (lambda_type == 'unif') {
        lambda <- rep(1, M)
      } else if (lambda_type == 'zero_conc') {
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
                              true_Z = true_Z)
    pi_current <- em_out$pi_new
    Z_current <- em_out$Z_new
    llike_current <- em_out$llike

    ## Random init for other EM algorithms.
    if(num_em_runs > 1) {
        for(em_index in 2:num_em_runs) {
            em_out <- uniform_succ_em(pi_init = pi_init, Z_init = Z_init,
                                      a_seq = a_seq, b_seq = b_seq,
                                      lambda = lambda, alpha = alpha,
                                      Y = Y, sig_diag = sig_diag,
                                      print_ziter = print_ziter,
                                      print_progress = print_progress,
                                      em_itermax = em_itermax,
                                      em_tol = em_tol,
                                      pi_init_type = "random",
                                      em_z_start_sd = em_z_start_sd,
                                      true_Z = true_Z)
            pi_new <- em_out$pi_new
            Z_new <- em_out$Z_new
            llike_new <- em_out$llike

            if(llike_new > llike_current) {
                pi_current <- pi_new
                Z_current <- Z_new
                llike_current <- llike_new
            }
        }
    }

    mix_fit <- ashr::unimix(pi = pi_current, a = c(a_seq, rep(0, length(b_seq) + 1)),
                            b = c(rep(0, length(a_seq) + 1), b_seq))

    az <- alpha %*% matrix(Z_current, ncol = 1)

    betahat <- ashr::postmean(m = mix_fit, betahat = c(Y - az), sebetahat = sqrt(sig_diag),
                              v = rep(1000, p))

    probs <- ashr::comppostprob(m = mix_fit, x = c(Y - az), s = sqrt(sig_diag), v = rep(1000, p))
    lfdr <- probs[length(a_seq) + 1,]
    qvals <- ashr::qval.from.lfdr(lfdr)

    pi0 <- pi_current[length(a_seq) + 1]

    NegativeProb = rep(0, length = p)
    NegativeProb <-
      ashr::cdf_post(mix_fit, 0,
               betahat = c(Y - az),
               sebetahat = sqrt(sig_diag), v = rep(1000, p)) -
      lfdr
    lfsr = ifelse(NegativeProb > 0.5 * (1 - lfdr), 1 - NegativeProb,
                  NegativeProb + lfdr)


    #        sq_out <-
    #            SQUAREM::fpiter(par = pi_Z, lambda = lambda, alpha = alpha,
    #                            Y = Y, a_seq = a_seq, b_seq = b_seq,
    #                            sig_diag = sig_diag,
    #                            fixptfn = succotash_unif_fixed)

    return(list(Z = Z_current, pi_vals = pi_current, a_seq = a_seq, b_seq = b_seq,
                lfdr = lfdr, lfsr = lfsr, betahat = betahat, qvals = qvals, pi0 = pi0))
  }

#' EM algorithm for second step of SUCCOTASH
#'
#'
#'
#' @inheritParams uniform_succ_given_alpha
#'
#'
uniform_succ_em <- function(pi_init, Z_init, a_seq, b_seq, lambda, alpha, Y, sig_diag,
                            print_ziter, print_progress, em_z_start_sd, em_itermax = 200,
                            em_tol = 10^-6,
                            pi_init_type = "random", true_Z = NULL) {

    M <- length(a_seq) + length(b_seq) + 1

    p <- nrow(Y)
    k <- ncol(alpha)

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
        pi_init[1:M] <- min(1 / p, 1 / M)
        pi_init[length(a_seq) + 1] <- 1 - sum(pi_init[2:M])
      } else {
          stop("pi_init_type is bad")
      }
    }

    if (is.null(Z_init)) {
      Z_init <- matrix(rnorm(k, sd = em_z_start_sd), nrow = k)
    } else if (length(Z_init) != k) {
      Z_init <- matrix(rnorm(k, sd = em_z_start_sd), nrow = k)
    }


    pi_Z <- c(pi_init, Z_init)

    pi_new <- pi_Z[1:M]
    plot(c(a_seq, 0, b_seq), pi_new, type = "h", ylab = expression(pi), xlab = "a or b",
         ylim = c(0,1))
    llike_current <- succotash_llike_unif(pi_Z, lambda, alpha, Y, a_seq, b_seq, sig_diag)
    mtext(side = 3, paste("llike =", round(llike_current)))

    em_index <- 1
    ldiff <- em_tol + 1
    zdiff <- 1
    Z_new <- Z_init
    while(em_index < em_itermax & ldiff > em_tol ) {
      llike_old <- llike_current
      Z_old <- Z_new
      succ_fixed_out <- succotash_unif_fixed(pi_Z, lambda, alpha, Y, a_seq, b_seq,
                                             sig_diag, print_ziter = print_ziter)
      pi_Z <- succ_fixed_out$pi_Z
      llike_current <- succ_fixed_out$llike_new
      eval_hess <- succ_fixed_out$eval_hess
      pi_new <- pi_Z[1:M]
      Z_new <- pi_Z[(M+1):length(pi_Z)]
      ldiff <- abs(llike_current / llike_old - 1)
      zdiff <- sum(abs(Z_old - Z_new))

      if(print_progress) {
        cat(" Iter =", em_index, "\n")
        cat("ldiff =", ldiff, "\n")
        cat("zdiff =", zdiff, "\n\n")
        par(mgp = c(1.5, 0.5, 0), mar = c(3, 3, 2, 0.5))
        if(!is.null(true_Z)) {
          par(mfrow = c(2,1))
          lm_Z <- lm(Z_new ~ true_Z)
          plot(true_Z, Z_new, xlab = "True Z", ylab = "Estimated Z")
          abline(lm_Z, lty = 2, col = 2)
          abline(0, 1)
          legend("bottomright", c("Y = X", "Regression Line"),
                 lty = c(1,2), col = c(1,2))
        }
        plot(c(a_seq, 0, b_seq), pi_new, type = "h", ylab = expression(hat(pi)[k]),
             xlab = "a or b", ylim = c(0,max(c(0.5, pi_new))))
        mtext(side = 3, paste0("Iteration = ", em_index, ", llike = ",
                               round(llike_current), ", Min/Max Eigenvalue = ",
                               format(min(eval_hess), digits = 3), "/",
                               format(max(eval_hess), digits = 3)))
        par(mfrow = c(1,1))
      }
      em_index <- em_index + 1
    }

    return(list(pi_new = pi_new, Z_new = Z_new, llike = llike_current))
}


## p <- 1000
## k <- 10

## Z <- matrix(rnorm(k), ncol = 1)
## beta <- c(matrix(rnorm(p / 2, mean = 1), ncol = 1), rep(0, p/2))
## alpha <- matrix(rnorm(p * k), nrow = p)
## sig_diag <- rep(1, length = p)
## df <- 2
## Y <- matrix(rt(p, df = df) + beta + alpha %*% Z, ncol = 1)

## em_tol = 10^-6
## em_itermax = 1000
## num_em_runs = 1
## em_z_start_sd = 1
