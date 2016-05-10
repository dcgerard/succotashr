library(succotashr)
context("Estimating scale val in uniform mixtures")

test_that("succotash_llike will actually run with var_scale = TRUE",{
    set.seed(12054)
    p <- 7
    k <- 2
    m <- 11
    lambda <- rep(1, length = m) ## nullpi weight
    scale_val <- 1
    pi_vals <- abs(rnorm(m))
    pi_vals <- pi_vals / sum(pi_vals)
    tau_seq <- seq(0, 4, length = m)
    sig_diag <- abs(rnorm(p))
    plot_new_ests <- FALSE
    Z <- matrix(rnorm(k))
    pi_Z <- c(pi_vals, Z, scale_val)
    var_scale <- TRUE
    alpha <- matrix(rnorm(p * k), nrow = p)
    Y <- 2 * rnorm(p) + alpha %*% Z

    a_seq <- seq(-1, 0.1, length = 5)
    b_seq <- seq(0.1, 1, length = 5)

    llike_out <- succotash_llike_unif(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,
                                      a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag,
                                      var_scale = TRUE)

    pi_Z <- c(pi_vals, Z)
    llike_out2 <- succotash_llike_unif(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,
                                      a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag,
                                      var_scale = FALSE)

    expect_equal(llike_out, llike_out2)
}
)


test_that("pnorm diffs are calculated correctly in the functions I split up",{
    set.seed(12054)
    p <- 7
    k <- 2
    m <- 11
    lambda <- rep(1, length = m) ## nullpi weight
    scale_val <- 1
    pi_vals <- abs(rnorm(m))
    pi_vals <- pi_vals / sum(pi_vals)
    tau_seq <- seq(0, 4, length = m)
    sig_diag <- abs(rnorm(p))
    plot_new_ests <- FALSE
    Z <- matrix(rnorm(k))
    pi_Z <- c(pi_vals, Z)
    var_scale <- FALSE
    alpha <- matrix(rnorm(p * k), nrow = p)
    Y <- 2 * rnorm(p) + alpha %*% Z

    a_seq <- seq(-1, -0.1, length = 5)
    b_seq <- seq(0.1, 1, length = 5)

    M <- length(a_seq) + length(b_seq) + 1
    p <- nrow(Y)
    k <- length(pi_Z) - M
    pi_old <- pi_Z[1:M]
    if (k != 0) {
        Z_old <- matrix(pi_Z[(M + 1):(M + k)], nrow = k)
    }

    pdiff_out <- get_pdiffs(resid_vec = c(Y - alpha %*% Z_old), sig_diag = sig_diag,
                            a_seq = a_seq, b_seq = b_seq)


    ## old way of calculating Tkj --------------------------------------------


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

    ## calculate new pi values
    fkj_left <- pnorm_diff_left %*% diag(1 / abs(a_seq))
    fkj_right <- pnorm_diff_right %*% diag(1 / b_seq)
    f0j <- dnorm(Y, az, sqrt(sig_diag))
    fkj <- cbind(fkj_left, f0j, fkj_right)
    fkj_pi <- fkj %*% diag(pi_old)
    Tkj <- diag(1 / rowSums(fkj_pi)) %*% fkj_pi

    llike_old <- sum(log(rowSums(fkj_pi)))

    pi_new <- (colSums(Tkj) + lambda - 1) / (p - M + sum(lambda))
    ## -------------------------------------------------------------

    unif_out <- unif_update_pi(pi_old = pi_old, Y = Y, sig_diag = sig_diag,
                               alpha = alpha, a_seq = a_seq,
                               b_seq = b_seq, Z_old = Z_old, lambda = lambda)

    expect_equal(unif_out$Tkj, Tkj)



    expect_equal(pnorm_diff_left, pdiff_out[, 1:length(a_seq)])
    expect_equal(pnorm_diff_right, pdiff_out[, (length(a_seq) + 1):ncol(pdiff_out)])
}
)



test_that("succotash_unif_fixed will actually run", {
    set.seed(1454)
    p <- 7
    k <- 2
    m <- 11
    lambda <- rep(1, length = m) ## nullpi weight
    scale_val <- 1
    pi_vals <- abs(rnorm(m))
    pi_vals <- pi_vals / sum(pi_vals)
    tau_seq <- seq(0, 4, length = m)
    sig_diag <- abs(rnorm(p))
    plot_new_ests <- FALSE
    Z <- matrix(rnorm(k))
    pi_Z <- c(pi_vals, Z)
    var_scale <- FALSE
    alpha <- matrix(rnorm(p * k), nrow = p)
    Y <- 2 * rnorm(p) + alpha %*% Z

    a_seq <- seq(-1, -0.1, length = 5)
    b_seq <- seq(0.1, 1, length = 5)

    succ_out <- succotash_unif_fixed(pi_Z = pi_Z, lambda = lambda, alpha = alpha,
                                     Y = Y, a_seq = a_seq, b_seq = b_seq,
                                     sig_diag = sig_diag, var_scale = FALSE)
    pi_Z <- c(pi_vals, Z, 1)
    succ_out2 <- succotash_unif_fixed(pi_Z = pi_Z, lambda = lambda, alpha = alpha,
                                     Y = Y, a_seq = a_seq, b_seq = b_seq,
                                     sig_diag = sig_diag, var_scale = TRUE)
    expect_equal(succ_out, succ_out2[1:length(succ_out)], tol = 0.001)
}
)

test_that("uniform_succ_given_alpha will actually run", {
  set.seed(145)
  p <- 100
  k <- 20
  m <- 11
  lambda <- rep(1, length = m) ## nullpi weight
  scale_val <- 1
  pi_vals <- abs(rnorm(m))
  pi_vals <- pi_vals / sum(pi_vals)
  tau_seq <- seq(0, 4, length = m)
  sig_diag <- abs(rnorm(p))
  plot_new_ests <- FALSE
  Z <- matrix(rnorm(k))
  pi_Z <- c(pi_vals, Z)
  var_scale <- FALSE
  alpha <- matrix(rnorm(p * k), nrow = p)
  Y <- 2 * rnorm(p) + alpha %*% Z

  a_seq <- seq(-1, -0.1, length = 5)
  b_seq <- seq(0.1, 1, length = 5)

  succ_out <- uniform_succ_given_alpha(Y = Y, alpha = alpha, sig_diag = sig_diag,
                                       var_scale = TRUE, print_ziter = TRUE, em_itermax = 50)

  ## plot(succ_out$Z, Z)
  ## abline(c(0,1))
}
)


test_that("fit_succotash_unif_coord will actually run", {
  set.seed(245)
  p <- 100
  k <- 10
  m <- 17
  lambda <- rep(1, length = m) ## nullpi weight
  scale_val <- 1
  pi_vals <- abs(rnorm(m))
  pi_vals <- pi_vals / sum(pi_vals)
  tau_seq <- seq(0, 4, length = m)
  sig_diag <- abs(rnorm(p))
  plot_new_ests <- FALSE
  Z <- matrix(rnorm(k))
  Z_init <- rnorm(k) / 10
  pi_Z <- c(pi_vals, Z_init)
  var_scale <- FALSE
  alpha <- matrix(rnorm(p * k), nrow = p)
  Y <- 2 * rnorm(p) + alpha %*% Z

  a_seq <- seq(-10, -0.1, length = 8)
  b_seq <- seq(0.1, 10, length = 8)

  fit_out <- fit_succotash_unif_coord(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,
                                      a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag,
                                      var_scale = var_scale, newt_itermax = 4)
  zhat <- fit_out$pi_Z[(length(pi_vals) + 1):length(fit_out$pi_Z)]
  pihat <- fit_out$pi_Z[1:(length(pi_vals))]
  ## em_out <- uniform_succ_em(Y = Y, alpha = alpha, sig_diag = sig_diag, a_seq = a_seq,
  ##                           b_seq = b_seq, pi_init = pihat, Z_init = zhat,
  ##                           lambda = lambda,
  ##                           var_scale = var_scale)


  ## em_piZ <- c(em_out$pi_new, em_out$Z_new)
  ## succotash_llike_unif(pi_Z = em_piZ, lambda = lambda,
  ##                      Y = Y, alpha = alpha, a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag,
  ##                      var_scale = FALSE)
  ## succotash_llike_unif(pi_Z = fit_out, lambda = lambda,
  ##                      Y = Y, alpha = alpha, a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag,
  ##                      var_scale = FALSE)


  ## plot(Z, zhat)
  ## plot(em_out$Z_new, zhat)
  ## abline(0, 1)

  expect_true(cor(zhat, c(Z)) > 0.9)

  expect_true(all(fit_out$llike_vec[1:(length(fit_out$llike_vec) - 1)] <=
                  fit_out$llike_vec[2:length(fit_out$llike_vec)]))

  pi_Z <- c(pi_Z, 1)
  var_scale <- TRUE

  fit_out2 <- fit_succotash_unif_coord(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,
                                      a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag,
                                      var_scale = var_scale)
  zhat2 <- fit_out2$pi_Z[(length(pi_vals) + 1):(length(fit_out2$pi_Z) - 1)]

  expect_true(all(fit_out2$llike_vec[1:(length(fit_out2$llike_vec) - 1)] <=
                  fit_out2$llike_vec[2:length(fit_out2$llike_vec)]))



  ## em_out <- uniform_succ_em(Y = Y, alpha = alpha, sig_diag = sig_diag, a_seq = a_seq,
  ##                           b_seq = b_seq, pi_init = pi_vals, Z_init = Z, lambda = lambda,
  ##                           var_scale = var_scale)

  ## em_piZ <- c(em_out$pi_new, em_out$Z_new, em_out$scale_val)

  ## fit_out3 <- fit_succotash_unif_coord(pi_Z = em_piZ, lambda = lambda, alpha = alpha, Y = Y,
  ##                                     a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag,
  ##                                     var_scale = var_scale)

  ## succotash_llike_unif(pi_Z = em_piZ, lambda = lambda,
  ##                      Y = Y, alpha = alpha, a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag,
  ##                     var_scale = TRUE)
  ## succotash_llike_unif(pi_Z = fit_out2, lambda = lambda,
  ##                      Y = Y, alpha = alpha, a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag,
  ##                      var_scale = TRUE)

  expect_true(cor(zhat2, Z) > 0.85)
}
)


test_that("two-step actually works for uniform mixtures", {
    set.seed(536)
    p <- 11
    n <- 5
    k <- 1
    q <- 2
    pi_vals <- c(0.5, 0.3, 0.2)
    tau_seq <- c(0, 1, 2)

    beta0 <- matrix(rnorm((q - 1) * p), nrow = q - 1)
    beta1 <- draw_beta(pi_vals = pi_vals, tau_seq = tau_seq, p = p)
    beta  <- rbind(beta0, beta1)
    X     <- matrix(rnorm(n * q), nrow = n)
    Z     <- matrix(rnorm(n * k), nrow = n)
    alpha <- matrix(rnorm(k * p), nrow = k)
    E     <- matrix(rnorm(n * p), nrow = n)
    Y     <- X %*% beta + Z %*% alpha + E

    suc0 <- succotash(Y = Y, X = X, k = k, two_step = FALSE, var_scale = TRUE,
                      mix_type = "uniform")
    new_scale <- suc0$scale_val * n / (n - k - q)
    suc1 <- succotash(Y = Y, X = X, k = k, two_step = FALSE, inflate_var = new_scale,
                      var_scale = FALSE, mix_type = "uniform")
    suc2 <- succotash(Y = Y, X = X, k = k, two_step = TRUE, var_scale = TRUE,
                      mix_type = "uniform")

    expect_equal(suc0$sig_diag_scaled * n / (n - k - q),
                 suc1$sig_diag_scaled)
    expect_equal(suc2$sig_diag_scaled, suc1$sig_diag_scaled, tol = 10 ^ -4)
}
)
