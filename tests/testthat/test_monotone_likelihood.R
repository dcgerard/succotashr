library(succotashr)
context("Make sure likelihood increases monotonically")

test_that("succotash_fixed will increases likelihood",{
    set.seed(692)
    p <- 7
    k <- 2
    m <- 11
    lambda <- 10 ## nullpi weight
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

    itermax <- 15
    llike_vec <- rep(NA, length = itermax)

    llike_vec[1] <- succotash_llike(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,
                                        tau_seq = tau_seq, sig_diag = sig_diag,
                                        plot_new_ests = plot_new_ests,
                                        var_scale = TRUE)

    for (index in 2:itermax) {
        pi_Z <- succotash_fixed(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,
                                 tau_seq = tau_seq, sig_diag = sig_diag,
                                 plot_new_ests = plot_new_ests,
                                 var_scale = TRUE)
        llike_vec[index] <- succotash_llike(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,
                                            tau_seq = tau_seq, sig_diag = sig_diag,
                                            plot_new_ests = plot_new_ests,
                                            var_scale = TRUE)
    }

    expect_true(all(llike_vec[2:length(llike_vec)] > llike_vec[1:(length(llike_vec) - 1)]))
}
)


test_that("succotash_unif_fixed will increase likelihood",{
  set.seed(211)
  p <- 7
  k <- 2

  sig_diag <- abs(rnorm(p))
  Z <- matrix(rnorm(k))
  alpha <- matrix(rnorm(p * k), nrow = p)
  Y <- 2 * rnorm(p) + alpha %*% Z

  a_seq <- -4:-1
  b_seq <- 1:4
  pi_vals <- abs(rnorm(length(a_seq) + length(b_seq) + 1))
  pi_vals <- pi_vals / sum(pi_vals)
  pi_Z <- c(pi_vals, Z, 1)

  lambda <- rep(1, length = length(pi_vals))


  itermax <- 15
  llike_vec <- rep(NA, length = itermax)
  llike_vec[1] <- succotash_llike_unif(pi_Z = pi_Z, lambda = lambda,
                                       alpha = alpha, Y = Y,
                                       a_seq = a_seq, b_seq = b_seq,
                                       sig_diag = sig_diag,
                                       var_scale = TRUE)

  for (index in 2:itermax) {
      pi_Z <- succotash_unif_fixed(pi_Z = pi_Z, lambda = lambda,
                                   alpha = alpha, Y = Y,
                                   a_seq = a_seq, b_seq = b_seq,
                                   sig_diag = sig_diag,
                                   var_scale = TRUE)
      llike_vec[index] <- succotash_llike_unif(pi_Z = pi_Z,
                                               lambda = lambda,
                                               alpha = alpha, Y = Y,
                                               a_seq = a_seq,
                                               b_seq = b_seq,
                                               sig_diag = sig_diag,
                                               var_scale = TRUE)
  }
  expect_true(all(llike_vec[2:length(llike_vec)] > llike_vec[1:(length(llike_vec) - 1)]))

}
)



test_that("llike_unif_simp provides same likelihood values as succotash_llike_unif", {
    set.seed(775)
    p <- 100
    k <- 20

    sig_diag <- abs(rnorm(p))
    Z        <- matrix(rnorm(k))
    alpha    <- matrix(rnorm(p * k), nrow = p)
    Y        <- 2 * rnorm(p) + alpha %*% Z

    a_seq <- -4:-1
    b_seq <- 1:4
    left_seq  <- c(a_seq, rep(0, length(b_seq) + 1))
    right_seq <- c(rep(0, length(a_seq) + 1), b_seq)
    pi_vals   <- abs(rnorm(length(left_seq)))
    pi_vals   <- pi_vals / sum(pi_vals)
    scale_val <- 1

    llike1 <- llike_unif_simp(Y = Y, Z = Z, pi_vals = pi_vals, alpha = alpha,
                               sig_diag = sig_diag,
                               left_seq = left_seq, right_seq = right_seq,
                               scale_val = scale_val)


    pi_Z <- c(pi_vals, c(Z), scale_val)
    lambda <- rep(1, length = length(pi_vals))


    llike2 <- succotash_llike_unif(pi_Z = pi_Z, lambda = lambda, alpha = alpha,
                                   Y = Y, a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag,
                                   var_scale = TRUE)

    expect_equal(llike1, llike2)

}
)



test_that("unif_grad_simp provides same gradient values as unif_grad_llike", {
    skip(message = "unif_grad_llike is broken and should not be used")
    set.seed(775)
    p <- 20
    k <- 5

    sig_diag <- abs(rnorm(p))
    Z        <- matrix(rnorm(k))
    alpha    <- matrix(rnorm(p * k), nrow = p)
    Y        <- 2 * rnorm(p) + alpha %*% Z

    a_seq <- -4:-1
    b_seq <- 1:4
    left_seq  <- c(a_seq, rep(0, length(b_seq) + 1))
    right_seq <- c(rep(0, length(a_seq) + 1), b_seq)
    pi_vals   <- abs(rnorm(length(left_seq)))
    pi_vals   <- pi_vals / sum(pi_vals)
    scale_val <- 1

    pi_Z <- c(pi_vals, c(Z), scale_val)
    lambda <- rep(1, length = length(pi_vals))

    grad1 <- unif_grad_simp(pi_Z = pi_Z, lambda = lambda, alpha = alpha,
                            Y = Y, left_seq = left_seq, right_seq = right_seq,
                            sig_diag = sig_diag, var_scale = TRUE)

    grad2 <- unif_grad_llike(pi_Z = pi_Z, lambda = lambda, alpha = alpha,
                             Y = Y, a_seq = a_seq, b_seq = b_seq,
                             sig_diag = sig_diag, var_scale = TRUE)

    expect_equal(grad1, grad2)

}
)
