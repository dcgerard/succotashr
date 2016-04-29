library(succotashr)
context("Make sure uniform mixtures run")

test_that("uniform_succotash_em will actually run", {
    set.seed(10)
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


    a_seq <- -4:-1
    b_seq <- 1:4
    
    pzout <- uniform_succ_em(Y = Y, alpha = alpha, sig_diag = sig_diag, a_seq = a_seq,
                             b_seq = b_seq)
    
    expect_true(all(pzout$pi_new > 0))
    expect_true(all(pzout$pi_new < 1))
    expect_equal(sum(pzout$pi_new), 1)
}
)


test_that("succotash_unif_fixed will actually run",{
  set.seed(1200)
  p <- 7
  k <- 2
  lambda = 10
  
  sig_diag <- abs(rnorm(p))
  Z <- matrix(rnorm(k))
  alpha <- matrix(rnorm(p * k), nrow = p)
  Y <- 2 * rnorm(p) + alpha %*% Z
  
  a_seq <- -4:-1
  b_seq <- 1:4
  pi_vals <- abs(rnorm(length(a_seq) + length(b_seq) + 1))
  pi_vals <- pi_vals / sum(pi_vals)
  pi_Z <- c(pi_vals, Z)
  
  
  pzout <- succotash_unif_fixed(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,
                           a_seq = a_seq, b_seq = b_seq, sig_diag = sig_diag)
  expect_equal(sum(pzout[1:length(pi_vals)]), 1)
  expect_true(all(pzout[1:length(pi_vals)] < 1))
  expect_true(all(pzout[1:length(pi_vals)] > 0))
}
)
