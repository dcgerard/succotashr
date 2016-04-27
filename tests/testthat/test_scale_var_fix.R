context("fixed point iteration of scaling variance")

library(succotashr)
context("var_scale tests")

test_that("succotash_fixed will actually run with var_scale = TRUE",{
    set.seed(1200)
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

    pzout <- succotash_fixed(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,
                             tau_seq = tau_seq, sig_diag = sig_diag, plot_new_ests = plot_new_ests,
                             var_scale = TRUE)
    expect_equal(length(pzout), m + k + 1)

    pi_Z <- c(pi_vals, Z)
    pzout <- succotash_fixed(pi_Z = pi_Z, lambda = lambda, alpha = alpha, Y = Y,
                             tau_seq = tau_seq, sig_diag = sig_diag, plot_new_ests = plot_new_ests,
                             var_scale = FALSE)
    expect_equal(length(pzout), m + k)
}
)


test_that("succotash_em will actually run with var_scale = TRUE",{
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

    pzout <- succotash_em(Y, alpha, sig_diag, var_scale = TRUE)

    pzout <- succotash_em(Y, alpha, sig_diag, var_scale = FALSE)
}
)


test_that("succotash_given_alpha will actually run with var_scale = TRUE",{
    set.seed(12)
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

    pzout <- succotash_given_alpha(Y, alpha, sig_diag, var_scale = TRUE)

}
)
