context("fixed point iteration of scaling variance")

library(succotashr)
context("var_scale tests")

test_that("succotash_fixed will actually run with var_scale = TRUE", {
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


test_that("succotash_em will actually run with var_scale = TRUE", {
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


test_that("succotash_given_alpha will actually run with var_scale = TRUE", {
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


test_that("two-step actually works", {
    set.seed(888)
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

    suc0 <- succotash(Y = Y, X = X, k = k, two_step = FALSE, var_scale = TRUE)
    new_scale <- suc0$scale_val * n / (n - k - q)
    suc1 <- succotash(Y = Y, X = X, k = k, two_step = FALSE, inflate_var = new_scale,
                      var_scale = FALSE)
    suc2 <- succotash(Y = Y, X = X, k = k, two_step = TRUE, var_scale = TRUE)

    expect_equal(suc0$sig_diag_scaled * n / (n - k - q),
                 suc1$sig_diag_scaled)
    expect_equal(suc1$sig_diag_scaled, suc2$sig_diag_scaled, tol = 10^-4)
}
)
