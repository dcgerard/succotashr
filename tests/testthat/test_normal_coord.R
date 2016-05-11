library(succotashr)
context("Make sure likelihood increases monotonically")

test_that("normal_coord will run", {
    skip("normal_coord not working yet")
    set.seed(7999)
    p <- 101
    k <- 20
    m <- 10
    lambda <- rep(1, length = m) ## nullpi weight
    scale_val <- 1
    pi_vals <- c(0.5, 0.3, 0.2)
    tau_seq <- seq(0, 4, length = m)
    sig_diag <- abs(rnorm(p))
    plot_new_ests <- FALSE
    Z <- matrix(rnorm(k))
    Z_init <- rnorm(k) / 10
    pi_Z <- c(pi_vals, Z_init)
    var_scale <- FALSE
    alpha <- matrix(rnorm(p * k), nrow = p)
    beta <- draw_beta(pi_vals = pi_vals, tau_seq = tau_seq, p = p)
    Y <- beta  + alpha %*% Z

    coord_out <- normal_coord(pi_Z = pi_Z, lambda = lambda,
                              alpha = alpha, Y = Y, tau_seq = tau_seq,
                              sig_diag = sig_diag,
                              var_scale = var_scale)

    em_out <- succotash_em(Y = Y, alpha = alpha, sig_diag = sig_diag,
                           tau_seq = tau_seq, pi_init = pi_vals,
                           lambda = lambda, Z_init = Z_init,
                           var_scale = var_scale)

    succotash_llike(pi_Z = coord_out$pi_Z, lambda = lambda,
                    alpha = alpha, Y = Y, tau_seq = tau_seq,
                    sig_diag = sig_diag, var_scale = var_scale)

    succotash_llike(pi_Z = c(em_out$pi_vals, em_out$Z),
                    lambda = lambda,
                    alpha = alpha, Y = Y, tau_seq = tau_seq,
                    sig_diag = sig_diag, var_scale = var_scale)


}
)

test_that("normal_llike_simp is same as succotahs_llike", {
    set.seed(22722)
    p <- 1001
    k <- 20
    m <- 3
    lambda <- rep(1, length = m) ## nullpi weight
    scale_val <- 1
    pi_vals <- c(0.5, 0.3, 0.2)
    tau_seq <- seq(0, 4, length = m)
    sig_diag <- abs(rnorm(p))
    plot_new_ests <- FALSE
    Z <- matrix(rnorm(k))
    pi_Z <- c(pi_vals, Z)
    var_scale <- FALSE
    alpha <- matrix(rnorm(p * k), nrow = p)
    beta <- draw_beta(pi_vals = pi_vals, tau_seq = tau_seq, p = p)
    Y <- beta  + alpha %*% Z

    llike1 <- succotash_llike(pi_Z = pi_Z, lambda = lambda,
                              alpha = alpha, Y = Y, tau_seq = tau_seq,
                              sig_diag = sig_diag,
                              var_scale = var_scale)

    llike2 <- normal_llike_simp(Z = Z, Y = Y, alpha = alpha,
                                sig_diag = sig_diag,
                                tau_seq = tau_seq, scale_val = 1,
                                pi_vals = pi_vals)
    expect_equal(llike1, llike2)
}
)
