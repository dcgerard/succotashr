library(succotashr)
context("Make sure it actually runs")

test_that("succotash_fixed will actually run with var_scale = TRUE", {
    set.seed(992)
    ## generate random data for succotash
    p <- 23
    n <- 11
    k <- 5
    q <- 2
    X <- matrix(rnorm(q * n), nrow = n, ncol = q)
    beta <- matrix(rnorm(q * p), nrow = q, ncol = p)
    Z <- matrix(rnorm(k * n), nrow = n, ncol = k)
    alpha <- matrix(rnorm(k * p), nrow = k, ncol = p)
    E <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- X %*% beta + Z %*% alpha + E

    suc_out <- succotash(Y = Y, X = X, k = k, fa_method = "pca", optmethod = "em")
    expect_true(suc_out$llike > suc_out$null_llike)

    ## generate random data for succotash_given_alpgha
    beta1 <- matrix(rnorm(p), ncol = 1)
    Z1 <- matrix(rnorm(k), ncol = 1)
    E1 <- matrix(rnorm(p), ncol = 1)
    Y1 <- beta1 + t(alpha) %*% Z1 + E1

    suca <- succotash_given_alpha(Y = Y1, alpha = t(alpha),
                                  sig_diag = rep(1, p),
                                  lambda_type = "zero_conc", lambda0 = 1)

    expect_equal(sum(suc_out$pi_vals), 1)
    expect_true(all(suc_out$pi_vals > -10 ^ -14))
    expect_true(all(suc_out$lfdr > -10 ^ -14))
    expect_true(all(suc_out$qvals > -10 ^ -14))

    expect_true(all(suc_out$pi_vals < 1 + 10 ^ -14))
    expect_true(all(suc_out$lfdr < 1 + 10 ^ -14))
    expect_true(all(suc_out$qvals < 1 + 10 ^ -14))
    expect_true(suca$llike > suca$null_llike)
}
)


test_that("k estimated when not provided", {
    set.seed(491)
    ## generate random data for succotash
    p <- 23
    n <- 11
    k <- 5
    q <- 2
    X <- matrix(rnorm(q * n), nrow = n, ncol = q)
    beta <- matrix(rnorm(q * p), nrow = q, ncol = p)
    Z <- matrix(rnorm(k * n), nrow = n, ncol = k)
    alpha <- matrix(rnorm(k * p), nrow = k, ncol = p)
    E <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- X %*% beta + Z %*% alpha + E

    suc_out <- succotash(Y = Y, X = X, optmethod = "em")
    expect_true(nrow(suc_out$Z) > 0)
}
)


test_that("all null give same likelihood in succotash_llike and succotash_llike_unif", {
    set.seed(92235)
    ## generate random data for succotash_given_alpha
    p <- 23
    n <- 11
    k <- 5
    q <- 2
    beta1 <- matrix(rnorm(p), ncol = 1)
    Z1 <- matrix(rnorm(k), ncol = 1)
    E1 <- matrix(rnorm(p), ncol = 1)
    alpha <- matrix(rnorm(k * p), nrow = k, ncol = p)
    Y1 <- beta1 + t(alpha) %*% Z1 + E1

    llike_norm <- succotash_llike(pi_Z = c(1, Z1, 1), lambda = 1,
                                  alpha = t(alpha), Y = Y1,
                                  tau_seq = 0, sig_diag = rep(1, p),
                                  var_scale = TRUE)

    llike_unif <- succotash_llike_unif(pi_Z = c(1, Z1, 1), lambda = 1,
                                       alpha = t(alpha), Y = Y1,
                                       a_seq = NULL, b_seq = NULL,
                                       sig_diag = rep(1, p),
                                       var_scale = TRUE)

    expect_equal(llike_norm, llike_unif)
}
)
