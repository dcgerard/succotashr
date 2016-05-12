library(succotashr)
context("OLS estimates and SE's are correctly calculated")

test_that("get out ols estimates and ols standard errors from certain versions of SUCCOTASH", {
    set.seed(122)
    ## generate random data for succotash
    p <- 23
    n <- 11
    k <- 5
    q <- 1
    X <- matrix(rnorm(q * n), nrow = n, ncol = q)
    beta <- matrix(rnorm(q * p), nrow = q, ncol = p)
    Z <- matrix(rnorm(k * n), nrow = n, ncol = k)
    alpha <- matrix(rnorm(k * p), nrow = k, ncol = p)
    E <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Y <- X %*% beta + Z %*% alpha + E

    suc_out <- succotash(Y = Y, X = X, k = k, use_ols_se = TRUE, var_scale = FALSE,
                         two_step = FALSE, optmethod = "em")

    lm_out <- limma::lmFit(object = t(Y), design = X)

    ols_sebetahat <- lm_out$stdev.unscaled[, q] * lm_out$sigma

    expect_equal(c(suc_out$Y1_scaled), lm_out$coefficients[, q])
    expect_equal(suc_out$sig_diag_scaled, ols_sebetahat ^ 2)

}
)
