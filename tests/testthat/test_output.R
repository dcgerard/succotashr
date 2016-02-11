library(succotashr)
context("Make sure it actually runs")

## generate random data for succotash
p <- 20
n <- 10
k <- 5
q <- 2
X <- matrix(rnorm(q * n), nrow = n, ncol = q)
beta <- matrix(rnorm(q *p), nrow = q, ncol = p)
Z <- matrix(rnorm(k * n), nrow = n, ncol = k)
alpha <- matrix(rnorm(k * p), nrow = k, ncol = p)
E <- matrix(rnorm(n * p), nrow = n, ncol = p)
Y <- X %*% beta + Z %*% alpha + E

suc_out <- succotash(Y = Y, X = X, k = k, fa_method = "quasi_mle")

## generate random data for succotash_given_alpgha
beta1 <- matrix(rnorm(p), ncol = 1)
Z1 <- matrix(rnorm(k), ncol = 1)
E1 <- matrix(rnorm(p), ncol = 1)
Y1 <- beta1 + t(alpha) %*% Z1 + E1

suca <- succotash_given_alpha(Y = Y1, alpha = t(alpha),
                              sig_diag = rep(1, p),
                              lambda_type = "zero_conc", lambda0 = 1)
plot(suca$tau_seq, suca$pi_vals, type = "h")
