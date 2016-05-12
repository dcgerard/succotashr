library(succotashr)
context("function wrappers")

test_that("dt and pt work", {
    mean_val = 5
    sd_val = 2
    df = 3
    x_seq <- seq(-5, 5, length = 20)

    int_seq <- rep(NA, length = length(x_seq))
    pt_seq  <- rep(NA, length = length(x_seq))
    for (x_index in 1:length(x_seq)) {
        x <- x_seq[x_index]
        int_seq[x_index] <- stats::integrate(f = dt_wrap, lower = -100, upper = x, df = df,
                                             mean = mean_val, sd = sd_val)$value

        pt_seq[x_index] <- pt_wrap(x = x, df = df, mean = mean_val, sd = sd_val)
    }
    expect_equal(int_seq, pt_seq, tol = 10 ^ -5)
}
)
