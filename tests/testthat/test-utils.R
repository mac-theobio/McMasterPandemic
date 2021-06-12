library(McMasterPandemic)
library(testthat)

context("invlink_trans")

test_that("invlink_trans", {
    ## recapitulating invlink_trans examples
    expect_equal(
        invlink_trans(c(log_p1 = 0, logit_p2 = 0)),
        c(p1 = 1, p2 = 0.5)
    )
    expect_equal(
        invlink_trans(list(log_p1 = c(0, 0), logit_p2 = c(0, 0, 0))),
        list(p1 = c(1, 1), p2 = rep(0.5, 3))
    )
    expect_equal(
        invlink_trans(list(p1 = c(log_a = 0, log_b = 0), p2 = 4)),
        list(p1 = c(a = 1, b = 1), p2 = 4)
    )
    tst <- list(
        params = c(log_beta0 = 0.693147180559945),
        log_nb_disp = 4.60517018598809
    )
    expect_equal(
        invlink_trans(tst),
        list(params = c(beta0 = 2), nb_disp = 100)
    )
    tst2 <- list(
        params = c(log_E0 = 1, log_beta0 = 1),
        log_nb_disp = c(H = 0, report = 0, death = 0)
    )
    expect_equal(
        invlink_trans(tst2),
        list(params = c(E0 = 2.71828182845905, beta0 = 2.71828182845905), nb_disp = c(H = 1, report = 1, death = 1))
    )

    expect_equal(
        invlink_trans(list(time_beta = 1, log_time_beta = 0)),
        list(time_beta = 1, time_beta = 1)
    )
    expect_equal(
        invlink_trans(list(time_beta = c(2, 3), log_time_beta = c(a = 0, b = 0))),
        list(time_beta = c(2, 3), time_beta = c(a = 1, b = 1))
    )
})
