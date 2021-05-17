library(testthat)
library(McMasterPandemic)
do_slow <- FALSE
library(dplyr)
params <- fix_pars(read_params("ICU1.csv"))
dd0 <- (ont_all
    %>% trans_state_vars()
    %>% filter(var %in% c("report"))
)

context("calibrate")

opt_pars0 <- c(get_opt_pars(params, vars = "report"))
system.time(cal0 <- suppressWarnings(calibrate(
    data = dd0,
    base_params = params,
    opt_pars = opt_pars0,
    time_args = NULL,
    mle2_control = list(maxit = 5),
    use_DEoptim = FALSE
)))


test_that("truncated calibration 'works'", {
    expect_is(cal0, "fit_pansim")
    expect_equal(
        coef(cal0, "fitted"),
        list(params = c(E0 = 5, beta0 = 0.937671249368027), nb_disp = c(report = 1))
    )
})

dd_nonint <- (dd0 %>%
    mutate(value = value + 0.1))

test_that("non-integer values throw error", {
    expect_error(
        calibrate(
            data = dd_nonint,
            base_params = params,
            opt_pars = opt_pars0,
            time_args = NULL
        ),
        "need integer values"
    )
})

opt_pars <- c(
    get_opt_pars(params),
    list(logit_rel_beta0 = c(-1, -1))
)
dd <- (ont_all
    %>% trans_state_vars()
    %>% filter(var %in% c("report", "death", "H")))
if (do_slow) {
    cal1 <- suppressWarnings(calibrate(
        data = dd,
        base_params = params, opt_pars = opt_pars,
        mle2_control = list(maxit = 5),
        use_DEoptim = FALSE
    ))
}
