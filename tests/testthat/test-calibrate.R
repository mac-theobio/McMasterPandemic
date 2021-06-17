library(testthat)
library(McMasterPandemic)
do_slow <- FALSE
library(dplyr)
params <- fix_pars(read_params("ICU1.csv"))
dd0 <- (ont_all
    %>% trans_state_vars()
    %>% filter(var %in% c("report"))
)

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

test_that("calibrate parameters other than beta0", {
    params <- fix_pars(read_params("ICU1.csv"))
    state <- make_state(params = params)
    time_pars <- data.frame(
        Date = c("2020-05-01", "2020-07-01"),
        Symbol = c("phi2", "beta0"),
        Value = c(1e-3, 1e-3),
        Type = "rel_orig"
    )
    sd <- "2020-03-01"
    ed <- "2020-10-02"
    SS <- run_sim(params,
        state,
        start_date = sd,
        end_date = ed,
        params_timevar = time_pars,
        step_args = list(do_hazard = TRUE)
    )
    plot(SS, log = TRUE, condense = TRUE, log_lwr = -Inf, keep_states = c("I", "death"))
    ## FIXME: would like to be able to pass type="ICU1" through ... ?
    ## (run_sim_break uses make.state(params) without optional args
    pt <- time_pars
    pt$Value <- NA
    SS2 <- run_sim_break(params,
        sim_args = list(
            start_date = sd,
            end_date = ed,
            step_args = list(do_hazard = TRUE)
        ),
        extra_pars = list(time_params = c(1e-3, 1e-3)),
        time_args = list(params_timevar = pt)
    )
    ##
    expect_equal(c(SS), c(SS2))
    ## waldo::compare(c(SS),c(SS2),tolerance=1e-2)
    ## ## run_sim_break has hosp and X (condensed differently?)
    ## library(tidyverse)

    ## list(run_sim=SS, run_sim_break=SS2) %>%
    ##   map_dfr(condense, .id="method") %>%
    ##   select(method,date,death,I,H,ICU) %>%
    ##   pivot_longer(-c(method,date), names_to="var") %>%
    ##   ggplot(aes(x=date,y=value,linetype=method,colour=var)) + geom_line() +
    ##   scale_y_log10()

    ## NOW calibrate:
    dd <- pivot(condense(SS)[, c("date", "death", "report")])
    dd$value <- round(dd$value)
    ggplot(dd, aes(date, value, colour = var)) +
        geom_line() +
        scale_y_log10()
    if (FALSE) {

        ## too slow? not performing very well? but not insane
        calibrate(
            time_args = list(
                break_dates = time_pars$Date,
                Symbol = c("phi2", "beta0")
            ),
            base_params = params,
            data = dd,
            debug_plot = TRUE,
            opt_pars = list(
                params = c(
                    log_E0 = 4,
                    log_beta0 = -1
                ),
                log_value = c(-2.5, -2.5),
                log_nb_disp = NULL
            )
        )
    }
})
