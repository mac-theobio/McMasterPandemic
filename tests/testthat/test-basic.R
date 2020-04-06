library(testthat)
library(McMasterPandemic)

context("very basic simulation")

params <- read_params(system.file("params","ICU1.csv",
                                  package="McMasterPandemic"))
state <- make_state(params=params)

test_that("basic examples", {
    expect_is(params,"params_pansim")
    s0 <- run_sim_range(params,state, nt=100)
    expect_is(s0,"data.frame")
    expect_is(state,"state_pansim")
    s1 <- run_sim(params,state,start_date="1-March-2020",end_date="1-Jun-2020")
    expect_is(s1,"pansim")
})

test_that("params methods", {
    expect_equal(summary(params),
                 c(R0 = 6.521938, Gbar = 12.236827, r0 = 0.227825,
                   kappa = 0.45739, kappa_eff = 0.398295, dbl_time = 3.042448),
                 tolerance=1e-6)
})

test_that("time-varying example", {
    time_pars <- data.frame(Date=c("10-Mar-2020","25-Mar-2020"),
                            Symbol=c("beta0","beta0"),
                            Relative_value=c(0.5,0.1))
    resICU_t <- run_sim(params,state,
                        start_date="1-Mar-2020",
                        end_date="1-Jun-2020",
                        params_timevar=time_pars,
                        step_args=list(do_hazard=TRUE))
    expect_is(resICU_t,"pansim")
    plot(resICU_t)
    plot(resICU_t,aggregate=FALSE,log=TRUE,drop_vars=NULL)
})
