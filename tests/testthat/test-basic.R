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


test_that("time-varying example", {
    time_pars <- data.frame(Date=c("20-Mar-2020","25-Mar-2020"),
                            Symbol=c("beta0","beta0"),
                            Relative_value=c(0.7,0.0))
    resICU_t <- run_sim(params,state,
                        start_date="1-Mar-2020",
                        end_date="1-Jun-2020",
                        params_timevar=time_pars)
    expect_is(resICU_t,"pansim")
})
