library(testthat)
library(McMasterPandemic)
library(ggplot2)

local_edition(3)
## suppress warnings from snapshots

## test-basic with snapshots instead, essentially testing this test method
## better than hardcoding values

params <- read_params("ICU1.csv")
state <- make_state(params = params, type = "ICU1")

test_that("basic examples_snapshot", {
    s0 <- run_sim_range(params, state, nt = 100)
    expect_snapshot(tail(s0, 10))
    s1 <- run_sim(params, state, start_date = "2020-03-01", end_date = "2020-06-01")
    expect_snapshot(tail(s1, 10))
})

test_that("params methods_snapshot", {
    expect_snapshot(summary(params))
})


time_pars <- data.frame(
    Date = c("2020-03-10", "2020-03-25"),
    Symbol = c("beta0", "beta0"),
    Relative_value = c(0.5, 0.1)
)

test_that("time-varying snapshot", {
    resICU_t <- run_sim(params, state,
        start_date = "2020-03-01",
        end_date = "2020-06-01",
        params_timevar = time_pars,
        step_args = list(do_hazard = TRUE)
    )
    expect_snapshot(tail(resICU_t, 10))
    suppressWarnings(print(plot(resICU_t)))
    ## not showing foi because of log / <= 1 filter
    suppressWarnings(print(plot(resICU_t, condense = FALSE, log = TRUE, drop_states = c("t", "S", "R", "E"))))
    ## FIXME: test values!
})


test_that("time-varying with ndt>1", {
    resICU_t2 <- run_sim(params, state,
        start_date = "2020-03-01",
        end_date = "2020-06-01",
        params_timevar = time_pars,
        step_args = list(do_hazard = TRUE),
        ndt = 10
    )
    expect_snapshot(tail(resICU_t2, 10))
})

test_that("ndt>1", {
    s2 <- run_sim_range(params, state, nt = 100, dt = 0.2)
    expect_snapshot(tail(s2, 10))
    s3 <- run_sim(params, state,
        ndt = 20,
        start_date = "2020-03-01",
        end_date = "2020-03-20"
    )
    s3B <- run_sim(params, state,
        start_date = "2020-03-01",
        end_date = "2020-03-20"
    )
    expect_snapshot(tail(s3, 10))
    expect_snapshot(tail(s3B, 10))
})


test_that("var-specific obsdisp", {
    params <- update(params, obs_disp = 1, obs_disp_I = NA, obs_disp_E = 1000)
    set.seed(101)
    ## initial state *without* using eigvec, to match previous reference results
    ss <- make_state(params[["N"]], params[["E0"]], use_eigvec = FALSE)
    s0 <- run_sim(params, stoch = c(obs = TRUE, proc = FALSE), state = ss)
    plot(s0, keep_states = c("I", "E", "report"), log = TRUE)
    expect_snapshot(tail(s0, 10))
    s1 <- run_sim(params, stoch = c(obs = TRUE, proc = FALSE))
    expect_snapshot(tail(s1, 10))
})
