library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)

test_that("make state matches classic macpan without state rounding", {
  set_spec_version('0.1.1', '../../inst/tmb')
  options(MP_use_state_rounding = FALSE)

  params <- ("ICU1.csv"
    %>% read_params
    %>% expand_params_S0(1 - 1e-5)
  )
  state = make_state(params = params)
  state[] = 0 # make sure that we have to remake this on the c++ side
  model = make_base_model(
    params = params,
    state = state,
    start_date = "1900-01-01", end_date = "1900-01-01",
    do_make_state = TRUE)

  r_init_state = c(make_state(params = params))
  tmb_init_state = initial_state_vector(model)
  expect_equal(r_init_state, tmb_init_state)
})

test_that('make state works with a time-varying parameter', {
  # TODO: to make this pass, we need to engage R-based make_state ... I think
  #       https://github.com/mac-theobio/McMasterPandemic/issues/124

  set_spec_version('0.1.1', '../../inst/tmb')
  options(MP_use_state_rounding = FALSE)

  params <- ("ICU1.csv"
             %>% read_params
             %>% expand_params_S0(1 - 1e-5)
  )
  state = make_state(params = params)
  state[] = 0 # make sure that we have to remake this on the c++ side
  start_date = "2021-05-10"
  end_date = "2021-12-10"
  tv_dat <- data.frame(
    Date = c("2021-06-18", "2021-07-01", "2021-07-25"),
    Symbol = c("beta0", "beta0", "beta0"),
    Value = c(0.1, NA, 0.9),
    Type = c("rel_orig", "rel_orig", "rel_orig")
  )
  model <- make_base_model(params, state,
                           start_date = start_date, end_date = end_date,
                           params_timevar = tv_dat,
                           do_hazard = TRUE,
                           do_make_state = TRUE,
                           max_iters_eig_pow_meth = 100,
                           tol_eig_pow_meth = 1e-03
  )
  tv_dat_filled = tv_dat
  tv_dat_filled$Value = model$timevar$piece_wise$schedule$init_tv_mult

  r_init_state = c(make_state(params = params))
  tmb_init_state = c(initial_state_vector(model))
  expect_equal(tmb_init_state, r_init_state)

  sim_time_comparison = time_wrap(
    tmb_sim <- run_sim(
      params = params, #state = model$state,# make_state(params = params),
      start_date = start_date,
      end_date = end_date,
      params_timevar = tv_dat_filled,
      step_args = list(do_hazard = TRUE),
      condense = TRUE,
      use_flex = TRUE,
      flexmodel = model
    ),
    r_sim <- run_sim(
      params = params, state = make_state(params = params),
      start_date = start_date,
      end_date = end_date,
      params_timevar = tv_dat_filled,
      step_args = list(do_hazard = TRUE),
      condense = TRUE,
      use_flex = FALSE
    )
  )
  compare_sims(r_sim, tmb_sim)
})
