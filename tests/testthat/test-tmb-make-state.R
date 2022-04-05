library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(lubridate)

test_that("an error is thrown when disease-free state is missing when required", {
  expect_error({
    (flexmodel(
      params = c(a = 0.5, b = 0.25, c = 0.1),
      state = c(X = 1, Y = 2),
      start_date = "2000-01-01",
      end_date = "2000-02-05",
      do_make_state = TRUE
    ) %>% add_rate("X", "Y", ~ (a) * (c) * (X) + (c) * (b) * (Y))
    %>% add_outflow
    %>% update_tmb_indices
    %>% tmb_fun
    )
  }, regexp = "cannot make the initial state because a disease-free state was not supplied")
})

test_that("make state matches classic macpan without state rounding", {
  reset_spec_version()
  tmb_mode()

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
  reset_spec_version()
  tmb_mode()

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
                           max_iters_eig_pow_meth = 200,
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
  compare_sims(r_sim, tmb_sim, na_is_zero = TRUE)
})

test_that('make state matches vax/variant model without hazard intialization', {

  reset_spec_version()
  options(macpan_pfun_method = "grep")
  #r_tmb_comparable()

  base_params <- read_params("PHAC.csv")
  vax_params <- expand_params_vax(
    params = base_params,
    model_type = "twodose"
  )
  model_params <- expand_params_variant(
    vax_params,
    variant_prop = 1e-7,
    variant_advantage = 1.5,
    variant_vax_efficacy_dose1 = 0.3,
    variant_vax_efficacy_dose2 = 0.8
  ) %>% expand_params_S0(1 - 1e-5)

  model = make_vaccination_model(
    params = model_params,
    state = NULL,
    start_date = "2000-01-01", end_date = "2000-01-01",
    do_hazard = TRUE,
    do_hazard_lin = FALSE,
    do_approx_hazard = FALSE,
    do_approx_hazard_lin = FALSE,
    do_make_state = TRUE,
    max_iters_eig_pow_meth = 200,
    do_variant = TRUE)

  expect_equal(
    initial_state_vector(model),
    c(make_state(params = model_params)))
})

test_that('make state matches vax/variant model with realistic parameters', {
  reset_spec_version()
  #r_tmb_comparable()
  options(macpan_pfun_method = "grep")

  # Need to take more than 100 steps for the
  # eigenvector to converge in rExp
  options(MP_rexp_steps_default = 150)

  load("../../inst/testdata/ontario_flex_test.rda")

  start_date = min(params_timevar$Date)
  end_date = start_date

  r_state = make_state(params = model_params)
  mm = make_vaccination_model(
    params = model_params,
    state = r_state,
    start_date = start_date,
    end_date = end_date,
    do_hazard = TRUE,
    do_variant = TRUE
  )
  tmb_state = initial_state_vector(mm)
  expect_equal(c(r_state), c(tmb_state))
})
