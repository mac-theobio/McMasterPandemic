Sys.setenv(R_TESTS="")

library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)
library(lubridate)
library(tidyr)

testLevel <- if (nzchar(s <- Sys.getenv("MACPAN_TEST_LEVEL"))) as.numeric(s) else 1
skip_slow_tests = isTRUE(testLevel == 1)

test_that('time variation works for several parameters and levels of continuity', {
  test_fun = function(tv_dat) {
    reset_spec_version()
    r_tmb_comparable()
    params <- read_params("ICU1.csv")
    mm = make_base_model(
      expand_params_S0(params, 1-1e-5),
      start_date = "2021-09-10",
      end_date = "2021-10-10",
      params_timevar = tv_dat,
      step_args = list(do_hazard = TRUE)
    )

    # change beta0, which is time-varying, so that
    # we can check that the parameter updates are
    # happening correctly in the C++ side
    test_pars = expand_params_S0(params, 1-1e-5)
    test_pars[1] = 1

    tmb_sim <- run_sim(
      params = test_pars,
      start_date = mm$start_date,
      end_date = mm$end_date,
      params_timevar = tv_dat,
      condense = FALSE,
      step_args = list(do_hazard = TRUE),
      flexmodel = mm
    )

    r_sim <- run_sim(
      params = test_pars,
      start_date = mm$start_date,
      end_date = mm$end_date,
      params_timevar = tv_dat,
      condense = FALSE,
      step_args = list(do_hazard = TRUE)
    )
    print("tmb")
    print(attributes(tmb_sim)$params_timevar)
    print("r")
    print(attributes(r_sim)$params_timevar)
    compare_sims(r_sim, tmb_sim, compare_attr = FALSE)
  }

  test_fun(data.frame(
    Date = as.Date(c("2021-09-15", "2021-09-20", "2021-10-05", "2021-09-20", "2021-10-03")),
    Symbol = c("beta0", "beta0", "beta0", "Ca", "Ca"),
    Value = c(0.5, 0.1, 0.05, 0.1, 2.3),
    Type = c("rel_orig", "rel_orig", "rel_orig", "rel_orig", "rel_orig")
  ))

  test_fun(data.frame(
    Date = ymd(20210910) + days(0:1),
    Symbol = "beta0",
    Value = seq(from = 0.1, to = 2, length.out = 2),
    Type = "rel_orig"
  ))

  test_fun(data.frame(
   Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
   Symbol = c("beta0", "beta0", "beta0"),
   Value = c(0.5, 0.1, 0.05),
   Type = c("rel_orig", "rel_orig", "rel_orig")
  ))
})

test_that("time variation on the first two steps matches in r and tmb engines", {
  reset_spec_version()
  r_tmb_comparable()
  params <- read_params("PHAC.csv") %>% expand_params_S0(1-1e-5)

  test_fun = function(type) {

    tv_dat = data.frame(
      Date = c("2021-09-10", "2021-09-11"),
      Symbol = "beta0",
      Value = c(0.1, 0.2),
      Type = type
    )

    mm = make_base_model(
      params,
      start_date = "2021-09-10",
      end_date = "2021-09-12",
      params_timevar = tv_dat,
      step_args = list(do_hazard = TRUE)
    )

    params[1] = 0.8

    tmb_sim <- run_sim(
      params = params,
      start_date = mm$start_date,
      end_date = mm$end_date,
      condense = FALSE,
      params_timevar = tv_dat,
      step_args = list(do_hazard = TRUE),
      flexmodel = mm
    )

    r_sim <- run_sim(
      params = params,
      start_date = mm$start_date,
      end_date = mm$end_date,
      condense = FALSE,
      params_timevar = tv_dat,
      step_args = list(do_hazard = TRUE),
    )

    attributes(r_sim) = NULL
    attributes(tmb_sim) = NULL
    compare_sims(r_sim, tmb_sim, compare_attr = FALSE)
  }

  sapply(c('rel_orig', 'rel_prev', 'abs'), test_fun)
})

test_that('time variation works for vax models', {
  skip_if(skip_slow_tests)

  reset_spec_version()
  #r_tmb_comparable()
  options(macpan_pfun_method = "grep")
  options(MP_rexp_steps_default = 150)

  load(system.file('testdata', 'ontario_flex_test_better.Rdata', package = 'McMasterPandemic'))

  start_date = min(params_timevar$Date)
  end_date = max(params_timevar$Date) + days(1)

  r_sim <- run_sim(
    params = model_params,
    start_date = start_date,
    end_date = end_date,
    params_timevar = params_timevar,
    condense = FALSE,
    step_args = list(do_hazard = TRUE)
  )
  mm = make_vaccination_model(
    params = model_params,
    state = make_state(params = model_params),
    start_date = start_date,
    end_date = end_date,
    params_timevar = params_timevar,
    do_hazard = TRUE,
    do_variant = TRUE
  )

  test_fun = function(do_make_state, tolerance = testthat_tolerance()) {
    mm$do_make_state = do_make_state
    tmb_sim <- run_sim(
      params = unlist(model_params),
      state = mm$state,
      start_date = mm$start_date,
      end_date = mm$end_date,
      params_timevar = params_timevar,
      condense = FALSE,
      step_args = list(do_hazard = TRUE),
      flexmodel = mm
    )

    compare_sims(r_sim, tmb_sim, tolerance, compare_attr = FALSE)
  }
  test_fun(FALSE)
  test_fun(TRUE)
})

test_that('time variation works for a mix of types', {
  reset_spec_version()
  r_tmb_comparable()
  test_fun = function(tv_dat) {
    params <- read_params("ICU1.csv")
    mm = make_base_model(
      expand_params_S0(params, 1-1e-5),
      start_date = "2021-09-10",
      end_date = "2021-10-10",
      params_timevar = tv_dat,
      step_args = list(do_hazard = TRUE)
    )

    # change beta0, which is time-varying, so that
    # we can check that the parameter updates are
    # happening correctly in the C++ side
    test_pars = expand_params_S0(params, 1-1e-5)
    test_pars[1] = 1

    tmb_sim <- run_sim(
      params = test_pars,
      start_date = mm$start_date,
      end_date = mm$end_date,
      params_timevar = tv_dat,
      condense = FALSE,
      step_args = list(do_hazard = TRUE),
      flexmodel = mm
    )

    r_sim <- run_sim(
      params = test_pars,
      start_date = mm$start_date,
      end_date = mm$end_date,
      params_timevar = tv_dat,
      condense = FALSE,
      step_args = list(do_hazard = TRUE)
    )

    compare_sims(r_sim, tmb_sim, compare_attr = FALSE)
  }

  data.frame(
    Date = as.Date(c("2021-09-15", "2021-09-20", "2021-10-05", "2021-09-20", "2021-10-03")),
    Symbol = c("beta0", "beta0", "beta0", "Ca", "Ca"),
    Value = c(0.5, 0.1, 0.05, 0.1, 2.3),
    Type = c("rel_orig", "rel_prev", "rel_orig", "abs", "rel_prev")
  ) %>% test_fun

  data.frame(
    Date = ymd(20210910) + days(0:1),
    Symbol = "beta0",
    Value = seq(from = 0.1, to = 2, length.out = 2),
    Type = "rel_prev"
  ) %>% test_fun

  data.frame(
    Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
    Symbol = c("beta0", "beta0", "beta0"),
    Value = c(0.5, 0.1, 0.05),
    Type = c("rel_prev", "abs", "rel_orig")
  ) %>% test_fun
})

test_that("breakpoints outside of the simulation range cause a warning", {
  reset_spec_version()
  #r_tmb_comparable()
  options(macpan_pfun_method = "grep")

  start_date <- "2021-02-02"
  end_date <- "2021-02-05"

  params_timevar = data.frame(
    Date = as.Date(c("2021-02-02", "2021-02-03", "2021-02-05")),
    Symbol = c("beta0", "beta0", "beta0"),
    Value = c(0.5, 0.1, 2.2),
    Type = c("rel_orig")
  ) %>% arrange(Symbol, Date)

  base_params <- read_params("PHAC.csv")
  vax_params <- expand_params_vax(
    params = base_params,
    vax_doses_per_day = 1e4,
    model_type = "twodose"
  )
  model_params <- expand_params_variant(
    vax_params,
    variant_prop = 0.5,
    variant_advantage = 1.5,
    variant_vax_efficacy_dose1 = 0.3,
    variant_vax_efficacy_dose2 = 0.8
  ) %>% expand_params_S0(1 - 1e-5)

  expect_warning(
    make_vaccination_model(
      params = model_params,
      state = NULL,
      start_date = start_date, end_date = end_date,
      do_hazard = TRUE,
      do_approx_hazard = FALSE,
      do_make_state = TRUE,
      max_iters_eig_pow_meth = 100,
      tol_eig_pow_meth = 1e-6,
      params_timevar = params_timevar,
      do_variant = TRUE),
    regexp = "some time-varying parameters will not change"
  )

})

test_that("logit-scale time-variation works", {
  factory_fresh_macpan_options()
  options(MP_default_do_sim_constraint = TRUE)
  initial_gamma = 0.01
  gamma_multiplier = 0.6
  timevar = data.frame(
    Date = c("2020-04-01", "2020-05-01", "2020-06-01", "2020-07-01", "2020-08-01"),
    Symbol = "gamma",
    Value = gamma_multiplier,
    Type = "rel_prev_logit"
  )
  sir = (flexmodel(
      params = c(beta = 0.1, gamma = initial_gamma, N = 100),
      state = c(S = 99, I = 1, R = 0),
      start_date = "2020-03-11",
      end_date = "2020-12-01"
    )
    %>% update_piece_wise(timevar)
    %>% add_rate("S", "I", ~ (I) * (beta) * (1/N))
    %>% add_rate("I", "R", ~ (gamma))
  )
  expect_equal(
    simulation_history(sir)$I_to_R %>% unique,
    plogis(qlogis(initial_gamma) + (0:5) * gamma_multiplier)
  )
  synth_data = (sir
    %>% update_params(c(
      dist = 1
    ))
    %>% update_error_dist(
      S ~ negative_binomial('dist'),
      I ~ negative_binomial('dist'),
      R ~ negative_binomial('dist')
    )
    %>% simulation_history(obs_error = TRUE, include_initial_date = FALSE)
    %>% select(Date, S, I, R)
    %>% rename(date = Date)
    %>% pivot_longer(-date, names_to = 'var')
  )
  timevar_to_cal = (timevar
    %>% mutate(Value = NA)
  )
  sir_to_cal = (sir
    %>% update_piece_wise(timevar_to_cal)
    %>% update_observed(synth_data)
    %>% update_params(c(
      dist = 1
    ))
    %>% update_error_dist(
      S ~ negative_binomial('dist'),
      I ~ negative_binomial('dist'),
      R ~ negative_binomial('dist')
    )
    %>% update_opt_params(
      log_beta ~ log_flat(-2),
      log_gamma ~ log_flat(-2)
    )
    %>% update_opt_tv_params(tv_type = 'rel_prev_logit'
      , logit_gamma ~ logit_flat(0)
    )
  )
  # o = tmb_fun(sir_to_cal)
  # args = tmb_fun_args(sir_to_cal)
  # args$random = 'tv_mult'
  # o = do.call(MakeADFun, args)
  # opt = nlminb(o$par, o$fn, o$gr)
  # #o$env$parList()
  # opt$objective
  # o$env$random
  # o$env$retape()
  # o$fn()
  sir_cal = calibrate_flexmodel(sir_to_cal)
  sir_cal$opt_par
})

test_that("convolution paramters can be time-varying", {
  tv = data.frame(
    Date = "2016-09-02",
    Symbol = "c_prop",
    Value = 0.5,
    Type = 'abs'
  )
  sir = (make_hello_world_model()
    %>% update_params(
      c_prop = 1,
      c_delay_cv = 0.25,
      c_delay_mean = 11
    )
    %>% add_sim_report_expr("incidence", "^S_to_I$")
    %>% add_conv("^incidence$")
  )
  sir_with_tv = (sir
    %>% add_piece_wise(tv)
  )
  tv_conv = na.omit(simulation_history(sir_with_tv)$conv_incidence)
  conv = na.omit(simulation_history(sir)$conv_incidence)
  if (interactive()) {
    plot(conv, type = 'l')
    lines(tv_conv, lty = 2, col = 'red')
  }
  expect_equal(sort(unique(tv_conv / conv)), c(0.5, 1.0))
})
