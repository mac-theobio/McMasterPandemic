Sys.setenv(R_TESTS="")

library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)
library(lubridate)

test_that('time variation works for several parameters and levels of continuity', {
  reset_spec_version()
  tmb_mode()
  params <- read_params("ICU1.csv")
  test_fun = function(tv_dat) {
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

    compare_sims(r_sim, tmb_sim)
  }

  data.frame(
    Date = as.Date(c("2021-09-15", "2021-09-20", "2021-10-05", "2021-09-20", "2021-10-03")),
    Symbol = c("beta0", "beta0", "beta0", "Ca", "Ca"),
    Value = c(0.5, 0.1, 0.05, 0.1, 2.3),
    Type = c("rel_orig", "rel_orig", "rel_orig", "rel_orig", "rel_orig")
  ) %>% test_fun

  data.frame(
    Date = ymd(20210910) + days(0:1),
    Symbol = "beta0",
    Value = seq(from = 0.1, to = 2, length.out = 2),
    Type = "rel_orig"
  ) %>% test_fun

  data.frame(
   Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
   Symbol = c("beta0", "beta0", "beta0"),
   Value = c(0.5, 0.1, 0.05),
   Type = c("rel_orig", "rel_orig", "rel_orig")
  ) %>% test_fun
})

test_that("time variation on the first two steps matches in r and tmb engines", {
  reset_spec_version()
  tmb_mode()
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
    compare_sims(r_sim, tmb_sim)
  }

  sapply(c('rel_orig', 'rel_prev', 'abs'), test_fun)
})

test_that('time variation works for vax models', {
  reset_spec_version()
  tmb_mode()
  options(macpan_pfun_method = "grep")
  options(MP_rexp_steps_default = 150)

  load("../../inst/testdata/ontario_flex_test.rda")

  start_date = min(params_timevar$Date)
  end_date = max(params_timevar$Date) + days(1)

  # (params_timevar
  #   %>% ggplot2::ggplot(ggplot2::aes(x = Date, y = Value))
  #    + ggplot2::geom_line()
  #    + ggplot2::facet_wrap(~Symbol, scales = "free_y")
  # )

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
      params = model_params,
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
  tmb_mode()
  params <- read_params("ICU1.csv")
  test_fun = function(tv_dat) {
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

    compare_sims(r_sim, tmb_sim)
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