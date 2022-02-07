library(McMasterPandemic)
library(dplyr)
library(lubridate)
library(testthat)

test_that("macpan ontario forecasts work the same regardless of engine", {
  rerun_r_engine_calibrations = FALSE
  rerun_r_engine_forecasts = FALSE
  run_simulation_comparison = FALSE
  reset_spec_version()
  options(macpan_pfun_method = "grep")
  options(MP_rexp_steps_default = 150)
  options(MP_warn_bad_breaks = FALSE)
  r_tmb_comparable()
  load("../../inst/testdata/ontario_flex_test_better.Rdata")

  model = make_vaccination_model(
    params = model_params,
    state = NULL,
    start_date = ymd("20200915") - days(start_date_offset),
    end_date = ymd("20211214"),
    params_timevar = params_timevar,
    do_hazard = TRUE,
    do_make_state = TRUE,
    do_variant = TRUE
  )

  tmb_init = initial_state_vector(model)
  r_init = make_state(params = model_params)
  all.equal(c(r_init), c(tmb_init))

  params_timevar_sim = (params_timevar
                        %>% within({Value[is.na(Value)] = opt_pars$time_params})
  )

  if(run_simulation_comparison) {
    r_sim = run_sim(
      params = model_params,
      start_date = model$start_date,
      end_date = model$end_date,
      state = NULL,
      params_timevar = params_timevar_sim,
      step_args = list(do_hazard = TRUE)
    )
    tmb_sim = run_sim(
      params = model_params,
      start_date = model$start_date,
      end_date = model$end_date,
      state = NULL,
      params_timevar = params_timevar_sim,
      step_args = list(do_hazard = TRUE),
      flexmodel = model
    )
    compare_sims(r_sim, tmb_sim, compare_attr = FALSE)
  }

  # suppressing warning for now -- nelder-mead complaining
  # because it is a 1d problem (TODO: ask Mike if he is getting
  # this in production)
  fitted_model_tmb <- suppressWarnings(calibrate(
    base_params  = model_params
    , data       = fitdat
    , opt_pars   = opt_pars
    , sim_args   = list(ndt = 1,
                        step_args = list(do_hazard = TRUE),
                        flexmodel = model)
    , time_args = list(params_timevar = params_timevar)
    , mle2_control = list(maxit = 1e3,
                          reltol = calibration_params[["reltol"]],
                          ## here beta and gamma are nelder-mead parameters,
                          ## not macpan model parameters!!!
                          beta = 0.5,
                          gamma = 2
    )
    , start_date_offset = start_date_offset
    , use_DEoptim = FALSE
    , debug = FALSE
    , debug_plot = FALSE
  ))

  if (rerun_r_engine_calibrations) {
    # takes a really long time so best to avoid it
    fitted_model_r <- calibrate(
      base_params  = model_params
      , data       = fitdat
      , opt_pars   = opt_pars
      , sim_args   = list(ndt = 1,
                          step_args = list(do_hazard = TRUE))
      , time_args = list(params_timevar = params_timevar)
      , mle2_control = list(maxit = 1e3,
                            reltol = calibration_params[["reltol"]],
                            ## here beta and gamma are nelder-mead parameters,
                            ## not macpan model parameters!!!
                            beta = 0.5,
                            gamma = 2
      )
      , start_date_offset = start_date_offset
      , use_DEoptim = FALSE
      , debug = FALSE
      , debug_plot = FALSE
    )
    save(fitted_model_r, file = "../../inst/testdata/ontario_flex_test_better_fit.rda")
  } else {
    load("../../inst/testdata/ontario_flex_test_better_fit.rda")
  }
  # [[1]]
  # Time difference of 4.138171 secs
  #
  # [[2]]
  # Time difference of 2174.062 secs
  #
  # [[3]]
  # [1] 525.3678
  # save(fitted_model_r, file = "../../inst/testdata/ontario_flex_test_better_fit.Rdata")

  expect_equal(fitted_model_r$mle2@coef, fitted_model_tmb$mle2@coef)
  expect_equal(fitted_model_r$mle2@min, fitted_model_tmb$mle2@min)
  expect_equal(fitted_model_r$mle2@vcov, fitted_model_tmb$mle2@vcov)

  forecast_tmb = forecast_ensemble(fitted_model_tmb, nsim = 2, seed = 1)
  if (rerun_r_engine_forecasts) {
    forecast_r = forecast_ensemble(fitted_model_r, nsim = 2, seed = 1)
    save(forecast_r, file = "../../inst/testdata/ontario_flex_test_better_forecasts.rda")
  } else {
    load("../../inst/testdata/ontario_flex_test_better_forecasts.rda")
  }
  expect_equal(forecast_tmb, forecast_r)
})
