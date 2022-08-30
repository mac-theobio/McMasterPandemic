library(McMasterPandemic)
library(dplyr)
library(lubridate)
library(testthat)

testLevel <- if (nzchar(s <- Sys.getenv("MACPAN_TEST_LEVEL"))) as.numeric(s) else 1
skip_slow_tests = isTRUE(testLevel == 1)

test_that("macpan ontario forecasts work the same regardless of engine", {
  skip_if(skip_slow_tests)

  rerun_r_engine_calibrations = FALSE
  rerun_r_engine_forecasts = FALSE
  run_simulation_comparison = FALSE
  #reset_spec_version()

  set_spec_version("0.2.0", system.file('tmb', package = 'McMasterPandemic'))
  options(macpan_pfun_method = "grep")
  options(MP_rexp_steps_default = 150)
  options(MP_warn_bad_breaks = FALSE)
  r_tmb_comparable()
  load(system.file('testdata', 'ontario_flex_test_better.Rdata', package = 'McMasterPandemic'))

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
    compare_sims(r_sim, tmb_sim, compare_attr = FALSE, na_is_zero = TRUE)
  }

  cntrl = list(maxit = 1e3,
       reltol = calibration_params[["reltol"]],
       ## here beta and gamma are nelder-mead parameters,
       ## not macpan model parameters!!!
       beta = 0.5,
       gamma = 2
  )

  # suppressing warning for now -- nelder-mead complaining
  # because it is a 1d problem (TODO: ask Mike if he is getting
  # this in production)
  model$do_sim_constraint = TRUE
  opt_pars_alt = opt_pars
  opt_pars_alt$time_params = 1
  fitted_model_tmb <- suppressWarnings(calibrate(
    base_params  = model_params
    , data       = fitdat
    , opt_pars   = opt_pars_alt
    , sim_args   = list(ndt = 1,
                        step_args = list(do_hazard = TRUE),
                        flexmodel = model)
    , time_args = list(params_timevar = params_timevar)
    #, mle2_control = cntrl
    , start_date_offset = start_date_offset
    , use_DEoptim = FALSE
    , debug = FALSE
    , debug_plot = FALSE
  ))


  model_to_fit = (model
    #%>% update_piece_wise(params_timevar_sim)
    %>% update_observed(fitdat)
    #%>% update_opt_params(
    #  log_nb_disp_report ~ log_normal(0, 5)
    #)
    %>% update_opt_tv_params("rel_orig"
      , beta0 ~ flat(1)
    )
    %>% update_tmb_indices
  )

  expect_true(
    isTRUE(
      compare_grads(
        model_to_fit,
        tmb_pars = tmb_params_trans(model_to_fit)
      )
    )
  )
  fitted_model_tmb_condense = calibrate_flexmodel(model_to_fit)

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
    save(
      fitted_model_r,
      file = system.file('testdata', 'ontario_flex_test_better_fit.Rdata', package = 'McMasterPandemic')
     )
      #file = file.path(pkg_home, "inst/testdata/ontario_flex_test_better_fit.rda"))
  } else {
    load(system.file('testdata', 'ontario_flex_test_better_fit.Rdata', package = 'McMasterPandemic'))
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

  expect_equal(
    unname(fitted_model_tmb_condense$opt_par),
    unname(fitted_model_tmb$mle2@coef),
    tolerance = 1e-3
  )
  expect_equal(
    fitted_model_r$mle2@coef,
    fitted_model_tmb$mle2@coef,
    tolerance = 1e-4
  )
  expect_equal(fitted_model_r$mle2@min, fitted_model_tmb$mle2@min)
  expect_equal(fitted_model_r$mle2@vcov, fitted_model_tmb$mle2@vcov)

  forecast_tmb = forecast_ensemble(fitted_model_tmb, nsim = 2, seed = 1)
  if (rerun_r_engine_forecasts) {
    forecast_r = forecast_ensemble(fitted_model_r, nsim = 2, seed = 1)
    save(forecast_r, file = system.file('testdata', 'ontario_flex_test_better_forecasts.rda', package = 'McMasterPandemic'))
  } else {
    load(system.file('testdata', 'ontario_flex_test_better_forecasts.rda', package = 'McMasterPandemic'))
  }

  Rprof("~/testing.out")
  forecast_tmb_condense = simulate_ensemble(fitted_model_tmb_condense, n = 200)
  Rprof(NULL)

  expect_equal(forecast_tmb, forecast_r, tolerance = 1e-6)

  if (FALSE) {
    # these look similar, but the variation is very small in tmb??
    # although it does make sense that there is not any variation until
    # the breakpoint of the only (time-varying) parameter
    library(ggplot2)
    (forecast_tmb_condense
      %>% filter(name == "report")
      %>% ggplot
      + geom_ribbon(aes(Date, ymin = lwr, ymax = upr), alpha = 0.2)
      + geom_line(aes(Date, value))
    )
    (forecast_tmb
      %>% filter(var == "report")
      %>% ggplot
      + geom_ribbon(aes(date, ymin = lwr, ymax = upr), alpha = 0.2)
      + geom_line(aes(date, value))
    )
  }
})
