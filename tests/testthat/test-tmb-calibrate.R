library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)
library(lubridate)

test_that('v0.1.0 simple models calibrate the same regardless of engine', {
  set_spec_version("0.1.0", "../../inst/tmb/")

  params <- read_params("ICU1.csv")
  start_date = "2021-05-10"
  end_date = "2021-12-10"
  tv_dat <- data.frame(
    Date = c("2021-06-18", "2021-07-01", "2021-07-25"),
    Symbol = c("beta0", "beta0", "beta0"),
    Value = c(0.1, 0.5, 0.9),
    Type = c("rel_orig", "rel_orig", "rel_orig")
  )

  model <- make_base_model(params,
                           state = make_state(params = params),
                           start_date = start_date, end_date = end_date,
                           params_timevar = tv_dat,
                           do_hazard = TRUE)

  obj <- tmb_fun(model)

  tmb_sim <- run_sim(
    params = params, state = model$state,
    start_date = start_date,
    end_date = end_date,
    params_timevar = tv_dat,
    step_args = list(do_hazard = TRUE),
    condense = TRUE,
    use_flex = TRUE,
    flexmodel = model,
    obj_fun = obj
  )

  r_sim <- run_sim(
    params = params, state = model$state,
    start_date = start_date,
    end_date = end_date,
    params_timevar = tv_dat,
    step_args = list(do_hazard = TRUE),
    condense = TRUE,
    use_flex = FALSE
  )

  compare_sims(r_sim, tmb_sim)

  report_data <- (r_sim
                  %>% mutate(value = round(report), var = "report")
                  %>% select(date, value, var)
                  %>% na.omit()
  )
  if(FALSE) plot(report_data$value, type = "l")

  params[["beta0"]] = 0.5

  time_wrap(
    fitted_r <- calibrate(
      data = report_data,
      time_args = list(params_timevar = tv_dat),
      base_params = params,
      opt_pars = list(params = c(beta0 = params[["beta0"]])),
      sim_args = list(
        step_args = list(do_hazard = TRUE),
        use_flex = FALSE
      )
    ),
    fitted_tmb <- calibrate(
      data = report_data,
      time_args = list(params_timevar = tv_dat),
      base_params = params,
      opt_pars = list(params = c(beta0 = params[["beta0"]])),
      sim_args = list(
        step_args = list(do_hazard = TRUE),
        use_flex = TRUE,
        flexmodel = model
      )
    )
  )

  expect_equal(fitted_tmb$mle2@details, fitted_r$mle2@details)
})

test_that('v0.1.1 simple models can calibrate time varying multipliers', {

  reset_spec_version()
  tmb_mode()
  options(MP_use_state_rounding = FALSE)

  params <- ("ICU1.csv"
             %>% read_params
             %>% expand_params_S0(1 - 1e-5)
  )
  start_date = "2021-05-10"
  end_date = "2021-12-10"
  tv_dat <- data.frame(
    Date = c("2021-06-18", "2021-07-01", "2021-07-25"),
    Symbol = c("beta0", "beta0", "beta0"),
    Value = c(0.1, NA, 0.9),
    Type = c("rel_orig", "rel_orig", "rel_orig")
  )
  model <- make_base_model(params,
                           state = make_state(params = params),
                           start_date = start_date, end_date = end_date,
                           params_timevar = tv_dat,
                           do_hazard = TRUE,
                           do_make_state = TRUE,
                           max_iters_eig_pow_meth = 100,
                           tol_eig_pow_meth = 1e-03
  )

  tv_dat_filled = tv_dat
  tv_dat_filled$Value = model$timevar$piece_wise$schedule$init_tv_mult

  tmb_sim <- run_sim(
    params = params, state = model$state,
    start_date = start_date,
    end_date = end_date,
    params_timevar = tv_dat_filled,
    step_args = list(do_hazard = TRUE),
    condense = TRUE,
    flexmodel = model
  )

  report_data <- (tmb_sim
                  %>% mutate(value = round(report), var = "report")
                  %>% select(date, value, var)
                  %>% na.omit()
  )

  params[["beta0"]] = 0.5

  time_wrap(
    fitted_r <- calibrate(
      data = report_data,
      time_args = list(params_timevar = tv_dat),
      base_params = params,
      debug = TRUE,
      opt_pars = list(params = c(beta0 = params[["beta0"]]),
                      time_params = c(0.5)),
      sim_args = list(
        step_args = list(do_hazard = TRUE),
        use_flex = FALSE
      )
    ),
    fitted_tmb <- calibrate(
      data = report_data,
      time_args = list(params_timevar = tv_dat),
      base_params = params,
      debug = TRUE,
      opt_pars = list(params = c(beta0 = params[["beta0"]]),
                      time_params = c(0.5)),
      sim_args = list(
        step_args = list(do_hazard = TRUE),
        use_flex = TRUE,
        flexmodel = model,
        obj_fun = tmb_fun(model)
      )
    )
  )

  expect_equal(fitted_r$mle2@coef, fitted_tmb$mle2@coef)
  expect_equal(fitted_r$mle2@min, fitted_tmb$mle2@min)
  expect_equal(fitted_r$mle2@vcov, fitted_tmb$mle2@vcov)
})

test_that("v0.1.1 vaccination models calibrate the same regardless of engine", {

  reset_spec_version()
  tmb_mode()
  options(macpan_pfun_method = "grep")
  options(MP_rexp_steps_default = 150)

  start_date <- "2020-02-01"
  end_date <- "2020-09-01"
  params <- read_params("ICU1.csv")
  state <- make_state(params = params)
  vax_params <- expand_params_vax(
    params = params,
    model_type = "twodose"
  )
  vax_state <- expand_state_vax(
    x = state,
    model_type = "twodose",
    unif = FALSE
  )

  ## set up time-varying parameters
  params_timevar <- data.frame(
    Date = c(as.Date(start_date) + 30,
             as.Date(start_date) + 60),
    Symbol = c("vax_prop_first_dose", "beta0"),
    Value = rep(0.5, 2),
    Type = rep("rel_orig", 2)
  )

  ## generate reports from sim
  synth_reports <- (run_sim(
      params = vax_params,
      start_date = start_date,
      end_date = end_date,
      # do the same thing with the switch to second doses as above, but now also cut the original transmission rate to 50% of its value 60 days after the simulation start date
      params_timevar = params_timevar
    )
    ## reshape into the correct format for input data passed to calibrate()
    %>% mutate(value=round(report), var="report")
    %>% select(date, value, var)
    %>% na.omit()
  )

  ## set up optimization parameters
  ## (base parameter values)
  opt_pars <- list(
    params = c(beta0 = 0.6) ## set initial guess for beta0
  )

  state_nms =
    c("S_unvax", "E_unvax", "Ia_unvax", "Ip_unvax", "Im_unvax", "Is_unvax",
    "H_unvax", "H2_unvax", "ICUs_unvax", "ICUd_unvax", "D_unvax",
    "R_unvax", "X_unvax", "V_unvax", "S_vaxdose1", "E_vaxdose1",
    "Ia_vaxdose1", "Ip_vaxdose1", "Im_vaxdose1", "Is_vaxdose1", "H_vaxdose1",
    "H2_vaxdose1", "ICUs_vaxdose1", "ICUd_vaxdose1", "D_vaxdose1",
    "R_vaxdose1", "X_vaxdose1", "V_vaxdose1", "S_vaxprotect1", "E_vaxprotect1",
    "Ia_vaxprotect1", "Ip_vaxprotect1", "Im_vaxprotect1", "Is_vaxprotect1",
    "H_vaxprotect1", "H2_vaxprotect1", "ICUs_vaxprotect1", "ICUd_vaxprotect1",
    "D_vaxprotect1", "R_vaxprotect1", "X_vaxprotect1", "V_vaxprotect1",
    "S_vaxdose2", "E_vaxdose2", "Ia_vaxdose2", "Ip_vaxdose2", "Im_vaxdose2",
    "Is_vaxdose2", "H_vaxdose2", "H2_vaxdose2", "ICUs_vaxdose2",
    "ICUd_vaxdose2", "D_vaxdose2", "R_vaxdose2", "X_vaxdose2", "V_vaxdose2",
    "S_vaxprotect2", "E_vaxprotect2", "Ia_vaxprotect2", "Ip_vaxprotect2",
    "Im_vaxprotect2", "Is_vaxprotect2", "H_vaxprotect2", "H2_vaxprotect2",
    "ICUs_vaxprotect2", "ICUd_vaxprotect2", "D_vaxprotect2", "R_vaxprotect2",
    "X_vaxprotect2", "V_vaxprotect2")

  test_model <- make_vaccination_model(
    params = vax_params,
    state = make_state(params = vax_params),
    params_timevar = params_timevar,
    start_date = start_date, end_date = end_date,
    do_hazard = TRUE,
    do_hazard_lin = FALSE,
    do_approx_hazard = FALSE,
    do_approx_hazard_lin = FALSE,
    do_make_state = TRUE,
    max_iters_eig_pow_meth = 1000,
    tol_eig_pow_meth = 1e-12
  )

  simulate_timings = time_wrap(
    sims_tmb <- run_sim(
      params = test_model$params,
      params_timevar = params_timevar,
      start_date = start_date, end_date = end_date,
      step_args = list(do_hazard = TRUE),
      condense = FALSE,
      flexmodel = test_model
    ),
    sims_r <- run_sim(
      params = vax_params,
      params_timevar = params_timevar,
      start_date = start_date, end_date = end_date,
      step_args = list(do_hazard = TRUE),
      condense = FALSE
    )
  )
  compare_sims(sims_r, sims_tmb, compare_attr = FALSE)

  calibrate_timings = time_wrap(
    fitted_mod_tmb <- calibrate(
      base_params = test_model$params,
      data = synth_reports,
      opt_pars = opt_pars,
      debug = TRUE,
      time_args = list(
        params_timevar = params_timevar
      ),
      sim_args = list(
        ndt = 1,
        step_args = list(do_hazard = TRUE),
        use_flex = TRUE,
        flexmodel = test_model
      )
    ),
    fitted_mod <- calibrate(
      base_params = vax_params,
      data = synth_reports,
      opt_pars = opt_pars,
      debug = TRUE,
      time_args = list(
        params_timevar = params_timevar
      ),
      sim_args = list(
        ndt = 1,
        step_args = list(do_hazard = TRUE)
      )
    )
  )

  calibrate(
    base_params = test_model$params,
    data = synth_reports,
    opt_pars = opt_pars,
    debug = TRUE,
    time_args = list(
      params_timevar = params_timevar
    ),
    sim_args = list(
      ndt = 1,
      step_args = list(do_hazard = TRUE),
      use_flex = TRUE,
      flexmodel = test_model,
      obj_fun = tmb_fun(test_model)
    )
  )

  expect_equal(fitted_mod$mle2@coef, fitted_mod_tmb$mle2@coef)
  expect_equal(fitted_mod$mle2@min, fitted_mod_tmb$mle2@min)
  expect_equal(fitted_mod$mle2@vcov, fitted_mod_tmb$mle2@vcov)
})

test_that("v0.1.1 vaccination models can calibrate time varying multipliers", {

  reset_spec_version()
  tmb_mode()
  options(macpan_pfun_method = "grep")
  options(MP_rexp_steps_default = 150)

  start_date <- "2020-02-01"
  end_date <- "2020-09-01"

  params <- read_params("ICU1.csv")
  state <- make_state(params = params)
  vax_params <- expand_params_vax(
    params = params,
    model_type = "twodose"
  )
  vax_state <- expand_state_vax(
    x = state,
    model_type = "twodose",
    unif = FALSE
  )

  ## set up time-varying parameters
  params_timevar <- data.frame(
    Date = c(as.Date(start_date) + 30,
             as.Date(start_date) + 60),
    Symbol = c("vax_prop_first_dose", "beta0"),
    Value = rep(0.5, 2),
    Type = rep("rel_orig", 2)
  )

  ## generate reports from sim
  synth_reports <- (run_sim(
    params = vax_params,
    start_date = start_date,
    end_date = end_date,
    # do the same thing with the switch to second doses as above, but now also cut the original transmission rate to 50% of its value 60 days after the simulation start date
    params_timevar = params_timevar
  )
    ## reshape into the correct format for input data passed to calibrate()
    %>% mutate(value=round(report), var="report")
    %>% select(date, value, var)
    %>% na.omit()
  )

  ## set up optimization parameters
  ## (base parameter values)
  opt_pars <- list(
    params = c(beta0 = 0.6), ## set initial guess for beta0
    time_params = c(0.8) ## initial guess for change in beta0 on the one and only break date (guess 80% of the base value)
  )
  ## (time-varying relative values: insert an NA wherever you want a parameter to be calibrated)
  params_timevar_calib <- (params_timevar
                           %>% mutate(Value = ifelse(Symbol == "beta0",
                                                     NA,
                                                     Value))
  )

  fitted_mod <- calibrate(
    base_params = vax_params,
    data = synth_reports,
    opt_pars = opt_pars,
    debug = TRUE,
    time_args = list(
      params_timevar = params_timevar_calib
    ),
    sim_args = list(
      ndt = 1,
      step_args = list(do_hazard = TRUE)
    )
  )

  test_model <- make_vaccination_model(
    params = vax_params,
    state = vax_state,
    params_timevar = params_timevar,
    start_date = start_date, end_date = end_date,
    step_args = list(do_hazard = TRUE)
  )

  fitted_mod_tmb <- calibrate(
    base_params = test_model$params,
    data = synth_reports,
    opt_pars = opt_pars,
    debug = TRUE,
    time_args = list(
      params_timevar = params_timevar_calib
    ),
    sim_args = list(
      ndt = 1,
      step_args = list(do_hazard = TRUE),
      use_flex = TRUE,
      flexmodel = test_model
    )
  )

  expect_equal(fitted_mod$mle2@coef, fitted_mod_tmb$mle2@coef)
  expect_equal(fitted_mod$mle2@min, fitted_mod_tmb$mle2@min)
  expect_equal(fitted_mod$mle2@vcov, fitted_mod_tmb$mle2@vcov)
})
