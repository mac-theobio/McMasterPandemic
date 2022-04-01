library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)
library(lubridate)

test_that('simple models calibrate the same regardless of engine', {
  reset_spec_version()
  tmb_mode()

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

  tmb_sim <- run_sim(
    params = model$params,
    state = model$state,
    start_date = start_date,
    end_date = end_date,
    params_timevar = tv_dat,
    step_args = list(do_hazard = TRUE),
    condense = TRUE,
    flexmodel = model
  )

  r_sim <- run_sim(
    params = model$params,
    state = model$state,
    start_date = start_date,
    end_date = end_date,
    params_timevar = tv_dat,
    step_args = list(do_hazard = TRUE),
    condense = TRUE
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
                           max_iters_eig_pow_meth = 300,
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
      debug = FALSE,
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
      debug = FALSE,
      opt_pars = list(params = c(beta0 = params[["beta0"]]),
                      time_params = c(0.5)),
      sim_args = list(
        flexmodel = model
      )
    )
  )

  expect_equal(fitted_r$mle2@coef, fitted_tmb$mle2@coef)
  expect_equal(fitted_r$mle2@min, fitted_tmb$mle2@min)
  expect_equal(fitted_r$mle2@vcov, fitted_tmb$mle2@vcov)
})

if(FALSE) {
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
    tol_eig_pow_meth = 1e-4
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
      debug = FALSE,
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
      debug = FALSE,
      time_args = list(
        params_timevar = params_timevar
      ),
      sim_args = list(
        ndt = 1,
        step_args = list(do_hazard = TRUE)
      )
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
    debug = FALSE,
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
    debug = FALSE,
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
}

test_that("macpan ontario calibration example works the same regardless of engine", {
  rerun_r_engine_calibrations = FALSE
  run_simulation_comparison = FALSE
  reset_spec_version()
  #r_tmb_comparable()
  options(macpan_pfun_method = "grep")
  options(MP_rexp_steps_default = 150)
  load("../../inst/testdata/ontario_flex_test_better.Rdata")

  model = make_vaccination_model(
    params = model_params,
    state = NULL,
    start_date = ymd("20200915") - days(start_date_offset),
    end_date = ymd("20211215"),
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

  # suppressing warning for now -- nelder-mead compaining
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
    save(fitted_model_r, "../../inst/testdata/ontario_flex_test_better_fit.rda")
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
})

test_that("tmb engine calibrates correctly to multiple data streams", {

  library(McMasterPandemic)
  library(lubridate)
  library(tidyr)
  library(dplyr)

  refit_no_flex_model = FALSE # slow if TRUE

  reset_spec_version()
  tmb_mode()

  params1 <- read_params("PHAC.csv")

  # update default parameters for the Yukon case
  #   beta0 -- baseline transmission rate
  #   N -- population of Yukon
  #   mu -- proportion of cases that are mild
  #   phi1 -- proportion of hospitalization cases that are not in the ICU
  params1[c("beta0", "N", "mu", "phi1")] <- c(0.1, 42507, 0.95, 0.98)

  # initial state of the simulation
  state1 <- make_state(params=params1)

  # start and end dates
  sdate <- as.Date("2021-11-01")
  edate <- as.Date("2022-01-19")
  initial_date = as.Date("2021-08-03")
  start_date_offset = as.integer(sdate - initial_date)

  # read and process data
  covid_data <- ("../../sandbox/yukon/report_data_yukon_h_and_i.csv"
                 %>% read.csv
                 %>% mutate(date = as.Date(date))
                 %>% filter(date >= ymd(20210803))
                 %>% filter(between(as.Date(date), sdate, edate))
                 %>% select(date, report, death, hosp, ICU)
                 %>% pivot_longer(names_to = "var", -date)
                 %>% mutate(value=round(value))
  )

  # establish schedule of time variation of parameters
  params_timevar = data.frame(
    Date = ymd(
      20211115, # nov 15
      20211215, # dec 15
      20220101  # jan 01
    ),
    Symbol = "beta0",
    Value = NA,
    Type = "abs"
  )

  # declare parameters to be fitted
  #   - two types of parameters: params and time_params
  #   - params is a named vector with names referring to
  #     the parameter to be optimized
  #   - optionally change the scale on which optimization occur
  #     by prefixing these names with log_ or logit_
  #   - values of this vector are initial values fed to the
  #     optimizer
  #   - time_params is a vector of values in the order of
  #     of the params_timevar schedule above
  opt_pars <- list(
    params = c(

      # baseline transmission rate
      # (fit on log scale to avoid negative rates)
      log_beta0 = log(params1[["beta0"]]),

      # proportion of cases that are mild
      # (fit on logit scale to keep proprotions between zero and one)
      logit_mu = qlogis(params1[["mu"]]),

      # proportion of hospitalizations that are not in the ICU
      # (fit on logit scale to keep proprotions between zero and one)
      logit_phi1 = qlogis(params1[["phi1"]])
    ),

    # time varying transmissions rates
    # (see params_time_var above for schedule of changes in these rates)
    # (fit on log scale to avoid negative rates)
    log_time_params = rep(
      log(params1[["beta0"]]),
      nrow(params_timevar)
    )
  )

  mm = make_base_model(
    params = params1,
    start_date = sdate - start_date_offset,
    end_date = edate,
    params_timevar = params_timevar,
    do_make_state = TRUE,
    do_hazard = TRUE
  )
  sim_args_flex = list(flexmodel = mm)
  sim_args = list()

  # fit the models
  if (refit_no_flex_model) {
    fit_no_flex <- calibrate(
      data = covid_data,
      time_args = list(params_timevar = params_timevar),
      start_date_offset = start_date_offset,
      base_params = params1,
      opt_pars = opt_pars,
      debug = FALSE,
      sim_args = sim_args
    )
    save(fit_no_flex, file = '../../inst/testdata/yukon_no_flex.rda')
  } else {
    load('../../inst/testdata/yukon_no_flex.rda')
  }
  fit_flex <- calibrate(
    data = covid_data,
    time_args = list(params_timevar = params_timevar),
    start_date_offset = start_date_offset,
    base_params = params1,
    opt_pars = opt_pars,
    debug = FALSE,
    sim_args = sim_args_flex
  )

  expect_equal(fit_no_flex$mle2@coef, fit_flex$mle2@coef)
  expect_equal(fit_no_flex$mle2@min, fit_flex$mle2@min)
  expect_equal(fit_no_flex$mle2@vcov, fit_flex$mle2@vcov)
})
