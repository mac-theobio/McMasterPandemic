library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)
library(lubridate)

test_that('simple models calibrate the same regardless of engine', {
  set_spec_version("0.1.0", "../../inst/tmb/")

  params <- read_params("ICU1.csv")
  #state <- make_state(params = params)
  start_date = "2021-05-10"
  end_date = "2021-12-10"
  tv_dat <- data.frame(
    Date = c("2021-06-18", "2021-07-01", "2021-07-25"),
    Symbol = c("beta0", "beta0", "beta0"),
    Value = c(0.1, 0.5, 0.9),
    Type = c("rel_orig", "rel_orig", "rel_orig")
  )


  model <- make_base_model(params, start_date = start_date, end_date = end_date,
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
  head(report_data)
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
        #obj_fun = tmb_fun(model)
      )
    )
  )

  expect_equal(fitted_tmb$mle2@details, fitted_r$mle2@details)
})

test_that('simple models can calibrate time varying multipliers', {

  # TODO: to make this pass, we need to engage R-based make_state ... I think
  #       https://github.com/mac-theobio/McMasterPandemic/issues/124

  set_spec_version('0.1.1', '../../inst/tmb')
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
    use_flex = TRUE,
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

if(FALSE) {

  ## make synthetic data to fit to

  spec_version = "0.1.0"
  options(MP_flex_spec_version = spec_version)
  cpp <- file.path(test_files, spec_version, "macpan.cpp")
  dll <- tools::file_path_sans_ext(cpp)
  options(MP_flex_spec_dll = basename(dll))

  start_date <- "2020-02-01"
  end_date <- "2020-09-01"
  options(macpan_pfun_method = "grep")
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
    #time_params = c(0.8) ## initial guess for change in beta0 on the one and only break date (guess 80% of the base value)
  )
  ## (time-varying relative values: insert an NA wherever you want a parameter to be calibrated)
  params_timevar_calib <- (params_timevar
                           %>% mutate(Value = ifelse(Symbol == "beta0",
                                                     NA,  # tmb does not have this capability
                                                     Value))
  )
  ## fit!
  r_strt = Sys.time()
  fitted_mod <- calibrate(
    base_params = vax_params,
    data = synth_reports,
    opt_pars = opt_pars,
    debug = TRUE,
    time_args = list(
      params_timevar = params_timevar_calib #-- tmb can't yet handle optimizing changed parameters
    ),
    sim_args = list(
      ndt = 1,
      step_args = list(do_hazard = TRUE)
    ) ## there are both the defaults currently but i'm putting them here to emphasize that for the vaxified model's current implementation,
    ## we need ndt = 1 because the vax rate is specified as doses *per day*, implicitly assuming that the simulation algorithm takes daily steps
    ## and we need do_hazard = TRUE to avoid accidentally stepping into negative state space when vaccination occurs relatively quickly
  )
  r_nd = Sys.time()
  r_speed = as.numeric(r_nd - r_strt)

  test_model <- make_vaccination_model(
    params = vax_params, state = vax_state,
    params_timevar = params_timevar,
    start_date = start_date, end_date = end_date,
    step_args = list(do_hazard = TRUE)
  )

  tmb_strt = Sys.time()
  fitted_mod_tmb <- calibrate(
    base_params = vax_params,
    data = synth_reports,
    debug = TRUE,
    opt_pars = opt_pars,
    time_args = list(
      params_timevar = params_timevar#_calib -- tmb can't yet handle optimizing changed parameters
    ),
    sim_args = list(
      ndt = 1,
      step_args = list(do_hazard = TRUE),
      use_flex = TRUE,
      flexmodel = test_model,
      obj_fun = tmb_fun(test_model)
    ) ## there are both the defaults currently but i'm putting them here to emphasize that for the vaxified model's current implementation,
    ## we need ndt = 1 because the vax rate is specified as doses *per day*, implicitly assuming that the simulation algorithm takes daily steps
    ## and we need do_hazard = TRUE to avoid accidentally stepping into negative state space when vaccination occurs relatively quickly
  )
  tmb_nd = Sys.time()
  tmb_speed = as.numeric(tmb_nd - tmb_strt)
  r_speed / tmb_speed

  fitted_mod
  fitted_mod_tmb

}
