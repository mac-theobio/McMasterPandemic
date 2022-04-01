library(McMasterPandemic)
library(dplyr)
library(testthat)
library(lubridate)

test_that('foi can be expressed as model structure', {
  reset_spec_version()
  tmb_mode()

  params <- read_params("ICU1.csv")

  mm = (flexmodel(
    params = params,
    state = make_state(params = params),
    do_make_state = FALSE,
    do_hazard = TRUE,
    start_date = "2021-09-10", end_date = "2021-10-10",
  )
    # force of infection (FOI)
    %>% add_rate("S", "E",
      sum(
        vec('Ia', 'Ip', 'Im', 'Is')
      * vec('Ca', 'Cp', '(Cm) * (1 - iso_m)', '(Cs) * (1 - iso_s)')
      * vec('beta0')
      * vec('1/N')
      )
    )

    # flows out of the exposed class
    %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
    %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))

    # flows out of the pre-symptomatic compartment
    %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
    %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))

    # flows out of the severe symptomatic compartment
    %>% add_rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
    %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
    %>% add_rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
    %>% add_rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))

    # flows out of the ICU
    %>% add_rate("ICUs", "H2", ~ (psi1))
    %>% add_rate("ICUd", "D", ~ (psi2))

    # recovery
    %>% add_rate("Im", "R", ~ (gamma_m))
    %>% add_rate("Ia", "R", ~ (gamma_a))
    %>% add_rate("H2", "R", ~ (psi3))
    %>% add_rate("H", "R", ~ (rho))


    # accumulators
    %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
    %>% add_outflow('.+', all_except("X"))

    %>% update_tmb_indices()
  )

  tmb_sim <- run_sim(
    params = params, state = mm$state,
    start_date = "2021-09-10",
    end_date = "2021-10-10",
    condense = FALSE,
    step_args = list(do_hazard = TRUE),
    flexmodel = mm
  )

  r_sim <- run_sim(
    params = params, state = mm$state,
    start_date = "2021-09-10",
    end_date = "2021-10-10",
    condense = FALSE,
    step_args = list(do_hazard = TRUE),
  )

  compare_sims(r_sim, tmb_sim, compare_attr = FALSE)
})

test_that('simple models still work when structure is allowed', {

  reset_spec_version()
  tmb_mode()

  params <- read_params("ICU1.csv")

  tv_dat <- data.frame(
    Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
    Symbol = c("beta0", "beta0", "beta0"),
    Value = c(0.5, 0.1, 0.05),
    Type = c("rel_prev", "rel_orig", "rel_prev")
  )

  mm = make_base_model(
    params,
    state = make_state(params = params),
    start_date = "2021-09-10",
    end_date = "2021-10-10",
    params_timevar = tv_dat,
    do_hazard = TRUE, do_make_state = FALSE
  )

  # change beta0, which is time-varying, so that
  # we can check that the parameter updates are
  # happening correctly in the C++ side
  test_pars = params
  test_pars[1] = 3

  tmb_sim <- run_sim(
    params = test_pars, state = mm$state,
    start_date = "2021-09-10",
    end_date = "2021-10-10",
    params_timevar = tv_dat,
    condense = FALSE,
    flexmodel = mm
  )

  r_sim <- run_sim(
    params = test_pars, state = mm$state,
    start_date = "2021-09-10",
    end_date = "2021-10-10",
    params_timevar = tv_dat,
    condense = FALSE,
    step_args = list(do_hazard = TRUE)
  )
  compare_sims(r_sim, tmb_sim)
})

test_that("simple vaccination model in TMB matches and is faster than existing R model", {
  options(macpan_pfun_method = "grep")
  reset_spec_version()
  tmb_mode()

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
  vax_params = expand_params_S0(vax_params, 1-1e-5)
  test_model <- make_vaccination_model(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    do_make_state = FALSE,
    do_hazard = TRUE
  )

  tmb_strt = Sys.time()
  tmb_sim <- run_sim(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    condense = TRUE,
    step_args = list(do_hazard = TRUE),
    flexmodel = test_model,
    use_flex = TRUE
  )
  tmb_nd = Sys.time()
  r_strt = Sys.time()
  r_sim = run_sim(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    condense = TRUE,
    step_args = list(do_hazard = TRUE)
  )
  r_nd = Sys.time()
  tmb_speed = as.numeric(tmb_nd - tmb_strt)
  r_speed = as.numeric(r_nd - r_strt)
  print(paste0('tmb speed-up: ', round(r_speed/tmb_speed), 'x'))
  expect_lt(tmb_speed, r_speed)
  compare_sims(r_sim, tmb_sim)
})

test_that("vax_prop_first_dose can be != 1", {
  reset_spec_version()
  tmb_mode()
  options(macpan_pfun_method = "grep")

  params <- read_params("ICU1.csv")
  state <- make_state(params = params)
  vax_params <- expand_params_vax(
    params = params,
    model_type = "twodose"
  )
  vax_params[["vax_prop_first_dose"]] = 0.5
  vax_state <- expand_state_vax(
    x = state,
    model_type = "twodose",
    unif = FALSE
  )
  vax_state[["E_vaxprotect1"]] = 3
  vax_state[["S_unvax"]] = vax_state[["S_unvax"]] - 3

  params_timevar = data.frame(
    Date = c("2021-09-20"),
    Symbol = c("beta0"),
    Value = c(0.5),
    Type = c("rel_orig")
  )

  test_model <- make_vaccination_model(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    params_timevar = params_timevar,
    do_hazard = TRUE,
    do_make_state = FALSE
  )

  tmb_strt = Sys.time()
  tmb_sim <- run_sim(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    condense = FALSE,
    params_timevar = params_timevar,
    step_args = list(do_hazard = TRUE),
    flexmodel = test_model
  )
  tmb_nd = Sys.time()
  r_strt = Sys.time()
  r_sim = run_sim(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    params_timevar = params_timevar,
    condense = FALSE,
    step_args = list(do_hazard = TRUE)
  )
  r_nd = Sys.time()
  tmb_speed = as.numeric(tmb_nd - tmb_strt)
  r_speed = as.numeric(r_nd - r_strt)
  print(paste0('tmb speed-up: ', round(r_speed/tmb_speed), 'x'))
  expect_lt(tmb_speed, r_speed)
  compare_sims(r_sim, tmb_sim)
})

test_that("time-varying parameters can be used with a vaccination model", {

  reset_spec_version()
  tmb_mode()
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

  params_timevar = data.frame(
    Date = c("2021-09-20"),
    Symbol = c('vax_prop_first_dose'),
    Value = c(0.5),
    Type = c("rel_orig")
  )

  test_model <- make_vaccination_model(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    params_timevar = params_timevar,
    step_args = list(do_hazard = TRUE)
  )

  tmb_strt = Sys.time()
  tmb_sim <- run_sim(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    condense = FALSE,
    params_timevar = params_timevar,
    step_args = list(do_hazard = TRUE),
    flexmodel = test_model
  )
  tmb_nd = Sys.time()
  r_strt = Sys.time()
  r_sim = run_sim(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    params_timevar = params_timevar,
    condense = FALSE,
    step_args = list(do_hazard = TRUE)
  )
  r_nd = Sys.time()
  tmb_speed = as.numeric(tmb_nd - tmb_strt)
  r_speed = as.numeric(r_nd - r_strt)
  print(paste0('tmb speed-up: ', round(r_speed/tmb_speed), 'x'))
  expect_lt(tmb_speed, r_speed)
  compare_sims(r_sim, tmb_sim)
})

test_that("toy symbolic matrix multiplication examples give correct results", {
  reset_spec_version()
  a = struc_block(struc("(a) * (1-b) + (1-a) * (1-b)"),
                  row_times = 1, col_times = 3)
  b = struc(letters[3:4])
  c = vec(LETTERS[1:6])
  set.seed(1)
  e = runif(10) %>% as.list %>% setNames(c(letters[1:4], LETTERS[1:6]))

  r1 = struc_eval(kronecker(a, t(b)) %*% c, e)
  r2 = struc_eval(a, e) %x% t(struc_eval(b, e)) %*% struc_eval(c, e)

  expect_equal(r1, r2)
})

test_that("symbolic expansion gives correct results", {

  expect_equal(
    struc('(a) + (b)') * struc('(c) + (d)'),
    struc('(a) * (c) + (b) * (c) + (a) * (d) + (b) * (d)')
  )

  expect_equal(
    vec('(a) + (b)', '(c) + (d)') * vec('(e) + (f)', '(g) + (h)'),
    vec('(a) * (e) + (b) * (e) + (a) * (f) + (b) * (f)',
        '(c) * (g) + (d) * (g) + (c) * (h) + (d) * (h)')
  )

  expect_equal(
    struc('(a) + (b)') * vec('c', 'd', 'e'),
    vec('(c) * (a) + (c) * (b)',
        '(d) * (a) + (d) * (b)',
        '(e) * (a) + (e) * (b)')
  )

  expect_equal(
    struc('(a) * (x) * (y) + (b) * (z)') * vec('c', '(d) * (w)', 'e'),
    vec('(c) * (a) * (x) * (y) + (c) * (b) * (z)',
        '(d) * (w) * (a) * (x) * (y) + (d) * (w) * (b) * (z)',
        '(e) * (a) * (x) * (y) + (e) * (b) * (z)')
  )

  expect_equal(
    kronecker(struc('(a) + (b)'), vec('c', 'd')),
    vec('(a) * (c) + (b) * (c)',
        '(a) * (d) + (b) * (d)')
  )

})

test_that("vax/variant foi algebraic manipulations are correct", {

  reset_spec_version()
  options(macpan_pfun_method = "grep")

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

  state = make_state(params = model_params)

  (epi_states = c(attr(state, "epi_cat")))
  (asymp_cat = c("S", "E", "Ia", "Ip", "R"))
  (vax_cat = c(attr(state, "vax_cat")))
  (dose_from = rep(asymp_cat, 2))
  (dose_to = c(asymp_cat, rep("V", 5)))
  (accum = c("X", "V"))
  (non_accum = base::setdiff(epi_states, accum))
  (non_accum_non_S = non_accum[-1])

  # Specify structure of the force of infection calculation
  Istate = (c('Ia', 'Ip', 'Im', 'Is')
            %>% McMasterPandemic:::expand_names(vax_cat)
            %>% vec
  )

  vax_trans_red = struc_block(vec(
    '1',
    '1',
    '(1 - vax_efficacy_dose1) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose1) * (variant_advantage) * (variant_prop)',
    '(1 - vax_efficacy_dose1) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose1) * (variant_advantage) * (variant_prop)',
    '(1 - vax_efficacy_dose2) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose2) * (variant_advantage) * (variant_prop)'),
    row_times = 1, col_times = 5)
  baseline_trans_rates =
    vec(
      'Ca',
      'Cp',
      '(1 - iso_m) * (Cm)',
      '(1 - iso_s) * (Cs)') *
    struc('(beta0) * (1/N)')


  e = c(as.list(state), as.list(model_params))

  r1 = struc_eval(kronecker(vax_trans_red, t(baseline_trans_rates)) %*% Istate, e)
  r2 = struc_eval(vax_trans_red, e) %x% t(struc_eval(baseline_trans_rates, e)) %*% struc_eval(Istate, e)

  expect_equal(r1, r2)
})

test_that("variants and vaccination model simulation matches run_sim", {

  reset_spec_version()
  #r_tmb_comparable()
  options(macpan_pfun_method = "grep")

  start_date <- "2021-02-01"
  end_date <- "2021-09-01"

  params_timevar = data.frame(
    Date = as.Date(c("2021-04-20", "2021-06-20", "2021-08-20", "2021-05-03", "2021-07-08")),
    Symbol = c("beta0", "beta0", "beta0", "Ca", "Ca"),
    Value = c(0.5, 0.1, 2.2, 0.2, 0.01),
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

  state = make_state(params = model_params)

  model = make_vaccination_model(
    params = model_params,
    state = state,
    start_date = start_date, end_date = end_date,
    do_hazard = TRUE,
    do_approx_hazard = FALSE,
    do_make_state = TRUE,
    max_iters_eig_pow_meth = 200,
    tol_eig_pow_meth = 1e-3,
    params_timevar = params_timevar,
    do_variant = TRUE)

  time_wrap(
    tmb_sim <- run_sim(
      params = model_params,
      state = state,
      start_date = start_date, end_date = end_date,
      params_timevar = params_timevar,
      condense = FALSE,
      step_args = list(do_hazard = TRUE),
      use_flex = TRUE,
      flexmodel = model
    ),
    r_sim <- run_sim(
      params = model_params,
      state = state,
      start_date = start_date, end_date = end_date,
      params_timevar = params_timevar,
      condense = FALSE,
      step_args = list(do_hazard = TRUE),
      use_flex = FALSE
    )
  )

  # useful diagnostics if tests start failing
  if(FALSE) {
    compare_sims_cbind(r_sim, tmb_sim, "Ia_unvax")
    compare_sims_plot(r_sim, tmb_sim, "Ia_unvax")
  }

  compare_sims(r_sim, tmb_sim, tolerance = 1e-9)

  # test state construction on c++ side matches r-side version
  expect_equal(c(state), initial_state_vector(model))

})

test_that("vax/variants model simulation matches run_sim with many break points", {

  # > reports_all$var %>% unique
  # [1] "report"
  # > reports_all$date %>% range
  # [1] "2020-02-04" "2021-11-28"
  # > est_date
  # [1] "2021-08-30"

  # calibration params used to produce an (initial?) set of
  # base parameters using read_params target argument

  reset_spec_version()
  #r_tmb_comparable()
  options(macpan_pfun_method = "grep")
  options(MP_rexp_steps_default = 200)

  # From macpan_ontario -----------------------------

  # key dates
  start_date <- as.Date("2020-09-15")
  end_date <- as.Date("2021-11-28")
  start_date_offset <- 60
  date_vec <- lubridate:::as_date(start_date:end_date)
  est_date = as.Date("2021-08-30")
  break_dates = structure(c(18520, 18538, 18589, 18626, 18641, 18685, 18705,
                            18725, 18789, 18808, 18824, 18857, 18871), class = "Date")
  bd = as.Date("2021-09-01")
  first_vax_date = as.Date("2020-12-14")

  # optimization parameters --------------
  opt_pars = list(time_params = 0.1)

  # simulation arguments ----------
  sim_args   = list(ndt = 1, step_args = list(do_hazard = TRUE))

  # optimizer arguments -----------
  mle2_control = list(maxit = 1e3,
                      reltol = 1e-8, #calibration_params[["reltol"]],
                      ## here beta and gamma are nelder-mead parameters,
                      ## not macpan model parameters!!!
                      beta = 0.5,
                      gamma = 2
  )

  # load data that came from macpan ontario
  load("../../inst/testdata/ontario_flex_test.rda")

  # make flex vaccination model -----------
  model = make_vaccination_model(
    params = model_params,
    state = NULL,
    start_date = start_date - days(start_date_offset),
    end_date = end_date + days(1),
    params_timevar = params_timevar,
    do_hazard = TRUE,
    do_hazard_lin = FALSE,
    do_approx_hazard = FALSE,
    do_approx_hazard_lin = FALSE,
    do_make_state = TRUE,
    max_iters_eig_pow_meth = 1000,
    tol_eig_pow_meth = 1e-4,
    do_variant = TRUE)

  time_wrap(
    tmb_sim <- run_sim(
      params = model_params,
      start_date = model$start_date,
      end_date = model$end_date,
      params_timevar = params_timevar,
      condense = FALSE,
      step_args = list(do_hazard = TRUE),
      flexmodel = model
    ),
    r_sim <- run_sim(
      params = model_params,
      start_date = model$start_date,
      end_date = model$end_date,
      params_timevar = params_timevar,
      condense = FALSE,
      step_args = list(do_hazard = TRUE)
    )
  )

  compare_sims(r_sim, tmb_sim, compare_attr = FALSE)

  # options(MP_rexp_steps_default = 100)
  #
  # fitted_mod_tmb <- calibrate(
  #   base_params = model$params,
  #   data = fitdat,
  #   opt_pars = opt_pars,
  #   start_date_offset = start_date_offset,
  #   debug = TRUE,
  #   time_args = list(
  #     params_timevar = params_timevar
  #   ),
  #   mle2_control = mle2_control,
  #   sim_args = list(
  #     ndt = 1,
  #     step_args = list(do_hazard = TRUE),
  #     flexmodel = model
  #   )
  # )

  # fitted_mod_r <- calibrate(
  #   base_params = model$params,
  #   data = fitdat,
  #   opt_pars = opt_pars,
  #   start_date_offset = start_date_offset,
  #   debug = TRUE,
  #   time_args = list(
  #     params_timevar = params_timevar
  #   ),
  #   mle2_control = mle2_control,
  #   sim_args = list(
  #     ndt = 1,
  #     step_args = list(do_hazard = TRUE)
  #   )
  # )

})
