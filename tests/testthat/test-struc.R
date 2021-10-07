library(McMasterPandemic)
library(dplyr)
library(testthat)

test_that('foi can be expressed as model structure', {

  options(MP_flex_spec_version = "0.1.0")
  params <- read_params("ICU1.csv")
  state <- make_state(params = params)
  M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

  foi = sum(
      vec('Ia', 'Ip', 'Im', 'Is')
    * vec('Ca', 'Cp', '(Cm) * (1 - iso_m)', '(Cs) * (1 - iso_s)')
    * vec('beta0')
    * vec('1/N')
  )

  mm = (init_model(
    params, state,
    start_date = "2021-09-09", end_date = "2021-10-09",
  )
    %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
    %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
    %>% add_rate("Ia", "R", ~ (gamma_a))
    %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
    %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
    %>% add_rate("Im", "R", ~ (gamma_m))
    %>% add_rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
    %>% add_rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
    %>% add_rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
    %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
    %>% add_rate("ICUs", "H2", ~ (psi1))
    %>% add_rate("ICUd", "D", ~ (psi2))
    %>% add_rate("H2", "R", ~ (psi3))
    %>% add_rate("H", "R", ~ (rho))
    %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
    %>% add_rate("S", "E", foi)
    %>% add_parallel_accumulators(c("X", "N", "P", "V"))
    %>% add_tmb_indices()
  )

  tmb_sim <- run_sim(
    params = params, state = state,
    start_date = "2021-09-10",
    end_date = "2021-10-10",
    condense = FALSE,
    step_args = list(do_hazard = TRUE, flexmodel = mm),
    use_flex = TRUE
  )

  r_sim <- run_sim(
    params = params, state = state,
    start_date = "2021-09-10",
    end_date = "2021-10-10",
    condense = FALSE,
    step_args = list(do_hazard = TRUE),
    use_flex = FALSE
  )

  compare_sims(r_sim, tmb_sim)
})

spec_version <- "0.1.0"
print(spec_version)
options(MP_flex_spec_version = spec_version)

test_files <- "../../inst/tmb/"

# TODO: not yet using spec-version-specific dll because it is not yet created
#cpp <- file.path(test_files, spec_version, "macpan.cpp")
#dll <- file_path_sans_ext(cpp)
#options(MP_flex_spec_dll = basename(dll))
options(MP_flex_spec_dll = "McMasterPandemic")

#compile(cpp)
#dyn.load(dynlib(dll))

test_that('simple models still work when structure is allowed', {
  params <- read_params("ICU1.csv")
  state <- make_state(params = params)
  M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

  tv_dat <- data.frame(
    Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
    Symbol = c("beta0", "beta0", "beta0"),
    Value = c(0.5, 0.1, 0.05),
    Type = c("rel_prev", "rel_orig", "rel_prev")
  )

  mm = make_base_model(
    params, state,
    start_date = "2021-09-10",
    end_date = "2021-10-10",
    params_timevar = tv_dat,
    step_args = list(do_hazard = TRUE)
  )

  # change beta0, which is time-varying, so that
  # we can check that the parameter updates are
  # happening correctly in the C++ side
  test_pars = params
  test_pars[1] = 3

  tmb_sim <- run_sim(
    params = test_pars, state = state,
    start_date = "2021-09-10",
    end_date = "2021-10-10",
    params_timevar = tv_dat,
    condense = FALSE,
    step_args = list(do_hazard = TRUE, flexmodel = mm),
    use_flex = TRUE
  )

  r_sim <- run_sim(
    params = test_pars, state = state,
    start_date = "2021-09-10",
    end_date = "2021-10-10",
    params_timevar = tv_dat,
    condense = FALSE,
    step_args = list(do_hazard = TRUE),
    use_flex = FALSE
  )
  compare_sims(r_sim, tmb_sim)
})

test_that("simple vaccination model in TMB matches and is faster than existing R model", {
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

  test_model <- make_vaccination_model(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    step_args = list(do_hazard = TRUE)
  )

  tmb_strt = Sys.time()
  tmb_sim <- run_sim(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    condense = FALSE,
    step_args = list(do_hazard = TRUE),
    flexmodel = test_model,
    use_flex = TRUE
  )
  tmb_nd = Sys.time()
  r_strt = Sys.time()
  r_sim = run_sim(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
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

test_that("vax_prop_first_dose can be != 1", {
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
    step_args = list(do_hazard = TRUE)
  )

  tmb_strt = Sys.time()
  tmb_sim <- run_sim(
    params = vax_params, state = vax_state,
    start_date = "2021-09-09", end_date = "2021-10-09",
    condense = FALSE,
    params_timevar = params_timevar,
    step_args = list(do_hazard = TRUE),
    flexmodel = test_model,
    use_flex = TRUE
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
    flexmodel = test_model,
    use_flex = TRUE
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
