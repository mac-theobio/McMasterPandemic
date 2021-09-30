library(McMasterPandemic)
library(dplyr)
library(testthat)

test_that('foi can be expressed as model structure', {

  options(MP_flex_spec_version = "0.0.6")
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


library(McMasterPandemic)
library(dplyr)
library(testthat)

test_that('many rates can be expressed at the same time', {

  options(MP_flex_spec_version = "0.0.6")

  params <- read_params("ICU1.csv")
  state <- make_state(params = params)
  M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

  foi = sum(
    vec('Ia', 'Ip', 'Im', 'Is')
    * vec('Ca', 'Cp', '(Cm) * (1 - iso_m)', '(Cs) * (1 - iso_s)')
    * vec('beta0')
    * vec('1/N')
  )
  foi

  mm = (init_model(
    params, state,
    start_date = "2021-09-09", end_date = "2021-10-09",
  )
  %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
  %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
  %>% add_rate("Ia", "R", ~ (gamma_a))
  %>% add_rate("Im", "R", ~ (gamma_m))
  %>% add_rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
  %>% add_rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
  %>% add_rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
  %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
  %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
  %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
  %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
  %>% add_rate("ICUs", "H2", ~ (psi1))
  %>% add_rate("ICUd", "D", ~ (psi2))
  %>% add_rate("H2", "R", ~ (psi3))
  %>% add_rate("H", "R", ~ (rho))
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
