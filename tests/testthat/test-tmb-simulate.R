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

options(MP_default_do_sim_constraint = TRUE)

sir = (flexmodel(
    params = c(
      beta = 0.1, gamma = 0.01, N = 100,
      c_prop = 0.1, c_delay_mean = 1.1, c_delay_cv = 0.25
    ),
    state = c(S = 99, I = 1, R = 0),
    start_date = "2020-03-11",
    end_date = "2020-12-01"
  )
  %>% add_rate("S", "I", ~ (I) * (beta) * (1/N))
  %>% add_rate("I", "R", ~ (gamma))
  %>% add_sim_report_expr("incidence", ~ (S_to_I) * (S))
  %>% add_conv("incidence")
  %>% update_params(disp = 1)
  %>% add_error_dist(
    conv_incidence ~ negative_binomial("disp")
  )
)

synthetic_data = (sir
  %>% simulation_history(include_initial_date = FALSE, obs_error = TRUE)
  %>% select(Date, I)
  %>% pivot_longer(-Date, names_to = "var", values_to = "value")
  %>% rename(date = Date)
)

sir_to_calibrate = (sir
  %>% update_observed(synthetic_data)
  %>% update_error_dist(
    I ~ poisson()
  )
  %>% update_opt_params(
    log_beta ~ log_flat(-2),
    logit_gamma ~ logit_flat(-4)
  )
)

sir_calibrated = calibrate_flexmodel(sir_to_calibrate)

cbind(
  pars_base_sim(sir),
  pars_base_init(sir_calibrated),
  pars_base_opt(sir_calibrated)
)
