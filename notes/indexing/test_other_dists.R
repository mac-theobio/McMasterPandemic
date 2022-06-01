library(ggplot2)
library(McMasterPandemic)
library(dplyr)
library(tidyr)
library(lubridate)

set_spec_version("0.2.0", '../../inst/tmb')

state = c(S = 20000, I = 100, R = 0)
params = c(gamma = 0.06, beta = 0.15)
start_date = ymd(20000101)
end_date = ymd(20000501)
sir = make_sir_model(
  params = params, state = state,
  start_date = start_date,
  end_date = end_date
)

sir_with_obs_err = (sir
  %>% update_params(c(
    R_err = 1,
    S_err = 0.1,
    I_err = 1e4
  ))
  %>% update_error_dist(
    S ~ negative_binomial("S_err"),
    R ~ negative_binomial("R_err"),
    I ~ negative_binomial("I_err")
  )
)
sims = (sir_with_obs_err
  %>% simulation_history(obs_error = TRUE, include_initial_date = FALSE)
  %>% select(-N, -S_to_I)
  %>% rename(date = Date)
  %>% pivot_longer(-date, names_to = 'var')
  %>% mutate(var = factor(var, topological_sort(sir)))
)

(ggplot(sims)
  + facet_wrap(~var, scales = 'free')
  + geom_point(aes(date, value))
)

sir_with_other_dists = (sir
  %>% update_params(c(
    R_err = 1,
    S_err = 0.1,
    I_err = 1e4
  ))
  %>% update_observed(sims)
  %>% update_error_dist(
    S ~ normal("S_err"),
    R ~ negative_binomial("R_err"),
    I ~ beta("I_err", "R_err")
  )
  %>% update_opt_params(log_beta ~ log_normal(-1, 0.1))
)

sir_cal = calibrate_flexmodel(sir_with_other_dists)
sir_cal$params
sir_with_obs_err$params
sir_with_other_dists$params

ii = update_tmb_indices(sir_with_other_dists)$tmb_indices
ii$observed[c("variable_id", "loss_id", "spi_loss_param", "loss_param_count")]

(ggplot(mutate(fitted(sir_cal), var = factor(var, c("S", "I", "R"))))
  + facet_wrap(~var, scales = 'free')
  + geom_line(aes(date, value))
  + geom_line(aes(date, value_fitted), colour = 'red')
)
