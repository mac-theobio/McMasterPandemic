# test script for developing tmbstan_tools
library(McMasterPandemic)
library(tidyr)

setup_stan()

# define SEIR model
set.seed(15)

# initial state, in units of individuals
state = c(
  S = 2000,
  E = 0,
  I = 1,
  R = 0
)

# parameters
params = c(
  beta = 0.7,
  N = sum(state),
  alpha = 0.05,
  gamma = 0.06
)

# start and end dates specified as strings,
# they don't have to be Date types
start_date = "2020-03-01"
end_date = "2020-06-01"

seir = flexmodel(
  params = params,
  state = state,
  start_date = start_date,
  end_date = end_date,
  do_hazard = TRUE # an option to ensure computations are stable when any of the state variables approach zero (keep this switched on)
)

seir = (
  seir
  %>% add_rate("S", "E", ~ (1/N) * (beta) * (I))
  %>% add_rate("E", "I", ~ (alpha))
  %>% add_rate("I", "R", ~ (gamma))
)

seir_obs_err = (seir
  # attach distribution parameters first
  %>% update_params(c(
    nb_disp_S = 1e2,
    nb_disp_E = 1e2,
    nb_disp_I = 1e2,
    nb_disp_R = 1e2
  ))
  # attach models for observation error distributions
  # using the parameters attached above
  %>% update_error_dist(
    S ~ negative_binomial("nb_disp_S"),
    E ~ negative_binomial("nb_disp_E"),
    I ~ negative_binomial("nb_disp_I"),
    R ~ negative_binomial("nb_disp_R")
  )
)

# generate synthetic data
seir_obs_err_result = (seir_obs_err
  # turn on observation error in simluation
  %>% simulation_history(obs_error = TRUE)
)
observed = (
  seir_obs_err_result
  %>% select(Date, I)
  %>% pivot_longer(-Date, names_to = "var")
  %>% rename(date = Date)
  # lob off first observation
  # (fitting the initial value is technically difficult)
  %>% slice(-1)
)

# set up model to calibrate
seir_obs_err_to_calibrate = (seir_obs_err
  # attach observed data
  %>% update_observed(
   observed
  )
  # attach priors for parameters we're fitting
  # ("optimizing" over)
  %>% update_opt_params(
   # fitting log beta
   log_beta ~ log_normal(
     -1, # log mean, so mean beta is exp(-1) = 0.36
     0.5 # standard deviation of the log normal prior
   ),
   # fitting log neg-binom dispersion for I
   log_nb_disp_I ~ log_normal(4, 1)
  )
)

# calibrate with stan
model_fit = calibrate_stan(
  model = seir_obs_err,
  model_to_calibrate = seir_obs_err_to_calibrate,
  chains = 2
)

# view traceplot
traceplot_stan(model_fit)

# forecast (status quo) with stan
fcst = forecast_stan(
  model_fit,
  parallel = TRUE,
  n_cores = 7
)

fcst_summary = summarise_forecast_stan(
  fcst,
  var_order = topological_sort(model_fit$model)
)

library(ggplot2)
(ggplot(fcst_summary,
        aes(x = date))
  + geom_ribbon(aes(ymin = lwr, ymax = upr, fill = var), alpha = 0.3)
  + geom_line(aes(y = value, colour = var), linewidth = 1.25)
  + facet_wrap(~ var, ncol = 1)
  + theme_bw()
)
