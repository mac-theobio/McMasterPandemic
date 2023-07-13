# reprex
# update_params not affecting simulation_history output when
# model object is a flexmodel_to_calibrate object??
# also, simulation_history(as.flexmodel(mod)) crashes rstudio when
# mod is a flexmodel_to_calibrate object

devtools::load_all() # load from local copy, new functions added

# base, deterministic model
sir = (flexmodel(
  params = c(beta = 0.1, gamma = 0.01, N = 100),
  state = c(S = 99, I = 1, R = 0),
  start_date = "2020-03-11",
  end_date = "2020-12-01"
)
  %>% add_rate("S", "I", ~ (I) * (beta) * (1/N))
  %>% add_rate("I", "R", ~ (gamma))
)

# sim with different betas
sim1 = simulation_history(sir)
sim2 = simulation_history(update_params(sir, c(beta = 0.5)))

n = nrow(sim1)

sim1[n, "S"] != sim2[n, "S"] # should be TRUE because we've changed beta

# now make it a flexmodel_to_calibrate object by attaching observed
# and opt params
obs = (sim1
 %>% dplyr::select(
   Date, I
 )
 %>% tidyr::pivot_longer(-Date, names_to = "var")
 %>% dplyr::rename(date = Date)
 %>% slice(-1)
)

# one prior value for beta
sir_to_calibrate1 = (sir
  %>% update_observed(
  obs
)
  # attach priors for parameters we're fitting
  # ("optimizing" over)
  %>% update_opt_params(
    # fitting log beta
    log_beta ~ log_normal(
      -1, # log mean, so mean beta is exp(-1) = 0.36
      0.5 # standard deviation of the log normal prior
    )
  )
)

# another prior for beta
sir_to_calibrate2 = (sir
   %>% update_observed(
     obs
   )
   # attach priors for parameters we're fitting
   # ("optimizing" over)
   %>% update_opt_params(
     # fitting log beta
     log_beta ~ log_normal(
       -0.5, # log mean, so mean beta is exp(-1) = 0.36
       0.5 # standard deviation of the log normal prior
     )
   )
)
sim3 = simulation_history(sir_to_calibrate1)
sim4 = simulation_history(sir_to_calibrate2)
testthat::expect_equal(sim1[n, "S"], sim3[n, "S"]) # should be equal if simulating from base parameter list
testthat::expect_equal(sim3[n, "S"], sim4[n, "S"]) # these are equal but have diff priors so that's not where pars are being drawn from?? where then??

sim5 = simulation_history(sir_to_calibrate1,
                          sim_params = sir_to_calibrate1$params)
testthat::expect_equal(sim1[n, "S"], sim5[n, "S"]) # should be equal if simulating from base parameter list
