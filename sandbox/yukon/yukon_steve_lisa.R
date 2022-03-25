library(ggplot2)
library(McMasterPandemic)
library(lubridate)
library(tidyr)
library(dplyr)
if (dir.exists(lisa_dir <- "C:/Users/lkanary/OneDrive - Yukon College/Research/NSERC CANMOD/McMaster SEIR Model")) {
  setwd(lisa_dir)
}

load_model_from_file = TRUE
options(MP_force_symm_vcov = TRUE)

#"in acute care" means "in hospital but not in the intensive care unit (ICU)"
params1 <- read_params("PHAC.csv")

# update default parameters for the Yukon case
#   beta0 -- baseline transmission rate
#   N -- population of Yukon
#   mu -- proportion of cases that are mild
#   phi1 -- proportion of hospitalization cases that are not in the ICU
#
# note: some more of these values will be updated below
#       by running exploratory simulations and trajectory
#       matching 'by eye'
params1[c("N", "phi1")] <- c(42507, 0.98)

# initial state of the simulation
state1 <- make_state(params=params1)

# start and end dates
sdate <- as.Date("2021-11-01")
edate <- as.Date("2022-01-19")
initial_date = as.Date("2021-08-03")
start_date_offset = as.integer(sdate - initial_date)

# read and process data
covid_data <- ("report_data_yukon_h_and_i.csv"
  %>% read.csv
  %>% mutate(date = as.Date(date))
  %>% filter(date >= ymd(20210803))
  %>% filter(between(as.Date(date), sdate, edate))

  # report -- new reported cases on that day
  # hosp -- new hospital admissions on that day
  %>% select(date, report, hosp)

  %>% pivot_longer(names_to = "var", -date)
  %>% mutate(value=round(value))
)

head(covid_data, n=12)

# establish schedule of time variation of parameters
params_timevar = data.frame(
  Date = ymd(
    # estimate a new transmission rate on
    # these dates (i'm no expert but these
    # seemed to "work")
    20211115, # nov 15 beta0
    20211215, # dec 15 beta0
    20211215, # dec 15 mu
    20220101 # jan 01 beta0
  ),
  Symbol = c("beta0", "beta0", 'mu', 'beta0'),
  Value = NA,
  Type = "abs"
)

# NOTE: manual exploration of parameter space
#   - this is where i figured out good starting values for
#     mu and beta0
#   - two types of parameters: params and time_params
#   - params is a named vector with names referring to
#     the parameter to be optimized
#   - optionally change the scale on which optimization occur
#     by prefixing these names with log_ or logit_
#   - values of this vector are initial values fed to the
#     optimizer
#   - time_params is a vector of values in the order of
#     of the params_timevar schedule above
#     if you change the scale )
explore_pars = list(
  params = c(
    # baseline transmission rate
    beta0 = 0.235,
    # proportion of cases that are severe
    mu = 0.9
  ),
  time_params = c(
    ## NOTE: this corresponds to the schedule
    ##       set in params_timevar above. any
    ##       parameter can vary but you need
    ##       to declare it in params_timevar.
    0.07, 0.63, 0.99, 0.17
  )
)
# after simulations, plot the graph and check how good
# the fit is -- go back and adjust the parameter values
# if the fit is 'way off' -- repeat until you have good
# starting values to give to the opt_pars object below
(forecast_sim(
  unlist(explore_pars), explore_pars, base_params = params1,
  start_date = sdate - start_date_offset,
  end_date = edate,
  sim_args = list(condense = TRUE),
  time_args = list(params_timevar = params_timevar)
)
  %>% list(covid_data)
  %>% setNames(c("expected", "observed"))
  %>% bind_rows(.id = 'source')
  %>% filter(var %in% c("report", "hosp"))
  %>% ggplot
  +  facet_wrap(~var, scale = 'free')
  +  geom_line(aes(date, value, colour = source))
)

# declare parameters to be fitted
opt_pars <- explore_pars


# make use of c++ performance gains if you can
if (!is.null(getOption("MP_flex_spec_version"))) {
  mm = (make_base_model(
      params = params1,
      start_date = sdate - start_date_offset,
      end_date = edate,
      params_timevar = params_timevar,
      data = covid_data
    )
    #%>% update_condense_map(c(
    #  conv_Incidence = "report",
    #  lag_1_diff_X = 'hosp'
    #))
    %>% update_tmb_indices
  )
  sim_args = list(flexmodel = mm)
} else {
  sim_args = list()
}

if (load_model_from_file) {
  fitted.mod = readRDS(file = "fitted.mod.rda")
} else {
  # fit the model
  fitted.mod <- calibrate(
    data = covid_data,
    time_args = list(params_timevar = params_timevar),
    start_date_offset = start_date_offset,
    base_params = params1,
    opt_pars = opt_pars,
    debug = TRUE,
    debug_plot = FALSE,
    sim_args = sim_args,
    #mle2_args = list(use.ginv = FALSE)
  )
  saveRDS(fitted.mod, file = "fitted.mod.rda")
}

#fitted.mod.orig = fitted.mod
vc = fitted.mod$mle@vcov
#fitted.mod$mle@vcov = (vc + t(vc)) / 2

# fitted parameters, including dispersion parameters (nb_disp)
# that measure the variability around the fitted curves
coef(fitted.mod, 'fitted')
fitted.mod
plot(fitted.mod, data = covid_data)

# set up the forecasting scenario
# (this is a very basic status quo scenario for 30 days
# beyond the final date in the fitting data)
# (if you want to do more than status quo you will
# probably want to modify the schedule of parameter changes
# in forecast_args$time_args$params_timevar)
n_days_to_forecast = 30
forecast_args = fitted.mod$forecast_args
fit_end_date = as.Date(forecast_args$end_date)
forecast_args$end_date = fit_end_date + days(n_days_to_forecast)

# speed up forecasts by making use of c++ if it is available
if (!is.null(getOption("MP_flex_spec_version"))) {
  mm = make_base_model(
    params = forecast_args$base_params,
    start_date = mm$start_date,
    end_date = forecast_args$end_date,
    params_timevar = forecast_args$time_args$params_timevar,
    do_make_state = TRUE,
    do_hazard = TRUE
  )
  forecast_args$sim_args = list(flexmodel = mm)
}

# run the forecasts multiple times to get confidence intervals
#  - the result of this function has the following columns
#    - date
#    - var: variable being simulated
#    - lwr: 90% confidence interval lower bound
#    - value: median of the prediction distribution
#    - upr: 90% confidence interval upper bound
mod.forecasts = forecast_ensemble(
  fitted.mod,

  # forecast scenario is defined by forecast_args
  forecast_args = forecast_args
)

# here is a function for plotting the forecasts for a
# chosen variable (modified from work by Irena Papst).
# dashed horizontal line is where the data end.
plot_forecast = function(forecast, variable, data) {
  (forecast
    %>% filter(var == variable)
    %>% ggplot(aes(x = date))
    + geom_vline(aes(xintercept = fit_end_date),
                 colour = "grey60", linetype = "dashed")
    + geom_ribbon(aes(ymin = lwr,
                      ymax = upr),
                  alpha = 0.3)
    + geom_line(aes(y = value))
    + geom_point(data = filter(data, var == variable),
                 aes(x = date, y = value),
                 shape = 1)
    + labs(y = variable)
    + theme(axis.title.x = element_blank())
  )

}

# example plots
plot_forecast(mod.forecasts, "report", covid_data)
plot_forecast(mod.forecasts, "ICU", covid_data)
plot_forecast(mod.forecasts, "death", covid_data)
plot_forecast(mod.forecasts, "hosp", covid_data)

# if the variable isn't in the data set, then only the
# forecasts are plotted (no training data)
plot_forecast(mod.forecasts, "E", covid_data)
