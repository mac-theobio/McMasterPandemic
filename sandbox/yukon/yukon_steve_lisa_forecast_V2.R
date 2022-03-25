library(ggplot2)
library(McMasterPandemic)
library(lubridate)
library(tidyr)
library(dplyr)
library(ggplot2)
library(zoo)
if (dir.exists(lisa_dir <- "C:/Users/lkanary/OneDrive - Yukon College/Research/NSERC CANMOD/McMaster SEIR Model")) {
  setwd(lisa_dir)
}

#"in acute care" means "in hospital but not in the intensive care unit (ICU)"
params1 <- read_params("PHAC.csv")

# update default parameters for the Yukon case
#   beta0 -- baseline transmission rate
#   N -- population of Yukon
#   mu -- proportion of cases that are mild
#   phi1 -- proportion of hospitalization cases that are not in the ICU


# these are the suggested parameters from Yukon Gov Health and
# Social Services - they do some strange things to the fit...
params1[c("beta0", "N", "mu", "phi1")] <- c(0.25, 42507, 0.8, 0.22)



# initial state of the simulation
state1 <- make_state(params=params1)

# start and end dates
sdate <- as.Date("2021-11-01")
edate <- as.Date("2022-01-19")
initial_date = as.Date("2021-08-03")
start_date_offset = as.integer(sdate - initial_date)

#Steve, I wasn't sure how to smooth using the %>%
covidData <- read.csv(file = "report_data_yukon_h_and_i-1.csv")
smoothSpan <- 0.2
index <- 1:dim(covidData)[1]
covidData$report <- pmax(round(predict(loess(report ~ index, data=covidData, span=smoothSpan)), 0), 0)
covidData$R <- pmax(round(predict(loess(R ~ index, data=covidData, span=smoothSpan)), 0), 0)
covidData$hosp <- pmax(round(na.approx(covidData$hosp),0))
covidData$H <- pmax(round(na.approx(covidData$H),0))
covidData$ICU <- pmax(round(na.approx(covidData$ICU),0))

# read and process data
covid_data <- (covidData
  %>% mutate(date = as.Date(date))
  %>% filter(date >= ymd(20210803))
  %>% filter(between(as.Date(date), sdate, edate))
  %>% select(date, report, death, hosp, H, ICU) #  ICU)#, hosp) # R
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
    20211115, # nov 15
    20211215, # dec 15
    20220101  # jan 01 -- NEW: seems like transmission has declined again
  ),
  Symbol = "beta0",
  Value = NA,
  Type = "abs"
)

# manual exploration of parameter space
explore_pars = list(
  params = c(
    # baseline transmission rate
    beta0 = 0.235,
    # 1/time for severely symptomatic transition to hospital/death
    gamma_s = 1/14,
    # proportion of hospitalizations that are not in the ICU
    psi2 = 0.9,
    # Fraction of hospital cases to acute care
    phi1 = 0.05
  ),
  time_params = c(
    ## NOTE: this corresponds to the schedule
    ##       set in params_timevar above. any
    ##       parameter can vary but you need
    ##       to declare it in params_timevar.
    0.087, # beta0 changes to this value on nov 15
    0.4,  # beta0 changes to this value on dec 15
    0.25  # beta0 changes to this value on jan 01
  )
)
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
  +  facet_wrap(~var)
  +  geom_line(aes(date, value, colour = source))
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
# if you change the scale )
opt_pars <- list(
  params = c(

    # baseline transmission rate
    # (fit on log scale to avoid negative rates)
    log_beta0 = log(params1[["beta0"]]),

    # 1/time for severely symptomatic transition to hospital/death
    # (fit on logit scale to keep proportions between zero and one)
    log_gamma_s = log(params1[["gamma_s"]]),

    # proportion of hospitalizations that are not in the ICU
    # (fit on logit scale to keep proportions between zero and one)
    logit_psi2 = qlogis(params1[["psi2"]])
  ),


  # time varying transmissions rates
  # (see params_time_var above for schedule of changes in these rates)
  # (fit on log scale to avoid negative rates)
  log_time_params = rep(
    log(params1[["beta0"]]),
    nrow(params_timevar)
  )
)

# make use of c++ performance gains if you can
if (!is.null(getOption("MP_flex_spec_version"))) {
  mm = make_base_model(
    params = params1,
    start_date = sdate - start_date_offset,
    end_date = edate,
    params_timevar = params_timevar,
    do_make_state = TRUE,
    do_hazard = TRUE
  )
  sim_args = list(flexmodel = mm)
} else {
  sim_args = list()
}


# fit the model
fitted.mod <- calibrate(
  data = covid_data,
  time_args = list(params_timevar = params_timevar),
  start_date_offset = start_date_offset,
  base_params = params1,
  opt_pars = opt_pars,
  debug = TRUE,
  sim_args = sim_args
)

# fitted parameters, including dispersion parameters (nb_disp)
# that measure the variability around the fitted curves
coef(fitted.mod, 'fitted')
fitted.mod
plot(fitted.mod, data = covid_data)

# set up the forecasting scenario
# (this is a very basic status quo scenario for 30 days
# beyond the final date in the fitting data)
n_days_to_forecast = 30
forecast_args = fitted.mod$forecast_args
fit_end_date = as.Date(forecast_args$end_date)
forecast_args$end_date = fit_end_date + days(n_days_to_forecast)

# speed up forecasts by making use of c++ if it is available
if (!is.null(getOption("MP_flex_spec_version"))) {
  mm = make_base_model(
    params = params_forecast,
    start_date = mm$start_date,
    end_date = forecast_args$end_date,
    params_timevar = params_timevar,
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
