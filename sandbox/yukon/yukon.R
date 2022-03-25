# negative log likelihood

library(ggplot2)
library(McMasterPandemic)
library(lubridate)
library(tidyr)
library(dplyr)

# construct example ---------------------------

params <- read_params("PHAC.csv")
params[c("N", "phi1")] <- c(42507, 0.98)
params1 = params

state1 <- make_state(params=params1)

# start and end dates
sdate <- as.Date("2021-11-01")
edate <- as.Date("2022-01-19")
initial_date = as.Date("2021-08-03")
start_date_offset = as.integer(sdate - initial_date)

# read and process data
covid_data <- ("sandbox/yukon/report_data_yukon_h_and_i.csv"
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
    20211115, # nov 15 beta0 -- transmission rate
    20211215, # dec 15 beta0 -- transmission rate
    20211215, # dec 15 mu    -- prop mild cases
    20220101  # jan 01 beta0 -- transmission rate
  ),
  Symbol = c("beta0", "beta0", 'mu', 'beta0'),
  Value = c(NA, NA, NA, NA),
  Type = "rel_prev"
)

yukon_model = make_base_model(
    params = params1,
    state = state1,
    start_date = sdate - start_date_offset,
    end_date = edate,
    params_timevar = params_timevar,
    do_hazard = TRUE,
    do_make_state = FALSE,
    data = covid_data
  )

yukon_model = (yukon_model
  %>% update_opt_params(
    log_beta0 ~ log_flat(0),
    logit_mu ~ logit_flat(-0.04499737), # set to zero to see if it matters
    log_nb_disp_hosp ~ log_flat(0),
    log_nb_disp_report ~ log_flat(0)
  )
  %>% update_opt_tv_params(
    tv_type = 'rel_prev',
    log_beta0 ~ log_flat(0),
    log_mu ~ log_flat(0)
  )
)

yukon_fit = nlminb_flexmodel(yukon_model)

(yukon_fit
  %>% fitted
  %>% ggplot()
   +  facet_wrap( ~ var, scales = 'free')
   +  geom_point(aes(date, value))
   +  geom_line(aes(date, value_fitted))
)

obj_fun = tmb_fun(yukon_fit)
obj_fun$fn(yukon_fit$opt_par) # negative log posterior
obj_fun$gr(yukon_fit$opt_par) # gradient of the negative log posterior
obj_fun$he(yukon_fit$opt_par) # hessian of the negative log posterior

# full simulation history, including times without observed data,
# all states, all condensation summaries (e.g. lag_1_diff_X),
# and all time-varying rate-matrix elements (e.g. S_to_E)
sim_hist_yukon = simulation_history(
  yukon_fit,
  sim_params = yukon_fit$opt_par
)
