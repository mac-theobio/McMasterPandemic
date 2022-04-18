### THIS THING IS BROKEN

# negative log likelihood

library(ggplot2)
library(McMasterPandemic)
library(lubridate)
library(tidyr)
library(dplyr)

set_spec_version('0.2.0', 'inst/tmb')

r_tmb_comparable()

# construct example ---------------------------

params1 <- read_params("PHAC.csv")
params1[c("N", "phi1")] <- c(42507, 0.98)
(state1 <- make_state(params=params1))

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
    20211115, # nov 15 beta0
    20211215, # dec 15 beta0
    20220101 # jan 01 beta0
  ),
  Symbol = c("beta0", "beta0", 'beta0'),
  Value = c(NA, NA, NA),
  Type = "rel_prev"
)

mm = (make_base_model(
    params = params1,
    state = state1,
    start_date = sdate - start_date_offset,
    end_date = edate,
    params_timevar = params_timevar,
    do_hazard = TRUE,
    do_make_state = TRUE,
    tol_eig_pow_meth = 1e-03,
    data = covid_data
  )
  %>% update_opt_params(
    logit_mu ~ logit_flat(-0.04499737),
    log_beta0 ~ log_normal(log(1), 1),
    log_nb_disp_hosp ~ log_flat(0),
    log_nb_disp_report ~ log_flat(0)
  )
  %>% update_opt_tv_params(
    tv_type = 'rel_prev',
    log_beta0 ~ log_flat(0)
  )
  %>% update_tmb_indices
)

# compare objective function values --------------------------

obj_fun = tmb_fun(mm)
mm_fit = nlminb_flexmodel(mm, update_default_params = TRUE)
opt_par = mm_fit$opt_par
obj_fun$fn(opt_par) == mm_fit$opt_obj$objective

op = list(
  params = c(log_beta0 = opt_par[['log_beta0']], logit_mu = opt_par[['logit_mu']]),
  log_time_params = unname(opt_par[grep('^log_beta0_t[0-9]{3}$', names(opt_par))]),
  log_nb_disp = c(hosp = opt_par[['log_nb_disp_hosp']], report = opt_par[['log_nb_disp_report']])
)


# test objective function ---------------------------------

r_loss = mle_fun(
  unlist(op),
  mm$observed$data,
  start_date = mm$start_date,
  end_date = mm$end_date,
  opt_pars = op,
  base_params = mm$params,
  time_args = list(params_timevar = params_timevar),
  sim_args = list(step_args = list(do_hazard = TRUE), flexmodel = NULL),
  priors = list(
    ~ dnorm(log(params[1]),  0,   1)
  )
)
tmb_loss = mm_fit$opt_obj$objective
all.equal(unname(r_loss), unname(tmb_loss))

# test gradients ---------------------------

compare_grads(mm_fit, tolerance = 1e-5)
