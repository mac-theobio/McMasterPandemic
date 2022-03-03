# negative log likelihood

library(ggplot2)
library(McMasterPandemic)
library(lubridate)
library(tidyr)
library(dplyr)
set_spec_version("0.2.0", 'inst/tmb')
r_mode()
tmb_mode()


load_model_from_file = TRUE
options(MP_force_symm_vcov = TRUE)
options(MP_rexp_steps_default = 200)

# construct example ---------------------------

params1 <- read_params("PHAC.csv")
params1[c("N", "phi1")] <- c(42507, 0.98)
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
    20211115, # nov 15 beta0
    20211215, # dec 15 beta0
    20211215, # dec 15 mu
    20211217, # dec 17 mu
    20220101 # jan 01 beta0
  ),
  Symbol = c("beta0", "beta0", 'mu', 'mu', 'beta0'),
  Value = c(0.8, NA, NA, NA, NA),
  Type = "rel_prev"
)

mm = (make_base_model(
    params = params1,
    state = state1,
    start_date = sdate - start_date_offset,
    end_date = edate,
    params_timevar = params_timevar,
    do_hazard = TRUE,
    do_make_state = FALSE,
    data = covid_data
  )
  %>% update_opt_params(
    log_E0 ~ log_flat(log(5)),
    log_beta0 ~ log_normal(log(1), 1),
    log_mu ~ logit_normal(-0.04499737, 1),
    log_nb_disp_hosp ~ log_flat(0),
    log_nb_disp_report ~ log_flat(0)
  )
  %>% update_opt_tv_params(
    tv_type = 'rel_prev',
    log_beta0 + log_mu ~ log_flat(1)
  )
  %>% update_tmb_indices
)

# test gradients ---------------------------

compare_grads(mm, tolerance = 1e-5)
compare_grads(mm, tolerance = NA)

# optimization smell test ---------------------------

obj_fun = tmb_fun(mm)
opt = do.call("optim", obj_fun)
obj_fun$env$parList(opt$par)

# test objective function ---------------------------

op = get_opt_pars(mm$params, vars = c('hosp', 'report'))
op$params['log_E0'] = log(mm$params['E0'])
op$params['log_beta0'] = log(mm$params['beta0'])
op$params['log_mu'] = log(mm$params['mu'])
obj_fun$env$data  # items available to be read in to tmb using DATA_* macros

params_timevar2 = mutate(params_timevar, Value = c(0.8, 0.85, 0.9, 0.95, 0.99))
mm2 = (mm
  %>% update_piece_wise(params_timevar2)
  %>% update_opt_tv_params
  %>% update_tmb_indices
)
obj_fun2 = tmb_fun(mm2)
r_obj_fun = mle_fun(
  unlist(op),
  mm$observed$data,
  start_date = mm2$start_date,
  end_date = mm2$end_date,
  opt_pars = op,
  base_params = mm2$params,
  time_args = list(params_timevar = params_timevar2),
  priors = list(
    ~ dnorm(params[2], 0, 1),
    ~ dnorm(params[3], -0.04499737, 1)
  )
)

# use this test _before_ implementing transformations
# to get close to r_obj_fun
tmb_obj_fun_without_trans = obj_fun2$fn()
print(r_obj_fun)
print(tmb_obj_fun_without_trans)
print(all.equal(r_obj_fun, tmb_obj_fun_without_trans))

# use this test _after_ implementing transformations
# to match exactly with r_obj_fun
# NOTE: using this before implementing transformations
#       will result in NaN
#tmb_obj_fun_with_trans = obj_fun2$fn(tmb_params_init(mm2))
#all.equal(r_obj_fun, tmb_obj_fun_with_trans)
