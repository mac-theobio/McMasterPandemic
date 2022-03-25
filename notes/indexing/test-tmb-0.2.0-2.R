# negative log likelihood

library(ggplot2)
library(McMasterPandemic)
library(lubridate)
library(tidyr)
library(dplyr)
set_spec_version("0.2.0", 'inst/tmb')
r_mode()
tmb_mode()

options(MP_force_symm_vcov = TRUE)
options(MP_rexp_steps_default = 200)

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


View(rate_summary(yukon_model))
yukon_fit = nlminb_flexmodel(yukon_model)

(yukon_fit
  %>% fitted
  %>% ggplot()
   +  facet_wrap( ~ var, scales = 'free')
   +  geom_point(aes(date, value))
   +  geom_line(aes(date, value_fitted))
)


.









obj_fun = tmb_fun(yukon_fit)
obj_fun$fn(yukon_fit$opt_par)
obj_fun$gr(yukon_fit$opt_par)
obj_fun$he(yukon_fit$opt_par)





# sh = simulation_history(yukon_fit, sim_params = yukon_fit$opt_par)

# compare objective function values --------------------------

op = list(
  params = c(log_beta0 = 0, logit_mu = -0.04499737),
  time_params = c(1, 1, 1, 1),
  log_nb_disp = c(hosp = 0, report = 0)
)

mle_fun(
  unlist(op),
  mm$observed$data,
  start_date = mm$start_date,
  end_date = mm$end_date,
  opt_pars = op,
  base_params = mm$params,
  time_args = list(params_timevar = params_timevar)
  #priors = list(
  #  ~ dnorm(params[2],  0,          1),
  #  ~ dnorm(params[3], -0.04499737, 1)
  #)
)
obj_fun = tmb_fun(mm)
obj_fun$fn()
opt = nlminb(obj_fun$par, obj_fun$fn, obj_fun$gr, obj_fun$he)

MASS::mvrnorm(10, opt$par, sdreport(obj_fun, opt$par)$cov.fixed)


obj_fun$fn(opt$par)
obj_fun$gr(opt$par)
obj_fun$he(opt$par)
obj_fun$env$
opt_par = obj_fun$env$parList(opt$par)
op$params[c('log_beta0', 'logit_mu')] = opt_par$params[c("beta0", "mu")]
op$time_params[] = opt_par$tv_mult[2:5]
op$log_nb_disp[c('hosp', 'report')] = opt_par$params[c("nb_disp_hosp", "nb_disp_report")]

mle_fun(
  unlist(op),
  mm$observed$data,
  start_date = mm$start_date,
  end_date = mm$end_date,
  opt_pars = op,
  base_params = mm$params,
  time_args = list(params_timevar = params_timevar)
  #priors = list(
  #  ~ dnorm(params[2],  0,          1),
  #  ~ dnorm(params[3], -0.04499737, 1)
  #)
)


# test gradients ---------------------------

compare_grads(mm, tolerance = 1e-5)
compare_grads(mm, tolerance = NA)

# optimization smell test ---------------------------

obj_fun = tmb_fun(mm)
opt = do.call("optim", obj_fun)
obj_fun$fn(opt$par)
obj_fun$env$parList(opt$par)
obj_fun$par
obj_fun$fn(opt$par)


# test objective function ---------------------------

params_timevar2 = mutate(params_timevar, Value = c(0.8, 0.85, 0.9, 0.95, 0.99))
mm2 = (mm
  %>% update_piece_wise(params_timevar2)
  %>% update_opt_tv_params
  %>% update_tmb_indices
)
obj_fun2 = tmb_fun(mm2)
opt2 = do.call("optim", obj_fun2)
opt2

op = get_opt_pars(mm$params, vars = c('hosp', 'report'))
op$params['log_E0'] = opt2$par[1]
op$params['log_beta0'] = opt2$par[2]
op$params['logit_mu'] = opt2$par[3]
op$params = op$params[names(op$params) != 'log_mu']
op$log_nb_disp['hosp'] = opt2$par[5]
op$log_nb_disp['report'] = opt2$par[4]

r_obj_fun = mle_fun(
  unlist(op),
  mm2$observed$data,
  start_date = mm2$start_date,
  end_date = mm2$end_date,
  opt_pars = op,
  base_params = mm2$params,
  time_args = list(params_timevar = params_timevar2)
  #priors = list(
  #  ~ dnorm(params[2],  0,          1),
  #  ~ dnorm(params[3], -0.04499737, 1)
  #)
)

# use this test _before_ implementing transformations
# to get close to r_obj_fun
#tmb_obj_fun_without_trans = obj_fun2$fn()
#print(r_obj_fun)
#print(tmb_obj_fun_without_trans)
#print(all.equal(r_obj_fun, tmb_obj_fun_without_trans))

# use this test _after_ implementing transformations
# to match exactly with r_obj_fun
# NOTE: using this before implementing transformations
#       will result in NaN
tmb_obj_fun_with_trans = opt2$value
all.equal(r_obj_fun, tmb_obj_fun_with_trans)


if (FALSE) {

  mm3 = (mm2
     %>% update_opt_params(
       log_beta0 ~ flat(1),
       log_mu ~ flat(0.956),
       log_nb_disp_hosp ~ flat(1),
       log_nb_disp_report ~ flat(1)
     )
     %>% update_tmb_indices
  )
  obj_fun3 = tmb_fun(mm3)
  opt_with_hess = nlminb(tmb_params_trans(mm3), obj_fun3$fn, obj_fun3$gr, obj_fun3$he)
  opt_wout_hess = optim(tmb_params_trans(mm3), obj_fun3$fn, obj_fun3$gr)
  tmb_params_trans(mm3)
  obj_fun3$gr()

  opt_with_hess$par
  opt_wout_hess$par
  opt_with_hess$objective
  opt_wout_hess$value

  op3 = op
  op3$params = op3$params[2:3]
  op3$params[] = tmb_params_trans(mm3)[1:2]
  names(op3) = c('params', 'nb_disp')
  names(op3$params) = c('beta0', 'mu')
  op3$nb_disp[] = 1
  op3$nb_disp = NULL


  tt = time_wrap(
    cal_r <- calibrate(
      data = covid_data,
      time_args = list(params_timevar = params_timevar2),
      start_date_offset = start_date_offset,
      base_params = mm3$params,
      opt_pars = op3,
      debug = FALSE
    ),
    cal_t <- calibrate(
      data = covid_data,
      time_args = list(params_timevar = params_timevar2),
      start_date_offset = start_date_offset,
      base_params = mm3$params,
      opt_pars = op3,
      debug = FALSE,
      sim_args = list(flexmodel = mm3)
    )
  )

}
