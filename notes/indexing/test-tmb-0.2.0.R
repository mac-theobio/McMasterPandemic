library(McMasterPandemic)
library(dplyr)

set_spec_version("0.2.0", '../../inst/tmb')
tmb_mode()

params <- read_params("ICU1.csv")
model = (make_base_model(
    params = params,
    start_date = "2000-01-01",
    end_date = "2000-03-01",
  )
  #%>% add_sim_report_expr("H_total", ~ (H) + (H2))
  #%>% add_sim_report_expr("ICU", ~ (ICUs) + (ICUd))
  #%>% add_sim_report_expr("Incidence", ~ (S_to_E) * (S))
  #%>% add_lag_diff("^(X|D)$")
  #%>% add_conv("^Incidence$")
  #%>% update_tmb_indices
)
obj_fun = tmb_fun(model)

obj_fun$env$data  # items available to be read in to tmb using DATA_* macros
obj_fun$report()  # should contain the matrix with the simulation history
final_sim_report_names(model)  # R-side names of the columns that should be in the matrix with the simulation history


his = (obj_fun
  $ report()
  $ simulation_history
  %>% as.data.frame
  %>% setNames(final_sim_report_names(model))
)
sim = simulate_state_vector(model)

all.equal(c(his[names(sim)]), c(sim))


r_sim = run_sim(
  params,
  start_date = model$start_date,
  end_date = model$end_date
)

all.equal(r_sim$incidence, his$Incidence)
all.equal(r_sim$foi, his$S_to_E)
all.equal(r_sim$H, his$Htotal)
all.equal(r_sim$ICU, his$ICU)
all.equal(r_sim$hosp[-1], his$lag_1_diff_X[-1])
all.equal(r_sim$death[-1], his$lag_1_diff_D[-1])
all.equal(
  his$conv_Incidence[!is.na(r_sim$report)],
  r_sim$report[!is.na(r_sim$report)]
)
plot(
  his$conv_Incidence[!is.na(r_sim$report)],
  r_sim$report[!is.na(r_sim$report)]
)

cbind(
  his$conv_Incidence,
  r_sim$report
) %>% View

gamma_shape <- 1 / delay_cv^2
gamma_scale <- delay_mean / gamma_shape
conv_tester = function(x, prop, delay_cv, delay_mean, qmax) {
  q = 1:qmax
  shape = 1 / (delay_cv^2)
  scale = delay_mean * (delay_cv^2)
  gamma_ = pgamma(q, shape = shape, scale = scale)
  delta = diff(gamma_)
  delta
  kappa = prop * (delta / sum(delta))
  kappa
  stats::filter(x, kappa, sides = 1)
}
cor(conv_tester(his$Incidence, params[["c_prop"]], params[["c_delay_cv"]], params[["c_delay_mean"]], 17)[-(1:15)],
  his$conv_Incidence[-c(1:15)])

with(
  as.list(params),
  McMasterPandemic:::make_delay_kernel(
    c_prop,
    c_delay_mean,
    c_delay_cv
  )
)
all.equal(conv_tester(his$Incidence, params[["c_prop"]], params[["c_delay_cv"]], params[["c_delay_mean"]], 17), kappa)

kappa = c(2.11639e-08,
3.58483e-06,
9.06584e-05,
0.000793502,
0.00352339,
0.00978145,
0.0192683,
0.0292384,
0.0361502,
0.0379025,
0.0347053,
0.0283761,
0.0210772,
0.0144171,
0.00918086,
0.00549147)




# negative log likelihood


library(ggplot2)
library(McMasterPandemic)
library(lubridate)
library(tidyr)
library(dplyr)
set_spec_version("0.2.0", '../../inst/tmb')
tmb_mode()

load_model_from_file = TRUE
options(MP_force_symm_vcov = TRUE)

params1 <- read_params("PHAC.csv")
params1[c("N", "phi1")] <- c(42507, 0.98)
state1 <- make_state(params=params1)

# start and end dates
sdate <- as.Date("2021-11-01")
edate <- as.Date("2022-01-19")
initial_date = as.Date("2021-08-03")
start_date_offset = as.integer(sdate - initial_date)

# read and process data
covid_data <- ("../../sandbox/yukon/report_data_yukon_h_and_i.csv"
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


mm = (make_base_model(
    params = params1,
    start_date = sdate - start_date_offset,
    end_date = edate,
    params_timevar = params_timevar,
    data = covid_data
  )
  %>% update_tmb_indices
)
sim_args = list(flexmodel = mm)
mm$tmb_indices$observed

obj_fun = tmb_fun(mm)
obj_fun$env$data  # items available to be read in to tmb using DATA_* macros
