library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)
library(numDeriv)
library(lubridate)
library(tidyr)

set_spec_version('0.2.1', '../../inst/tmb')

uneven_dates = structure(
  c(10957, 10958, 10959, 10963, 10969, 10971, 10972,
  10973, 10975, 10977, 10980, 10981, 10984, 10986, 10987, 10988,
  10992, 10996, 11001, 11007, 11016, 11017, 11018, 11026, 11029,
  11030, 11039, 11042, 11044, 11045, 11055, 11056, 11057, 11061,
  11069, 11073, 11074, 11076, 11077), class = "Date")

state = c(S = 20000, I = 100, R = 0, X = 0)
sir_model = (
  flexmodel(
    params = c(
      gamma = 0.06,
      beta = 0.15,
      N = sum(state)
    ),
    state = state,
    start_date = "2000-01-01",
    end_date = "2000-05-01",
    do_hazard = TRUE,
    do_make_state = FALSE
  )
  %>% add_rate("S", "I", ~ (1/N) * (beta) * (I))
  %>% add_rate("S", "X", ~ (1/N) * (beta) * (I))
  %>% add_rate("I", "R", ~ (gamma))
  %>% add_outflow("S", "I")
  %>% add_outflow("I", "R")
  %>% add_lag_diff("^X$", 1)
  %>% add_lag_diff_uneven("X", "X_uneven_diff", uneven_dates)
  %>% update_tmb_indices
)
(sir_model
  %>% simulate
  %>% filter(variable %in% c("I", "R", "S", "X"))
  %>% ggplot()
   +  geom_line(aes(Date, value, colour = variable))
)

args(lag_diff_indices)
lag_diff_uneven_indices = function(model) {
  get_lag_info = function(x) {
    data.frame(
      i = which(simulation_dates(model) %in% x$lag_dates)[-1],
      x = as.integer(diff(x$lag_dates))
    )
  }
  get_sri_info = function(x) {
    sri = which(x$input_names == intermediate_sim_report_names(model))
  }
  d = (model$lag_diff_uneven
    %>% lapply(get_lag_info)
    %>% bind_rows(.id = "j")
  )
  sri = unlist(lapply(model$lag_diff_uneven, get_sri_info))
  # Matrix::sparseMatrix(
  #   i = d$i,
  #   j = d$j,
  #   x = d$x
  # )
  d
}
hh = lag_diff_uneven_indices(sir_model)
Matrix:::sparseMatrix(i = hh$i, j = hh$j, x = hh$x)
# lag_times = (sir_model
#   %>% simulate
#   %>% filter(variable == "X")
#   %>% filter(as.character(Date) %in% as.character(as_date(runif(50, sir_model$start_date, sir_model$end_date))))
#   %>% mutate(n = as.integer(Date - lag(Date)))
#   %>% select(Date, value, n)
# )
# (sir_model
#   %>% simulate
#   %>% left_join(lag_times, "Date")
# ) %>% View
#
# mutate(s, ifelse(variable == "lag_")) %>% View
#
# $lag_diff_uneven

