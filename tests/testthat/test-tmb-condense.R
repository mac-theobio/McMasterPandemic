library(McMasterPandemic)
library(dplyr)


# version 0.2.0 will run without error, but will not create the required columns
# in the simulation history.
# version 0.2.1 with run with error until the 0.2.1 specs are implemented
#set_spec_version('0.2.1', system.file('tmb', package = 'McMasterPandemic'))
#set_spec_version('0.2.1', '../../inst/tmb')

uneven_dates_X = structure(
  c(10957, 10958, 10959, 10963, 10969, 10971, 10972,
  10973, 10975, 10977, 10980, 10981, 10984, 10986, 10987, 10988,
  10992, 10996, 11001, 11007, 11016, 11017, 11018, 11026, 11029,
  11030, 11039, 11042, 11044, 11045, 11055, 11056, 11057, 11061,
  11069, 11073, 11074, 11076, 11077), class = "Date")
uneven_dates_Y = structure(
  c(10962, 10969, 10971, 10972, 10977,
  11069, 11073, 11074, 11076, 11077), class = "Date")


state = c(S = 20000, I = 100, R = 0, X = 0, Y = 0)
sir_model = (
  flexmodel(
    params = c(
      gamma = 0.06,
      beta = 0.15,
      N = sum(state),
      c_prop = 1/10, c_delay_mean = 11, c_delay_cv = 0.25
    ),
    state = state,
    start_date = "2000-01-01",
    end_date = "2000-05-01",
    do_hazard = TRUE,
    do_make_state = FALSE
  )
  %>% add_piece_wise(data.frame(
    Date = "2000-03-01",
    Symbol = "c_prop",
    Value = 2/10,
    Type = "abs"
  ))
  %>% add_rate("S", "I", ~ (1/N) * (beta) * (I))
  %>% add_rate("S", "X", ~ (1/N) * (beta) * (I))
  %>% add_rate("I", "R", ~ (gamma))
  %>% add_rate("I", "Y", ~ (gamma))
  %>% add_outflow("S", "I")
  %>% add_outflow("I", "R")
  #%>% add_lag_diff("^X$", 1)
  %>% add_conv("^X$")
  %>% add_lag_diff_uneven("X", "X_uneven_diff", uneven_dates_X)
  %>% add_lag_diff_uneven("Y", "Y_uneven_diff", uneven_dates_Y)
  %>% update_tmb_indices
)

# xx = numeric(sir_model$iters)
# bb = (as.integer(difftime(uneven_dates_X, sir_model$start_date, units = "days")) + 1)[-1L]
# dd = as.integer(diff(uneven_dates_X))
# bb = bb[-length(bb)]
# xx[bb] = as.integer(diff(uneven_dates_X))

# uneven_dates_X

sir_model$tmb_indices$lag_diff$sri       # lag_diff_sri
View(sir_model$tmb_indices$lag_diff$delay_n)   # lag_diff_delay_n

View(simulation_history(sir_model)[,c("Date", "X", "Y", "X_uneven_diff", "Y_uneven_diff")])
View(simulation_history(extend_end_date(sir_model, 10))[,c("Date", "X", "lag_1_diff_X")])

McMasterPandemic:::uneven_from_even(
  extend_end_date(sir_model, 10),
  sir_model$lag_diff_uneven[[1]]$input_names,
  sir_model$lag_diff_uneven[[1]]$output_names
)


View(simulation_history(extend_end_date(sir_model, 10)))
