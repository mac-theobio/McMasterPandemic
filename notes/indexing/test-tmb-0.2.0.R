library(McMasterPandemic)
library(dplyr)

set_spec_version("0.2.0", '../../inst/tmb')
tmb_mode()

params <- read_params("ICU1.csv")
model = (make_base_model(
    params = params,
    start_date = "2000-01-01",
    end_date = "2000-01-02",
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
