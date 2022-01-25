library(McMasterPandemic)
library(dplyr)

set_spec_version("0.2.0", '../../inst/tmb')
tmb_mode()

params <- read_params("ICU1.csv")
model = (make_base_model(
    params = params,
    start_date = "2000-01-01",
    end_date = "2000-01-10",
  )
  %>% add_sim_report_expr("H_total", ~ (H) + (H2))
  %>% add_sim_report_expr("ICU", ~ (ICUs) + (ICUd))
  %>% add_sim_report_expr("Incidence", ~ (S_to_E) * (S))
  %>% add_lag_diff("^(X|D)$")
  %>% add_conv("^Incidence$")
  %>% update_tmb_indices
)
obj_fun = tmb_fun(model)

obj_fun$env$data  # items available to be read in to tmb using DATA_* macros
obj_fun$report()  # should contain the matrix with the simulation history
final_sim_report_names(model)  # R-side names of the columns that should be in the matrix with the simulation history
