library(McMasterPandemic)
library(tools)
library(dplyr)

cpp <- file.path("../../inst/tmb/0.1.1", "macpan.cpp")
dll <- file_path_sans_ext(cpp)
options(MP_flex_spec_dll = basename(dll))

compile(cpp)
dyn.load(dynlib(dll))

params <- read_params("ICU1.csv")
state <- make_state(params = params)
start_date = "2021-05-10"
end_date = "2021-12-10"
model <- make_base_model(
  params, state,
  start_date = start_date, end_date = end_date
)
model$do_hazard = FALSE
model$tmb_indices$outflow
model$tmb_indices$linearized_outflow

# Reformat linearized_outflow
n = length(model$tmb_indices$linearized_outflow)
model$tmb_indices$linearized_outflow_row_count = seq(n)
model$tmb_indices$linearized_outflow_col_count = seq(n)
model$tmb_indices$linearized_outflow_rows = vector()
model$tmb_indices$linearized_outflow_cols = vector()

for(i in 1:n) {
  rows = model$tmb_indices$linearized_outflow[[i]]$state
  cols = model$tmb_indices$linearized_outflow[[i]]$flow

  model$tmb_indices$linearized_outflow_row_count[i] = length(rows)
  model$tmb_indices$linearized_outflow_col_count[i] = length(cols)

  model$tmb_indices$linearized_outflow_rows = c(model$tmb_indices$linearized_outflow_rows, rows)
  model$tmb_indices$linearized_outflow_cols = c(model$tmb_indices$linearized_outflow_cols, cols)  
}

# Reformat outflow in the same way as we reformat linearized_outflow
n = length(model$tmb_indices$outflow)
model$tmb_indices$outflow_row_count = seq(n)
model$tmb_indices$outflow_col_count = seq(n)
model$tmb_indices$outflow_rows = vector()
model$tmb_indices$outflow_cols = vector()

for(i in 1:n) {
  rows = model$tmb_indices$outflow[[i]]$state
  cols = model$tmb_indices$outflow[[i]]$flow

  model$tmb_indices$outflow_row_count[i] = length(rows)
  model$tmb_indices$outflow_col_count[i] = length(cols)

  model$tmb_indices$outflow_rows = c(model$tmb_indices$outflow_rows, rows)
  model$tmb_indices$outflow_cols = c(model$tmb_indices$outflow_cols, cols)
}

report <- tmb_fun(model)$report()
report$j

