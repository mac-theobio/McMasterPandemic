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

report <- tmb_fun(model)$report()
