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
i = model$tmb_indices$initialization_mapping$all_to_eigen_idx
jac = report$j[i, i]
evec = eigen(jac)$vectors[,1]
evec / sum(evec)


report$eigenvec / sum(report$eigenvec)




lapply(model$tmb_indices, names)


model$tmb_indices$make_ratemat_indices$from
drop_pattern = '(S|Is)'
nms = c(names(state), names(params)) # TODO -- need to add sums here too
new_nms = grep(drop_pattern, nms, value = TRUE, invert = TRUE)
McMasterPandemic:::find_vec_indices(
  nms[model$tmb_indices$make_ratemat_indices$spi],
  setNames(new_nms, new_nms))

model$tmb_indices$outflow
model$tmb_indices$linearized_outflow

lapply(outflow, function(o) {
  list(
    state = grep(o$state_patterns, rownames(ratemat)),
    flow = grep(o$flow_state_patterns, colnames(ratemat))
  )
})

length(model$linearized_outflow)
length(model$tmb_indices$linearized_outflow)
# Reformat linearized_outflow

