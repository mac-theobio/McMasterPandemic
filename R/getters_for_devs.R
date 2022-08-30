# getting information about rates and sums ----------------------

# Get From-State Names
#
# Get vector of from-state names associated with each rate in the model.
#
# @param model \code{\link{flexmodel}} object
#
# @family get_flexmodel_info_functions

get_rate_from = function(model) get_rate_info(model, 'from')

# Get To-State Names
#
# Get vector of to-state names associated with each rate in the model.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_rate_to = function(model) get_rate_info(model, 'to')

# Get Rate Info
#
# Get information on the rate associated with a particular state transition.
#
# @param what character describing the rate (i.e. \code{'state1_to_state2'})
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_rate_info = function(model, what) {
  if (inherits(model, "flexmodel")) {
    rates = model$rates
  } else {
    rates = model
  }
  lapply(rates, '[[', what)
}

# Get Factr Info
#
# Get information about an intermediate factor
#
# @param what name of the intermediate factor, as named by
# \code{\link{add_factr}} or \code{\link{vec_factr}}
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_factr_info = function(model, what) lapply(model$factrs, '[[', what)

# Get Rate Formulas
#
# Get vector of formulas associated with each rate in the model.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_rate_formula = function(model) get_rate_info(model, 'formula')

# Get Rate Factors
#
# Get the factor table associated with each rate in the model.
# TODO: define the factor table or link to a description.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_rate_factors = function(model) get_rate_info(model, 'factors')

# Get Rate State Dependence
#
# Get a logical vector indicating which rates in the model depend on
# state variables.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_rate_state_dependent = function(model) get_rate_info(model, 'state_dependent')

# Get Rate Time Variation
#
# Get a logical vector indicating which rates in the model are
# time-varying.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_rate_time_varying = function(model) get_rate_info(model, 'time_varying')

# Get Rate Sum Dependence
#
# Get a logical vector indicating which rates in the model depend
# on sums of parameters and state variables.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_rate_sum_dependent = function(model) get_rate_info(model, 'sum_dependent')

# Get Factr Formula
#
# Get vector of formulas for defining each intermediate factor in the model.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_factr_formula = function(model) get_factr_info(model, 'formula')

# Get the Number of Products
#
# Get a vector giving the number of products (in the multiplication sense)
# that are required to compute each rate in the model.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_n_products = function(model) {
  (model
   %>% get_rate_info('factors')
   %>% lapply('[[', 'prod_indx')
   %>% lapply(unique)
   %>% lapply(length)
  )
}

# Get Number of Variables
#
# Get a vector giving the numbers of variables (parameters, state variables,
# sums of parameters and state variables, and intermediate factors) required
# to compute each rate in the model.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_n_variables = function(model) {
  (model
   %>% get_rate_info('factors')
   %>% lapply('[[', 'var')
   %>% lapply(length)
  )
}

# Get Number of Factors
#
# Get the number of factors required to compute each rate in the
# model.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_n_factors = function(model) {
  (model
   %>% get_rate_info('factors')
   %>% lapply(nrow)
  )
}

# Get Sum Info
#
# Get information on each sum (of parameters and state variables)
# in the model.
#
# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_sum_info = function(model, what) lapply(model$sums, '[[', what)

# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_sum_summands = function(model) get_sum_info(model, 'summands')

# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_sum_indices = function(model) get_sum_info(model, 'sum_indices')

# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_sum_initial_value = function(model) get_sum_info(model, 'initial_value')

# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_factr_initial_value = function(model) get_factr_info(model, 'initial_value')

# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_rate_vars = function(model) {
  (model
   %>% get_rate_factors
   %>% lapply(getElement, "var")
  )
}

# @family get_flexmodel_info_functions
# @inheritParams get_rate_from

get_rates_with_vars = function(model, var_pattern) {
  ii = (model
    %>% get_rate_vars
    %>% lapply(grepl, pattern = var_pattern)
    %>% lapply(any)
    %>% sapply(isTRUE)
  )
  get_rates(model)[ii]
}

get_rates = function(model) {
  model$rates
}

get_schedule = function(model) {
  model$timevar$piece_wise$schedule
}

get_params_timevar_orig = function(model) {
  (model
    %>% get_schedule
    %>% select(Date, Symbol, Value, Type)
  )
}

get_params_timevar_impute = function(model) {
  (model
    %>% get_schedule
    %>% select(Date, Symbol, init_tv_mult, Type)
    %>% rename(Value = init_tv_mult)
  )
}

get_params_timevar_series = function(model) {
  (model
   %>% get_schedule
   %>% select(Date, Symbol, tv_val)
   %>% rename(Value = tv_val)
  )
}

get_time_varying_baseline_params = function(model) {
  get_schedule(model)$Symbol
  tv_pars = unique(get_schedule(model)$Symbol)
  pars_base_sim(model)[tv_pars, drop = FALSE]
}

get_tmb_report = function(model) {
  tmb_fun(model)$report()
}

get_tmb_simulate = function(model) {
  tmb_fun(model)$simulate()
}

get_tmb_data = function(model) {
  tmb_fun(model)$env$data
}

get_tmb_params = function(model) {
  get_tmb_report(model)$params
}

get_tmb_tv_mult = function(model) {
  get_tmb_report(model)$tv_mult
}

get_tmb_hist = function(model) {
  get_tmb_report(model)$simulation_history
}

get_tmb_hist_stoch = function(model) {
  get_tmb_simulate(model)$simulation_history
}
