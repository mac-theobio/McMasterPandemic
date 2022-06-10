# Convert between different types of flexmodels

#' @export
as.flexmodel = function(model) {
  UseMethod('as.flexmodel')
}

#' @export
as.flexmodel_to_calibrate = function(model) {
  UseMethod('as.flexmodel_to_calibrate')
}

#' @export
as.flexmodel_obs_error = function(model) {
  UseMethod('as.flexmodel_obs_error')
}

#' @export
as.flexmodel_one_step = function(model) {
  UseMethod('as.flexmodel_one_step')
}

#' @exportS3Method
as.flexmodel.flexmodel_to_calibrate = function(model) {
  as.flexmodel(as.flexmodel_obs_error(model))
}

#' @exportS3Method
as.flexmodel_obs_error.flexmodel_to_calibrate = function(model) {
  model$opt_params = NULL
  model$opt_tv_params = NULL
  model$observed = init_observed
  class(model) = c("flexmodel_obs_error", "flexmodel")
  model
}

#' @exportS3Method
as.flexmodel.flexmodel_obs_error = function(model) {
  reset_error_dist(model)
}

#' @exportS3Method
as.flexmodel_one_step.flexmodel = function(model) {
  model = initialize_piece_wise(model)
  model$observed$data = init_observed$data
  model = update_simulation_bounds(model, model$start_date, model$start_date)
  class(model) = c("flexmodel_one_step", "flexmodel")
  model
}
