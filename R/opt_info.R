# getting information about calibrated models, objective functions,
# and models to calibrate

#' Optimizer Object
#'
#' Get the object returned by the optimizer used to calibrate a
#' \code{flexmodel_to_calibrate} object. The \code{convergence_info}
#' function gets convergence information in case it is buried within
#' the \code{opt_obj}.
#'
#' @param model object of class \code{flexmodel_calibrated}
#' @export
opt_obj = function(model) {
  stopifnot(inherits(model, 'flexmodel_calibrated'))
  model$opt_obj
}

#' @rdname opt_obj
#' @export
convergence_info = function(model) {
  UseMethod('convergence_info')
}

#' @exportS3Method
convergence_info.default = function(model) {
  stop('this function only applies to flexmodel_calibrated objects.')
}

#' @exportS3Method
convergence_info.flexmodel_calibrated = function(model) {
  return(opt_obj(model))
}

#' @exportS3Method
convergence_info.flexmodel_failed_calibration = function(model) {
  return(model$opt_err)
}

#' @exportS3Method
convergence_info.flexmodel_bbmle = function(model) {
  return(opt_obj(model)@details)
}


# getting information about objective functions and models to calibrate --------------

profile_obj_fun = function(model, focal_param) {
  UseMethod('profile_obj_fun')
}

profile_obj_fun.flexmodel_to_calibrate = function(model, focal_param) {
  model = remove_opt_param(model, focal_param)
  # o = tmb_fun(model)
  # param_names = names(tmb_params_trans(model))
  # par_ind = which(param_names == focal_param)
  # par = o$par
  function(x) {
    model$params[focal_param] = x
    o = tmb_fun(model)
    o$fn(par)
  }
}

# TODO: a similar user-facing function should refer to parameters
# for observation error distributions
loss_params = function(model) {
  model$observed$loss_params
}
