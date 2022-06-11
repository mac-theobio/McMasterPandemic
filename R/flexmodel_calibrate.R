#' Wrapped Optimizers for Flexmodels
#'
#' Currently there are \code{optim_flexmodel}, \code{nlminb_flexmodel},
#' and \code{bbmle_flexmodel}
#'
#' @param model a \code{\link{flexmodel}} object
#' @param object a \code{\link{flexmodel}} object (for S3 method consistency)
#' @param ... additional arguments to pass to the wrapped optimizer
#' to their calibrated values?
#'
#' @export
optim_flexmodel = function(model, ...) {
  stopifnot(inherits(model, 'flexmodel_to_calibrate'))
  model_to_calibrate = model
  obj_fun = tmb_fun(model)
  model$opt_obj = optim(tmb_params_trans(model), obj_fun$fn, obj_fun$gr, ...)
  model$opt_par = model$opt_obj$par
  model$model_to_calibrate = model_to_calibrate
  class(model) = c('flexmodel_optim', 'flexmodel_calibrated', 'flexmodel')
  update_params_calibrated(model)
}

#' @rdname optim_flexmodel
#' @export
nlminb_flexmodel = function(model, ...) {
  stopifnot(inherits(model, 'flexmodel_to_calibrate'))
  model_to_calibrate = model
  obj_fun = tmb_fun(model)
  model$opt_obj = nlminb(tmb_params_trans(model), obj_fun$fn, obj_fun$gr, obj_fun$he, ...)
  model$opt_par = model$opt_obj$par
  model$model_to_calibrate = model_to_calibrate
  class(model) = c('flexmodel_nlminb', 'flexmodel_calibrated', 'flexmodel')
  update_params_calibrated(model)
}

#' @rdname optim_flexmodel
#' @export
bbmle_flexmodel = function(model, ...) {
  stopifnot(inherits(model, 'flexmodel_to_calibrate'))
  model_to_calibrate = model
  obj_fun = tmb_fun(model)
  start_par = tmb_params_trans(model)
  if (getOption("MP_get_bbmle_init_from_nlminb")) {
    start_par[] = nlminb_flexmodel(model)$opt_par
  }
  bbmle::parnames(obj_fun$fn) = names(start_par)
  bbmle::parnames(obj_fun$gr) = names(start_par)
  model$opt_obj = bbmle::mle2(
    obj_fun$fn,
    start_par,
    gr = obj_fun$gr,
    parnames = names(start_par),
    vecpar = TRUE,
    ...
  )
  model$opt_par = model$opt_obj@coef
  model$model_to_calibrate = model_to_calibrate
  class(model) = c('flexmodel_bbmle', 'flexmodel_calibrated', 'flexmodel')
  update_params_calibrated(model)
}

#' @exportS3Method
fitted.flexmodel_calibrated = function(object, ...) {
  model = object # for S3 consistency
  obs_var = unique(model$observed$data$var)
  fits = (model
    %>% simulation_history(obs_error = FALSE, condense = TRUE)
    %>% pivot_longer(-Date, names_to = "var")
    %>% filter(var %in% obs_var)
  )
  # fits = (model
  #   %>% simulate(do_condensation = TRUE)
  #   %>% filter(variable %in% obs_var)
  # )
  comparison_data = (model$observed$data
    %>% left_join(
      fits, by = c("date" = "Date", "var" = "var"),
      suffix = c('', '_fitted')
    )
  )
  comparison_data
}
