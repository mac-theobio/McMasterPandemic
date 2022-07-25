# user-facing functions for getting and setting
# information about parameters

# generic getters ------------------------

#' Get and Set Model Parameters
#'
#' @note \code{pars_infer_*} are on the transformed scale that was chosen by the
#' user (e.g. log, logit).
#'
#' @param model \code{\link{flexmodel}} object
#' @return parameter vector or data frame of time-variation of parameters
#' @name get_and_set_model_parameters
NULL

#' @describeIn get_and_set_model_parameters base parameters that would be used in a simulation.
#' @export
pars_base_sim = function(model) {
  UseMethod('pars_base_sim')
}

#' @describeIn get_and_set_model_parameters base parameters that would be used as an initial value for an optimizer.
#' @export
pars_base_init = function(model) {
  UseMethod('pars_base_init')
}

#' @describeIn get_and_set_model_parameters optimized base parameters if they exist.
#' @export
pars_base_opt = function(model) {
  UseMethod('pars_base_opt')
}

#' @describeIn get_and_set_model_parameters time-varying parameter schedule that would be used in a simulation
#' @export
pars_time_sim = function(model) {
  UseMethod('pars_time_sim')
}

#' @describeIn get_and_set_model_parameters time-varying parameter schedule that was entered by the user at a model definition step.
#' @export
pars_time_spec = function(model) {
  UseMethod('pars_time_spec')
}

#' @describeIn get_and_set_model_parameters time-varying parameter schedule with optimized values that were entered as NA by the user in the model definition.
#' @export
pars_time_opt = function(model) {
  UseMethod('pars_time_opt')
}

#' @describeIn get_and_set_model_parameters time-varying parameter schedule with resolved multiplication strategy so that all values are on the same scale as the associated base parameter.
#' @export
pars_time_series = function(model) {
  UseMethod('pars_time_series')
}

#' @describeIn get_and_set_model_parameters the initial value of the parameter vector involved in inference (e.g. passed to an objective function during calibration).
#' @export
pars_infer_init = function(model) {
  UseMethod('pars_infer_init')
}

#' @describeIn get_and_set_model_parameters the optimized value of the parameter vector involved in inference.
#' @export
pars_infer_opt = function(model) {
  UseMethod('pars_infer_opt')
}

#' @rdname get_and_set_model_parameters
#' @export
pars_infer_conf = function(model) {
  UseMethod('pars_infer_conf')
}

# getter methods ------------------------

#' @exportS3Method
pars_base_sim.flexmodel = function(model) {
  as_vector_no_attr(model$params)
}

#' @exportS3Method
pars_base_init.flexmodel_to_calibrate = function(model) {
  p = model$params
  p[] = get_tmb_params(as.flexmodel_one_step(model))
  as_vector_no_attr(p)
}

#' @exportS3Method
pars_base_init.flexmodel_calibrated = function(model) {
  p = model$model_to_calibrate$params
  p[] = get_tmb_params(as.flexmodel_one_step(model$model_to_calibrate))
  as_vector_no_attr(p)
}

#' @exportS3Method
pars_base_opt.flexmodel_calibrated = function(model) {
  as_vector_no_attr(model$params)
}

#' @exportS3Method
pars_time_spec.flexmodel = function(model) {
  as_data_frame_no_row_names(get_params_timevar_orig(model))
}

#' @exportS3Method
pars_time_sim.flexmodel = function(model) {
  as_data_frame_no_row_names(get_params_timevar_impute(model))
}

#' @exportS3Method
pars_time_series.flexmodel = function(model) {
  (model
    %>% initialize_piece_wise
    %>% update_piece_wise(get_params_timevar_impute(model))
    %>% get_params_timevar_series
    %>% as_data_frame_no_row_names
  )
}

#' @exportS3Method
pars_time_opt.flexmodel_calibrated = function(model) {
  as_data_frame_no_row_names(get_params_timevar_orig(model))
}

#' @exportS3Method
pars_infer_init.flexmodel_to_calibrate = function(model) {
  as_vector_no_attr(tmb_params_trans(model))
}

#' @exportS3Method
pars_infer_init.flexmodel_calibrated = function(model) {
  as_vector_no_attr(tmb_params_trans(model$model_to_calibrate))
}

#' @exportS3Method
pars_infer_opt.flexmodel_calibrated = function(model) {
  as_vector_no_attr(model$opt_par)
}


# generic setters ------------------------

#' @param value named vector of parameters
#' @rdname get_and_set_model_parameters
#' @export
`pars_base_sim<-` = function(model, value) {
  UseMethod("pars_base_sim<-")
}

# @export
# `pars_base_init<-` = function(model) {
#   UseMethod('pars_base_init<-')
# }

# @export
# `pars_base_opt<-` = function(model) {
#   UseMethod('pars_base_opt<-')
# }

# @export
# pars_time_sim = function(model) {
#   UseMethod('pars_time_sim')
# }

# @export
# pars_time_spec = function(model) {
#   UseMethod('pars_time_spec')
# }

# @export
# pars_time_opt = function(model) {
#   UseMethod('pars_time_opt')
# }

# @export
# pars_time_series = function(model) {
#   UseMethod('pars_time_series')
# }

# @export
# pars_infer_init = function(model) {
#   UseMethod('pars_infer_init')
# }

# @export
# pars_infer_opt = function(model) {
#   UseMethod('pars_infer_opt')
# }

# @exportS3Method
# pars_infer_conf = function(model) {
#   UseMethod('pars_infer_conf')
# }


# setter methods -------------------------

#' @exportS3Method "pars_base_sim<-" flexmodel
`pars_base_sim<-.flexmodel` = function(model, value) {
  if (!all(names(value) %in% names(model$params))) {
    stop('only currently existing parameters can be set')
  }
  update_params(model, value)
}
