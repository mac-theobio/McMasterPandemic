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

#' @describeIn get_and_set_model_parameters base simulation parameters with time-varying parameters at the end of the simulation
#' @export
pars_base_final = function(model) {
  UseMethod('pars_base_final')
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

#' @describeIn get_and_set_model_parameters time-series for each time-varying parameter over the entire simulation range with resolved multiplication strategy.
#' @export
pars_time_hist = function(model) {
  UseMethod('pars_time_hist')
}

#' @describeIn get_and_set_model_parameters normalized time-series for each time-varying parameter over the entire simulation range with resolved multiplication strategy, with normalization obtained by dividing all values by the initial value.
#' @export
pars_time_norm = function(model) {
  UseMethod('pars_time_norm')
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
pars_base_final.flexmodel = function(model) {
  tv_value_end = function(d) d$Value[which.max(d$Date)]
  pars = pars_base_sim(model)
  ts = pars_time_series(model)
  final_tv_params = (ts
    %>% split(ts$Symbol)
    %>% vapply(tv_value_end, numeric(1L))
  )
  pars[names(final_tv_params)] = unname(final_tv_params)
  pars
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

#' @importFrom zoo na.locf
#' @exportS3Method
pars_time_hist.flexmodel = function(model) {
  sd = simulation_dates(model)
  ts = pars_time_series(model)
  base_tv_pars = get_time_varying_baseline_params(model)
  pars = names(base_tv_pars)
  ll = data.frame(Date = sd)
  for (i in seq_along(pars)) {
    ll[[pars[i]]] = NA
    ll[1, pars[i]] = base_tv_pars[i]
    w_ts = ts$Symbol == pars[i]
    w_ll = ll$Date %in% ts$Date[w_ts]
    ll[w_ll, pars[i]] = ts$Value[w_ts]
    ll[[pars[i]]] = zoo::na.locf(ll[[pars[i]]])
  }
  as_data_frame_no_row_names(ll)
}

#' @exportS3Method
pars_time_norm.flexmodel = function(model) {
  pars_norm = pars_time_hist(model)
  if (ncol(pars_norm) == 1L) return(pars_norm)
  pars_tv = get_time_varying_baseline_params(model)
  # TODO: warn if any baseline parameters are numerically zero
  pars_norm[,-1] = sweep(
    pars_norm[,-1, drop = FALSE],  # remove date column
    2L,
    pars_tv,
    "/"
  )
  as_data_frame_no_row_names(pars_norm)
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
