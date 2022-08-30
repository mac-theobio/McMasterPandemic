#' Convert between Types of Flexmodels
#'
#' There are different types of \code{\link{flexmodel}} with structures that
#' are designed for specific purposes. There are often better ways to
#' create different types of \code{\link{flexmodel}}s, and these ways
#' are linked below in the Functions section.
#'
#' @param model \code{\link{flexmodel}} object
#' @return another \code{\link{flexmodel}} object of a different type
#' @export
as.flexmodel = function(model) {
  UseMethod("as.flexmodel")
}

#' @describeIn as.flexmodel If possible, convert a \code{\link{flexmodel}}
#' to one that can be calibrated. A more reliable way to generate a model
#' that can be calibrated is to configure what parameters should be
#' optimized using \code{\link{add_opt_params}} and associated functions,
#' and by adding observed data using \code{\link{update_observed}}.
#' @export
as.flexmodel_to_calibrate = function(model) {
  UseMethod("as.flexmodel_to_calibrate")
}

#' @describeIn as.flexmodel If possible, fit a \code{\link{flexmodel}}
#' to data using default optimization setting. To modify the defaults it is
#' better to use \code{\link{calibrate_flexmodel}}. Calibration is only
#' possible for models of class \code{flexmodel_to_calibrate}.
#' @export
as.flexmodel_calibrated = function(model) {
  UseMethod("as.flexmodel_calibrated")
}

#' @describeIn as.flexmodel If possible, convert a \code{\link{flexmodel}}
#' to one that has been configured to simulate observation error but cannot
#' be calibrated because it lacks parameters to be optimized. A more reliable
#' way to generate a model that can simulate observation error is to use
#' the \code{\link{add_error_dist}} and related functions.
#' @export
as.flexmodel_obs_error = function(model) {
  UseMethod("as.flexmodel_obs_error")
}

#' @describeIn as.flexmodel Convert a \code{\link{flexmodel}} into one that
#' takes one single step only. This is useful for getting information from
#' C++ that doesn't require running a full simulation (e.g. initial state
#' vector, back-transformed parameters).
#' @export
as.flexmodel_one_step = function(model) {
  UseMethod("as.flexmodel_one_step")
}

#' @describeIn as.flexmodel In possible, convert a \code{\link{flexmodel}}
#' object to an epidemic cohort model. A cohort model starts with a population
#' of a single exposed individual. The force of infection associated with that
#' initial exposed class can be interpreted as the kernel of the original model.
#' The kernel can then be used to calculate R0 and mean generation interval.
#' This conversion is only possible if the epidemiological summaries have been
#' configured using \code{\link{configure_epi_summaries}}.
#' @export
as.flexmodel_cohort = function(model) {
  UseMethod("as.flexmodel_cohort")
}

identity_as = function(model) model

#' @rdname as.flexmodel
#' @export
as.flexmodel.flexmodel = function(model) model

#' @rdname as.flexmodel
#' @export
as.flexmodel_to_calibrate.flexmodel_to_calibrate = identity_as

#' @rdname as.flexmodel
#' @export
as.flexmodel_calibrated.flexmodel_calibrated = identity_as

#' @rdname as.flexmodel
#' @export
as.flexmodel_calibrated.flexmodel_calibrated = identity_as

#' @rdname as.flexmodel
#' @export
as.flexmodel_obs_error.flexmodel_obs_error = identity_as

#' @rdname as.flexmodel
#' @export
as.flexmodel_one_step.flexmodel_one_step = identity_as

#' @rdname as.flexmodel
#' @export
as.flexmodel.flexmodel_to_calibrate = function(model) {
  as.flexmodel(as.flexmodel_obs_error(model))
}

#' @rdname as.flexmodel
#' @export
as.flexmodel_calibrated.flexmodel_to_calibrate = function(model) {
  calibrate_flexmodel(model)
}

#' @rdname as.flexmodel
#' @export
as.flexmodel.flexmodel_calibrated = function(model) {
  model$opt_obj = NULL
  model$opt_par = NULL
  model$model_to_calibrate = NULL
  class(model) = "flexmodel"
  model
}

#' @rdname as.flexmodel
#' @export
as.flexmodel.flexmodel_obs_error = function(model) {
  reset_error_dist(model)
}

#' @rdname as.flexmodel
#' @export
as.flexmodel.flexmodel_failed_calibration = function(model) {
  as.flexmodel(as.flexmodel_to_calibrate(model))
}

#' @rdname as.flexmodel
#' @export
as.flexmodel_to_calibrate.flexmodel_failed_calibration = function(model) {
  model$opt_err = NULL
  class(model) = c("flexmodel_to_calibrate", "flexmodel")
  model
}

#' @rdname as.flexmodel
#' @export
as.flexmodel_obs_error.flexmodel_to_calibrate = function(model) {
  model$opt_params = NULL
  model$opt_tv_params = NULL
  model$observed = init_observed
  class(model) = c("flexmodel_obs_error", "flexmodel")
  model
}

#' @rdname as.flexmodel
#' @export
as.flexmodel_one_step.flexmodel = function(model) {
  model = initialize_piece_wise(model)
  model$observed$data = init_observed$data
  model = update_simulation_bounds(model, model$start_date, model$start_date)
  class(model) = c("flexmodel_one_step", "flexmodel")
  model
}

#' @rdname as.flexmodel
#' @export
as.flexmodel_cohort.flexmodel = function(model) {
  unpack(model$summary_config)

  assert_len1_char(var_nms$exposed_state_nm)
  assert_len1_char(var_nms$foi_nm)

  model = (model
    %>% as.flexmodel  # make deterministic
    %>% initialize_piece_wise  # remove time-varying parameters
  )

  # start each simulation with a single exposed individual
  state_default(model) = 0
  state_default(model) = setNames(1, var_nms$exposed_state_nm)

  # update any parameters (e.g. N = 1)
  if (length(kernel_param_updates) > 0L) {
    par_updates = setNames(
      as.numeric(kernel_param_updates),
      names(kernel_param_updates)
    )
    pars_base_sim(model) = par_updates
  }

  # be sure to avoid using the eigenvector method of
  # initial state construction
  model$do_make_state = FALSE

  kernel_length = getOption("MP_kernel_len")
  end_date_extension = kernel_length - model$iters

  cohort_map = setNames(
    c('kernel', 'exposed'),
    c(var_nms$foi_nm, var_nms$exposed_state_nm)
  )

  model = (model
    %>% update_condense_map(cohort_map)
    %>% extend_end_date(end_date_extension)
  )
}
