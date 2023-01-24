#' Simulation History
#'
#' @param model a \code{\link{flexmodel}} object
#' @param add_dates should a column called \code{Date} be added?
#' @param sim_params vector to pass to a TMB AD function -- you can get
#' examples of these out of the \code{opt_par} element of a fitted
#' \code{\link{flexmodel}} or by calling \code{\link{tmb_params_trans}}
#' @param include_initial_date should the first row be the initial state
#' vector, or the first date after the initial?
#' @param obs_error should random observation error be added if it is
#' defined in the \code{model}? -- see \code{\link{update_error_dist}}
#' @param condense should a condensed set of simulation variables be
#' returned? -- see \code{\link{update_condense_map}}
#' @param summaries should epidemiological summaries (including \code{Rt} and
#' the \code{rel_trans_rate}) be computed? note that this is only
#' possible if the epidemiological summaries have been configured using
#' \code{\link{configure_epi_summaries}}. the relative transmission rate
#' is the transmission rate at each time divided by the transmission rate at
#' the beginning of the simulations and computed using
#' \code{\link{pars_time_norm}}. note also that the effective
#' proportion of susceptible individuals will also be in the simulation
#' history (with a user-defined name) if the summaries have been configured.
#'
#' @family simulation
#' @export
simulation_history = function(
    model,
    add_dates = TRUE,
    sim_params = NULL,
    include_initial_date = TRUE,
    obs_error = FALSE,
    condense = FALSE,
    summaries = FALSE
  ) {
  if (is.null(sim_params)) sim_params = tmb_params(model)
  if (obs_error) {
    sim_hist_raw = tmb_fun(model)$simulate()$simulation_history
  } else {
    sim_hist_raw = tmb_fun(model)$report()$simulation_history
  }

  if (summaries) {
    unpack(model$summary_config)
    assert_len1_char(var_nms$trans_rate_nm)
    assert_len1_char(var_nms$prop_susceptible_nm)
    kern = epi_kernel(model)
    R0 = epi_R0(kern)
    rel_trans_rate = pars_time_norm(model)[[var_nms$trans_rate_nm]]
    uncondensed_prop_S = (
      (model$condensation_map == var_nms$prop_susceptible_nm)
      %>% which
      %>% names
    )
  }

  sim_hist = (sim_hist_raw
   %>% as.data.frame
   %>% setNames(final_sim_report_names(model))
  )

  if (summaries) {
    prop_S = sim_hist[[uncondensed_prop_S]]
    Rt = R0 *  prop_S * rel_trans_rate
  }

  # FIXME: code smell -- should have different lag diff classes
  if (spec_ver_gt("0.2.0")) {
    ld = model$lag_diff_uneven
  } else {
    ld = model$lag_diff
  }
  sim_hist = (sim_hist
    %>% pad_lag_diffs(ld)
    %>% pad_convs(model$conv, model$tmb_indices$conv)
  )

  if (condense) {
    sim_hist = setNames(
      sim_hist[names(model$condensation_map)],
      model$condensation_map
    )
    # FIXME: need cumRep hack in condense_flexmodel?
  }

  if (add_dates) {
    sim_hist = (model
      %>% simulation_dates
      %>% data.frame
      %>% setNames("Date")
      %>% cbind(sim_hist)
    )
  }

  sim_hist = sim_hist[seq_len(model$iters + 1L),]

  if (summaries) {
    if (any(c("rel_trans_rate", "Rt") %in% names(sim_hist))) {
      stop(
        "\nat least one simulation variable has one of the following",
        "\nreserved names: rel_trans_rate, Rt"
      )
    }
    sim_hist$rel_trans_rate = rel_trans_rate
    sim_hist$Rt = Rt
  }

  if (!include_initial_date) {
    if (nrow(sim_hist) == 1L) {
      warning("No simulations to return. Check start_date and end_date.")
    }
    sim_hist = sim_hist[-1,]
  }

  sim_hist
}


#' Ensemble Simulation
#'
#' Given a sample of parameter vectors, simulate the condensed state simulation
#' history for each ... TODO
#'
#' @param model \code{\link{flexmodel}} object, preferably a
#' \code{flexmodel_calibrated} object that was fitted with
#' \code{\link{bbmle_flexmodel}} or \code{\link{calibrate_flexmodel}}
#' @param sim_params_matrix numberic matrix with rows giving simulation
#' iterations and columns giving model parameters
#' @param covmat not used
#' @param qvec named numeric vector giving quantiles to compute
#' @param use_progress_bar should a progress bar be printed?
#' @param ... arguments to pass on the \code{pop_pred_samp}, which
#' is actually doing the sampling from the distribution of parameters.
#' useful parameters here include the sample size, \code{n} (default 1000),
#' and \code{PDify = TRUE} that can solve issues with non-positive definite
#' matrices.
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @family simulation
#' @export
simulate_ensemble = function(
    model,
    sim_params_matrix = NULL,
    covmat = NULL,
    qvec = c(value = 0.5, lwr = 0.025, upr = 0.975),
    use_progress_bar = TRUE,
    ...
  ) {
  if (!is.null(covmat)) {
    stop("covmat argument is not currently being used")
  }

  msg = paste0(
    "\nsim_params_matrix missing, and unable to produce it. ",
    "try using bbmle_flexmodel when calibrating.",
    collapse = "\n"
  )
  if (inherits(model, "flexmodel_calibrated")) {
    if (is.null(sim_params_matrix)) {
      if (is_fitted_by_bbmle(model)) {
        sim_params_matrix = pop_pred_samp(
          model$opt_obj,
          Sigma = vcov(model),
          ...
        )
      } else {
        stop(msg)
      }
    }
    model = model$model_to_calibrate
  }
  if (is.null(sim_params_matrix)) {
    stop(msg)
  }

  stopifnot(!is.null(names(qvec)))

  # avoid the cost of computing the loss function,
  # and just do the simulations
  model$observed$data = init_observed$data

  o = tmb_fun(model)
  if (length(o$par) != ncol(sim_params_matrix)) {
    warning(
      "You might not be able to fit as many parameters as you are\n",
      "or maybe there are multiple changes to a time-varying parameter\n",
      "on the same day.")
  }
  trajectories = list()
  ii = seq_len(nrow(sim_params_matrix))

  cond_map = model$condensation_map
  cond_nms = names(cond_map)

  date_frame = (model
    %>% simulation_dates
    %>% data.frame
    %>% setNames("Date")
  )
  fsrn = final_sim_report_names(model)

  if (use_progress_bar) {
    pb = txtProgressBar(min = min(ii), max = max(ii), initial = min(ii), style = 3)
  }

  for(i in ii) {
    traj = as.data.frame(o$simulate(sim_params_matrix[i,])$simulation_history)
    names(traj) = fsrn
    traj = setNames(
      traj[cond_nms],
      cond_map
    )
    trajectories[[i]] = cbind(date_frame, traj)
    if (use_progress_bar) {
      setTxtProgressBar(pb, i)
    }
  }

  if (use_progress_bar) cat("", "summarising ensemble ... ", sep = "\n")
  names(trajectories) = ii
  summarised_traj = (trajectories
    %>% bind_rows(.id = 'simulation')
    %>% pivot_longer(c(-simulation, -Date), names_to = 'var')
    %>% group_by(Date, var)
    %>% do(setNames(data.frame(t(quantile(.$value, probs = qvec, na.rm = TRUE))), names(qvec)))
  )
  attr(summarised_traj, "qvec") = qvec
  summarised_traj
}

#' Deprecated Simulation Functions
#'
#' @family simulation
#' @name deprecated_simulation_functions
NULL

#' @param model flexmodel
#' @param sim_params parameter vector to pass to a TMB objective function
#' @rdname deprecated_simulation_functions
#' @export
changing_ratemat_elements = function(model, sim_params = NULL) {
  if(is.null(sim_params)) sim_params = tmb_params(model)
  tmb_fun(model)$simulate(sim_params)$concatenated_ratemat_nonzeros
}

#' @param add_dates should column with dates be added to the output
#' @rdname deprecated_simulation_functions
#' @export
simulate_changing_ratemat_elements = function(model, sim_params = NULL, add_dates = FALSE) {
  if (getOption("MP_auto_outflow")) {
     model = add_outflow(model)
  }
  if (getOption("MP_auto_tmb_index_update")) {
    model = update_tmb_indices(model)
  }
  updateidx = model$tmb_indices$updateidx
  ratemat_elements = (model
    %>% changing_ratemat_elements(sim_params)
    %>% matrix(nrow = model$iters + 1L,
               ncol = length(updateidx),
               byrow = TRUE)
  )

  colnames(ratemat_elements) = names(updateidx)
  ratemat_elements = as.data.frame(ratemat_elements)
  if(add_dates) {
    ratemat_elements = (model
                  %>% simulation_dates
                  %>% data.frame
                  %>% setNames("Date")
                  %>% cbind(ratemat_elements)
    )
  }
  ratemat_elements
}

#' @rdname deprecated_simulation_functions
#' @export
concatenated_state_vector = function(model, sim_params = NULL) {
  if(is.null(sim_params)) sim_params = tmb_params(model)
  tmb_fun(model)$simulate(sim_params)$concatenated_state_vector
}

structure_state_vector = function(x, iters, state_nms) {
  matrix(
    x,
    nrow = iters + 1L,
    ncol = length(state_nms),
    dimnames = list(1:(iters + 1), state_nms),
    byrow = TRUE
  )
}

#' @param format how to return the results
#' @rdname deprecated_simulation_functions
#' @export
simulate_state_vector = function(model, sim_params = NULL, add_dates = FALSE,
                                 format = c('wide', 'long', 'ggplot')) {
  format = match.arg(format)
  state_sims = (model
   %>% concatenated_state_vector(sim_params)
   %>% structure_state_vector(model$iters, names(model$state))
   %>% as.data.frame
  )
  if(add_dates | format == 'ggplot') {
    state_sims = (model
      %>% simulation_dates
      %>% data.frame
      %>% setNames("Date")
      %>% cbind(state_sims)
    )
    if (format %in% c('long', 'ggplot')) {
      state_sims = pivot_longer(
        state_sims, !Date,
        names_to = 'Compartment',
        values_to = 'State'
      )
    }
  } else if(format == 'long') {
    state_sims = pivot_longer(
      state_sims,
      names_to = 'Compartment',
      values_to = 'State'
    )
  }
  if (format == 'ggplot') {
    state_sims = (state_sims
      %>% ggplot
       +  geom_line(aes(x = Date, y = State, colour = Compartment))
    )
  }
  state_sims
}

#' @rdname deprecated_simulation_functions
#' @export
initial_state_vector = function(model, sim_params = NULL) {
  # FIXME: minor performance optimization could be made
  # here by restricting the number of iterations temporarily
  # so that the full simulation is not run just to get the
  # first value -- i think all that is needed is this (but
  # need to test):
  # model$iters = 1

  # FIXME: should this be deprecated or just made to be not
  # user facing. the as.flexmodel_one_step is basically what
  # we want here

  (model
   %>% concatenated_state_vector(sim_params)
   %>% head(length(model$state))
   %>% setNames(names(model$state))
  )
}

#' @rdname deprecated_simulation_functions
#' @export
final_state_vector = function(model, sim_params = NULL) {
  (model
   %>% concatenated_state_vector(sim_params)
   %>% tail(length(model$state))
   %>% setNames(names(model$state))
  )
}

#' @rdname deprecated_simulation_functions
#' @export
penultimate_state_vector = function(model, sim_params = NULL) {
  (model
   %>% concatenated_state_vector(sim_params)
   %>% tail(2 * length(model$state))
   %>% head(length(model$state))
   %>% setNames(names(model$state))
  )
}


#' @rdname deprecated_simulation_functions
#' @export
final_state_ratio = function(model, sim_params = NULL) {
  last_two = (model
   %>% concatenated_state_vector(sim_params)
   %>% tail(2 * length(model$state))
  )
  n = length(model$state)
  setNames(
    head(last_two, n) / tail(last_two, n),
    names(model$state))
}

#' @rdname deprecated_simulation_functions
#' @export
initial_ratemat = function(model, sim_params = NULL) {
  if(is.null(sim_params)) sim_params = tmb_params(model)
  tmb_fun(model)$simulate(sim_params)$ratemat
}


# @param cond_map condensation_map element of flexmodel
# @param sim_hist value of simulation_history in a TMB report or simulation
# apply_condensation = function(cond_map, sim_hist) {


#' @rdname deprecated_simulation_functions
#' @export
simulation_condensed = function(model, add_dates = TRUE, sim_params = NULL) {
  cond_map = model$condensation_map
  cond_nms = names(cond_map)
  if (add_dates) {
    cond_nms = c("Date", cond_nms)
    cond_map = c(Date = "Date", cond_map)
  }
  sims = simulation_history(model, add_dates = TRUE, sim_params = sim_params)[cond_nms]
  names(sims) = c(cond_map)
  sims
}


#' @rdname deprecated_simulation_functions
#' @export
condense_flexmodel = function(model) {
  spec_check(
    introduced_version = '0.2.0',
    feature = "condensation in c++"
  )
  condensed_simulation_history = setNames(
    simulation_history(model)[names(model$condensation_map)],
    model$condensation_map
  )

  # HACK! ultimately we want cumulative reports calculated
  # on the c++ side (https://github.com/mac-theobio/McMasterPandemic/issues/171)
  # also this assumes no observation error, and doesn't
  # compute D as cumulative sum of deaths (as is done in run_sim)
  if ('report' %in% names(condensed_simulation_history)) {
    condensed_simulation_history$cumRep = cumsum(
      ifelse(
        !is.na(unname(unlist(condensed_simulation_history$report))),
        unname(unlist(condensed_simulation_history$report)),
        0
      )
    )
  }

  cbind(data.frame(date = simulation_dates(model), condensed_simulation_history))
}


#' @param object \code{\link{flexmodel}} object
#' @param nsim number of simulations
#' @param seed random seed
#' @param do_condensation should condensed set of variables be returned if
#' available
#' @param ... pass arguments on
#' @importFrom tidyr pivot_longer
#' @importFrom stats simulate
#' @rdname deprecated_simulation_functions
#' @export
simulate.flexmodel = function(
  object, nsim = 1, seed = NULL,
  do_condensation = FALSE,
  format = c('long', 'wide'),
  add_dates = TRUE,
  sim_params = NULL, ...) {

  format = match.arg(format)
  if ((nsim != 1) | !is.null(seed)) {
    stop(
      "nsim and seed cannot be set to non-default values yet.\n",
      "they may become relevant in the future if it becomes possible\n",
      "to add stochasticity to the simulations"
    )
  }
  if (do_condensation) {
    sims = simulation_condensed(object, add_dates, sim_params)
  } else {
    sims = simulation_history(object, add_dates, sim_params)
  }
  if (format == 'long') {
    sims = pivot_longer(
      sims, -Date,
      names_to = "variable",
      values_to = "value"
    )
  } else if (format != 'wide') {
    stop('format must be either long or wide')
  }
  sims
}


#' @rdname deprecated_simulation_functions
#' @export
simulation_fitted = function(model) {
  obsvars = unique(model$observed$data$var)
  simulation_condensed(model)[obsvars]
}
