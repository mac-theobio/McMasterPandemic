#' Configure Epidemiological Summaries
#'
#' @param model \code{\link{flexmodel}} object
#' @param exposed_state_nm name of the exposed state
#' @param foi_nm name of the force of infection associated with
#' the exposed state
#' @param prop_susceptible_nm name of the variable that contains the
#' user-defined effective proportion of individuals who are susceptible
#' @param trans_rate_nm name of the parameter giving the baseline transmission
#' rate
#' @param ... named parameter values to update for the cohort
#' model -- a common example is \code{N = 1}, which allows
#' certain models to track relative population sizes
#'
#' @export
configure_epi_summaries = function(
  model,
  exposed_state_nm,
  foi_nm,
  prop_susceptible_nm,
  trans_rate_nm,
  ...
) {

  assert_len1_char(exposed_state_nm)
  assert_len1_char(foi_nm)
  assert_len1_char(prop_susceptible_nm)
  assert_len1_char(trans_rate_nm)

  model$summary_config = list(
    var_nms = nlist(
      exposed_state_nm,
      foi_nm,
      prop_susceptible_nm,
      trans_rate_nm
    ),
    kernel_param_updates = list(...)
  )

  model
}

# Epidemic Cohort Model and Kernel Methods
#
#
#
# @param model \code{\link{flexmodel}} object
# @param exposed_state_nm name of the exposed state
# @param foi_nm name of the force of infection associated with
# the exposed state
# @param ... named parameter values to update for the cohort
# model -- a common example is \code{N = 1}, which allows
# certain models to track relative population sizes
#
# @return \code{\link{flexmodel}} object that represents the
# cohort model associated with \code{model}


#' @describeIn configure_epi_summaries obtain the kernel from \code{model}
#' @export
epi_kernel = function(model) {
  cohort_sim = (model
    %>% as.flexmodel_cohort
    #%>% epi_cohort(exposed_state_nm, foi_nm, ...)
    %>% simulation_history(condense = TRUE, add_dates = FALSE)
  )

  tol = getOption("MP_kernel_tol")
  final_number_exposed = tail(cohort_sim$exposed, 1L)

  if (final_number_exposed > tol) {
    warning(
      "\nfinal fraction of exposed individuals, ",
      final_number_exposed,
      ",\nis not below tolerance, ",
      tol,
      ".\nyou might want to adjust the tolerance with",
      "\noptions(MP_kernel_tol = xx)"
    )
  }

  cohort_sim$kernel
}

# Kernel Summaries
#
# These functions summarize the epidemiological kernel

epi_R0 = function(kern) sum(kern)
epi_Gbar = function(kern) sum(seq_along(kern) * kern) / epi_R0(kern)
epi_r = function(kern) {
  ub = getOption("MP_kernel_r_upper_bound")
  uniroot(make_lotka_euler(kern), c(0, ub))$root
}

#' @describeIn configure_epi_summaries compute epidemic growth parameters for \code{model}
#' @export
epi_pars = function(model) {
  k = epi_kernel(model)
  c(
    R0 = epi_R0(k),
    Gbar = epi_Gbar(k),
    r = epi_r(k)
  )
}


  # kern = epi_kernel(model, exposed_state_nm, foi_nm, ...)
  # R0 = epi_R0(kern)
  # rel_beta = pars_time_norm(model)[[trans_rate_nm]]
  # sim_cols = c("Date", prop_susceptible_nm)
  # sims = setNames(
  #   simulation_history(model, condense = TRUE)[, sim_cols],
  #   c("Date", "prop_susceptibles")
  # )
  # sims$rel_trans_rate = rel_beta
  # sims$Rt = R0 * rel_beta * sims[[prop_susceptible_nm]]
  # sims

