##' Represent a Standard Unstructured Model as a flexmodel
##'
##' @inheritDotParams init_model
##' @family flexmodels
##' @family canned_models
##' @export
make_base_model <- function(...) {

  spec_check("0.0.5", "run_sim with TMB")

  model <- (init_model(...)
            %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
            %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
            %>% add_rate("Ia", "R", ~ (gamma_a))
            %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
            %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
            %>% add_rate("Im", "R", ~ (gamma_m))
            %>% add_rate("Is", "H", ~
                           (1 - nonhosp_mort) * (phi1) * (gamma_s))
            %>% add_rate("Is", "ICUs", ~
                           (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
            %>% add_rate("Is", "ICUd", ~
                           (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
            %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
            %>% add_rate("ICUs", "H2", ~ (psi1))
            %>% add_rate("ICUd", "D", ~ (psi2))
            %>% add_rate("H2", "R", ~ (psi3))
            %>% add_rate("H", "R", ~ (rho))
            %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
            %>% add_rate("S", "E", ~
                           (Ia) * (beta0) * (1 / N) * (Ca) +
                           (Ip) * (beta0) * (1 / N) * (Cp) +
                           (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
                           (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))

            # TODO: choice of add_parallel_accumulator vs add_outflow will
            # need to depend on the spec version
            %>% add_parallel_accumulators(c("X", "N", "P", "V"))
            %>% add_outflow(".+", "^(S|E|I|H|ICU|D|R)")

            # Update parameters for use with the linearized model
            %>% update_linearized_params('N', 1) # scale population to 1
            %>% update_linearized_params('E0', 1e-5)

            # Set the disease-free equilibrium of the linearized model
            %>% update_disease_free_state('S', 'N')

            # Perturb the disease-free equilibrium of the linearized model
            %>% update_disease_free_state('E', 'E0')

            # Define outflows for the linearized model
            %>% add_linearized_outflow("S", "S")
            %>% add_linearized_outflow("^(E|I|H|ICU|D|R)", "^(S|E|I|H|ICU|D|R)")

            # Define state mappings used to put the initial state values in
            # the correct positions of the initial state vector
            %>% add_state_mappings(

              # regular expression to find states to drop before computing
              # the eigenvector of the linearized system
              eigen_drop_pattern = '^(X|V)',

              # regular expression to find states to drop from the eigenvector
              # before distributing individuals among infected compartments
              infected_drop_pattern = '^(S|D|R)',

              # regular expression to find states in the initial population
              # of susceptibles
              initial_susceptible_pattern = '^S$'
            )

            # Set the total number of individuals and the total number of
            # infected individuals in the initial state vector
            %>% initial_population(total = 'N', infected = 'E0')

            %>% add_tmb_indices()
  )
  return(model)
}


#' Make a Two-Dose Vaccination Model
#'
#' @inheritDotParams init_model
#' @family flexmodels
#' @family canned_models
#' @export
make_vaccination_model = function(...) {

  spec_check("0.1.0", "model structure")

  args = list(...)
  unpack(args)

  stopifnot(has_vax(params))
  stopifnot(attr(params, "model_type") != 'twodose')
  if(is.null(state)) state = make_state(params = params)
  if(has_vax(state)) {
    stopifnot(attr(state, "model_type") != 'twodose')
  } else {
    state <- expand_state_vax(
      x = state,
      model_type = "twodose",
      unif = FALSE
    )
  }

  # problem dimensions
  (epi_states = c(attr(state, "epi_cat")))
  (asymp_cat = c("S", "E", "Ia", "Ip", "R"))
  (vax_cat = c(attr(state, "vax_cat")))
  (dose_from = rep(asymp_cat, 2))
  (dose_to = c(asymp_cat, rep("V", 5)))

  # Specify structure of the force of infection calculation
  Istate = (c('Ia', 'Ip', 'Im', 'Is')
    %>% expand_names(vax_cat)
    %>% vec
  )
  baseline_trans_rates =
    vec(
      'Ca',
      'Cp',
      '(1 - iso_m) * (Cm)',
      '(1 - iso_s) * (Cs)') *
    struc('(beta0) * (1/N)')
  vax_trans_red = struc_block(vec(
    '1',
    '1',
    '(1 - vax_efficacy_dose1)',
    '(1 - vax_efficacy_dose1)',
    '(1 - vax_efficacy_dose2)'),
    row_times = 1, col_times = 5)

  alpha   = c("alpha", "alpha", "vax_alpha_dose1", "vax_alpha_dose1", "vax_alpha_dose2")
  mu      = c("mu",    "mu",    "vax_mu_dose1",    "vax_mu_dose1",    "vax_mu_dose2")
  sigma   = struc("sigma")
  gamma_p = struc("gamma_p")
  E_to_Ia_rates  = vec(           alpha ) * sigma
  E_to_Ip_rates  = vec(complement(alpha)) * sigma
  Ip_to_Im_rates = vec(              mu ) * gamma_p
  Ip_to_Is_rates = vec(complement(   mu)) * gamma_p


  (init_model(...)

    # Flow within vaccination categories,
    # with constant rates across categories
    %>% rep_rate("Ia",   "R",    ~                      (gamma_a))
    %>% rep_rate("Im",   "R",    ~                      (gamma_m))
    %>% rep_rate("Is",   "D",    ~ (    nonhosp_mort) * (gamma_s))
    %>% rep_rate("Is",   "H",    ~ (1 - nonhosp_mort) * (gamma_s) * (    phi1))
    %>% rep_rate("Is",   "X",    ~ (1 - nonhosp_mort) * (gamma_s) * (    phi1))
    %>% rep_rate("Is",   "ICUs", ~ (1 - nonhosp_mort) * (gamma_s) * (1 - phi1) * (1 - phi2))
    %>% rep_rate("Is",   "ICUd", ~ (1 - nonhosp_mort) * (gamma_s) * (1 - phi1) * (    phi2))
    %>% rep_rate("ICUs", "H2",   ~                                  (    psi1))
    %>% rep_rate("ICUd", "D",    ~                                  (    psi2))
    %>% rep_rate("H2",   "R",    ~                                  (    psi3))
    %>% rep_rate("H",    "R",    ~ (rho))

    # Flow within vaccination categories,
    # with rates that depend on category
    # (see struc objects created above)
    %>% vec_rate("E", "Ia",  E_to_Ia_rates)
    %>% vec_rate("E", "Ip",  E_to_Ip_rates)
    %>% vec_rate("Ip", "Im", Ip_to_Im_rates)
    %>% vec_rate("Ip", "Is", Ip_to_Is_rates)

    # Vaccination Response Rates
    %>% add_rate("R_vaxdose1", "R_vaxprotect1",  ~ (vax_response_rate_R))
    %>% add_rate("R_vaxdose2", "R_vaxprotect2",  ~ (vax_response_rate_R))
    %>% add_rate("S_vaxdose1", "S_vaxprotect1",  ~ (vax_response_rate))
    %>% add_rate("S_vaxdose2", "S_vaxprotect2",  ~ (vax_response_rate))

    # Forces of Infection
    %>% vec_rate(
      "S" %_% vax_cat,
      "E" %_% vax_cat,
      kronecker(vax_trans_red, t(baseline_trans_rates)) %*% Istate
    )

    # Sums across vaccination categories
    %>% add_state_param_sum("asymp_unvax_N",       asymp_cat %_% "unvax")
    %>% add_state_param_sum("asymp_vaxprotect1_N", asymp_cat %_% "vaxprotect1")

    # Flow among vaccination categories
    # (see dose_* above for epi states that are involved)
    %>% rep_rate(
      dose_from %_% 'unvax',
      dose_to   %_% 'vaxdose1',
      ~ (    vax_prop_first_dose) * (vax_doses_per_day) * (1 / asymp_unvax_N))
    %>% rep_rate(
      dose_from %_% 'vaxprotect1',
      dose_to   %_% 'vaxdose2',
      ~ (1 - vax_prop_first_dose) * (vax_doses_per_day) * (1 / asymp_vaxprotect1_N))

    # Technical steps
    %>% add_parallel_accumulators(c('V' %_% vax_cat, 'X' %_% vax_cat))
    %>% add_tmb_indices()
  )
}
