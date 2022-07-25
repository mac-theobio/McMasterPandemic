##' Represent a Standard Unstructured Model as a flexmodel
##'
##' @inheritDotParams flexmodel
##' @family flexmodels
##' @family canned_models
##' @export
make_base_model <- function(...) {

  spec_check("0.0.5", "run_sim with TMB")

  model = flexmodel(...)
  if (spec_ver_gt('0.1.0')) {
    model$params = expand_params_S0(model$params, 1-1e-5)
  }
  model = (model

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
  )
  if(spec_ver_lt('0.1.1')) {
    # no deprecation period for add_parallel_accumulators
    model = add_parallel_accumulators(model, c("X", "N", "P", "V"))
  } else {
    model = (model
      %>% add_outflow(".+", "^(S|E|I|H|ICU|D|R)")

      # Update parameters for use with the linearized model
      %>% update_linearized_params('N', 1) # scale population to 1
      %>% update_linearized_params('E0', 1e-5) # perturbation

      # Set the disease-free equilibrium of the linearized model
      %>% update_disease_free_state('S', 'S0') # instead of N

      # Perturb the disease-free equilibrium of the linearized model
      %>% update_disease_free_state('E', 'E0')

      # Define outflows for the linearized model
      %>% add_linearized_outflow("^S$", "^S$")
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
    )
  }

  if (spec_ver_gt('0.1.2')) {
    model = (model
      %>% add_state_param_sum("Htotal", "^H2?$")
      %>% add_state_param_sum("ICU", "^ICU(s|d)$")
      %>% add_state_param_sum("Itotal", "^I(a|p|m|s)$")
      %>% add_sim_report_expr("Incidence", ~ (S_to_E) * (S))
      %>% add_lag_diff("^(X|D)$")
      %>% add_conv("^Incidence$")
      %>% update_condense_map(c(
        S = "S",
        E = "E",
        Itotal = "I",
        Htotal = 'H',
        ICU = 'ICU',
        R = 'R',
        lag_1_diff_X = 'hosp',
        X = 'X',
        lag_1_diff_D = 'death',
        D = 'D',
        S_to_E = 'foi',
        Incidence = 'incidence',
        conv_Incidence = 'report'
      ))
    )
  }

  model = update_tmb_indices(model)
  if (spec_ver_gt('0.1.0')) {
    model = update_initial_state(model, silent = TRUE)
  }
  model$classic_macpan_model = TRUE
  model
}


#' Make a Two-Dose Vaccination Model
#'
#' @inheritDotParams flexmodel
#' @param do_variant allow for different variants (TODO: more info obviously required here)
#' @param do_variant_mu different mu for different variants
#' @param do_wane including waning process where R boxes go back to S within vaccination strata
#' @param do_het make use of the zeta heterogeneity parameter in the force of infection
#' @family flexmodels
#' @family canned_models
#' @export
make_vaccination_model = function(..., do_variant = FALSE, do_variant_mu = FALSE, do_wane = FALSE, do_het = FALSE) {

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

  # ---------------------------
  # problem dimensions
  # ---------------------------
  (epi_states = c(attr(state, "epi_cat"))) # 14 base epidemiological categories
  (asymp_cat = c("S", "E", "Ia", "Ip", "R")) # 5 asymptomatic categories
  (vax_cat = c(attr(state, "vax_cat"))) # 5 vaccination categories/layers
  (accum = c("X", "V")) # two base parallel accumulator states
  (non_accum = base::setdiff(epi_states, accum)) # 12 base non-parallel accumulator states
  (non_accum_non_S = non_accum[-1]) # 11 base non-susceptible/non-accumulator states

  # dosing transitions across vaccination layers
  (dose_from = rep(asymp_cat, 2))
  (dose_to = c(asymp_cat, rep("V", 5)))

  # ---------------------------
  # Specify structure of the force of infection calculation
  # ---------------------------

  Istate = (c('Ia', 'Ip', 'Im', 'Is')
    %>% expand_names(vax_cat)
    %>% vec
  )
  if(do_het) {
    baseline_trans_rates =
      vec(
        'Ca',
        'Cp',
        '(1 - iso_m) * (Cm)',
        '(1 - iso_s) * (Cs)') *
      struc('(beta0) * (1/hetN)')
  } else {
    baseline_trans_rates =
      vec(
        'Ca',
        'Cp',
        '(1 - iso_m) * (Cm)',
        '(1 - iso_s) * (Cs)') *
      struc('(beta0) * (1/N)')
  }
  if(!do_variant) {
    vax_trans_red = struc_block(vec(
      '1',
      '1',
      '(1 - vax_efficacy_dose1)',
      '(1 - vax_efficacy_dose1)',
      '(1 - vax_efficacy_dose2)'),
      row_times = 1, col_times = 5)
  } else {
    vax_trans_red = struc_block(vec(
      '(1 - variant_prop) + (variant_advantage) * (variant_prop)',
      '(1 - variant_prop) + (variant_advantage) * (variant_prop)',
      '(1 - vax_efficacy_dose1) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose1) * (variant_advantage) * (variant_prop)',
      '(1 - vax_efficacy_dose1) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose1) * (variant_advantage) * (variant_prop)',
      '(1 - vax_efficacy_dose2) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose2) * (variant_advantage) * (variant_prop)'),
      row_times = 1, col_times = 5)
  }
  foi = (kronecker(vax_trans_red, t(baseline_trans_rates)) %*% Istate)
  if (do_het) {
    foi = vec("hetS" %_% vax_cat) * foi
  }

  alpha   = c("alpha", "alpha", "vax_alpha_dose1", "vax_alpha_dose1", "vax_alpha_dose2")


  if (!do_variant_mu) { # | getOption("MP_tmb_models_match_r")) {
    mu      = c("mu",    "mu",    "vax_mu_dose1",    "vax_mu_dose1",    "vax_mu_dose2")
  } else {
    ## variant-based mild-illness probability adjustment in vaccinated individuals
    mu      = c("mu",    "mu",
                "(1 - variant_prop) * (vax_mu_dose1) + (variant_prop) * (variant_vax_mu_dose1)",
                "(1 - variant_prop) * (vax_mu_dose1) + (variant_prop) * (variant_vax_mu_dose1)",
                "(1 - variant_prop) * (vax_mu_dose2) + (variant_prop) * (variant_vax_mu_dose2)")
  }

  sigma   = struc("sigma")
  gamma_p = struc("gamma_p")
  E_to_Ia_rates  = vec(           alpha ) * sigma
  E_to_Ip_rates  = vec(complement(alpha)) * sigma
  Ip_to_Im_rates = vec(              mu ) * gamma_p
  Ip_to_Is_rates = vec(complement(   mu)) * gamma_p

  model = flexmodel(...)
  if (spec_ver_gt('0.1.0')) {
    model$params = expand_params_S0(model$params, 1-1e-5)
  }

  # Sums across states
  model = (model
    %>% add_state_param_sum("asymp_unvax_N",       asymp_cat %_% "unvax")
    %>% add_state_param_sum("asymp_vaxprotect1_N", asymp_cat %_% "vaxprotect1")
    %>% add_state_param_sum("Stotal", "^S" %_% alt_group(vax_cat))
    %>% add_state_param_sum("Etotal", "^E" %_% alt_group(vax_cat))
    %>% add_state_param_sum("Itotal", "^I(a|s|p|m)" %_% alt_group(vax_cat))
    %>% add_state_param_sum("Htotal", "^H2?" %_% alt_group(vax_cat))
    %>% add_state_param_sum("ICU", "^ICU(s|d)" %_% alt_group(vax_cat))
    %>% add_state_param_sum("Rtotal", "^R" %_% alt_group(vax_cat))
    %>% add_state_param_sum("Xtotal", "^X" %_% alt_group(vax_cat))
    %>% add_state_param_sum("Dtotal", "^D" %_% alt_group(vax_cat))
  )

  if (do_het) {
    model = (model
      %>% update_params(one = 1)
      %>% add_factr("zeta_plus_1", ~ (zeta) + (one))
      %>% add_pow("hetN", "N", "zeta_plus_1", "one")
      %>% add_pow(
        "hetS" %_% vax_cat,
        "S" %_% vax_cat,
        "zeta",
        "one"
      )
    )
  }

  model = (model

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
      foi
      #t(baseline_trans_rates) %*% Imat %*% t(vax_trans_red)
    )

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
  )
  if (do_wane) {
    model = rep_rate(
      model,
      "R" %_% vax_cat,
      "S" %_% vax_cat,
      ~ (wane_rate)
    )
  }
  if (spec_ver_lt('0.1.1')) {
    # no deprecation period for add_parallel_accumulators
    model = add_parallel_accumulators(model, c('V' %_% vax_cat, 'X' %_% vax_cat))
  } else {
    model = (model

      # definitely need a better syntax here
      %>% add_outflow(
        ".+",
        "^" %+% alt_group(non_accum) %_% alt_group(vax_cat))

      # Update parameters for use with the linearized model
      # -- confirmed correct (TODO: check if params are missing? variant-related?)
      %>% update_linearized_params('^N$', 1) # scale population to 1
      %>% update_linearized_params('^E0$', 1e-5)
      %>% update_linearized_params('^vax_doses_per_day$', 0)
      %>% update_linearized_params('^vax_response_rate$', 0)
      %>% update_linearized_params('^vax_response_rate_R$', 0)
      # FIXED: vax_response_rate_R isn't necessary because it is regex-matched by vax_response_rate

      # Set the disease-free equilibrium of the linearized model
      %>% update_disease_free_state('S_unvax', 'S0')

      # Perturb the disease-free equilibrium of the linearized model
      %>% update_disease_free_state('E_unvax', 'E0')

      # Define outflows for the linearized model
      # -- confirmed that this is producing the correct indices
      %>% add_linearized_outflow("^S", "^S") # S_pos, S_pos
      %>% add_linearized_outflow(
        "^" %+% alt_group(non_accum_non_S) %_% alt_group(vax_cat), # notS_pos
        "^" %+% alt_group(non_accum)       %_% alt_group(vax_cat)) # p_states

      # Define state mappings used to put the initial state values in
      # the correct positions of the initial state vector
      %>% add_state_mappings(

        # regular expression to find states to drop before computing
        # the eigenvector of the linearized system
        # -- generated indices are correct
        eigen_drop_pattern = '^(X|V)',

        # regular expression to find states to drop from the eigenvector
        # before distributing individuals among infected compartments
        # -- generated indices are correct
        infected_drop_pattern = '^(S|D|R)',

        # regular expression to find states in the initial population
        # of susceptibles
        initial_susceptible_pattern = '^S_unvax$'
      )

      # Set the total number of individuals and the total number of
      # infected individuals in the initial state vector
      %>% initial_population(total = 'N', infected = 'E0')
    )
  }

  if (spec_ver_gt('0.1.2')) {
    foi_vec = vec("S" %_% vax_cat %_% "to" %_% "E" %_% vax_cat)
    S_vec = vec("S" %_% vax_cat)
    model = (model
      %>% add_sim_report_expr("Incidence_unvax", ~ (S_unvax_to_E_unvax) * (S_unvax))
      %>% add_sim_report_expr("Incidence_vaxdose1", ~ (S_vaxdose1_to_E_vaxdose1) * (S_vaxdose1))
      %>% add_sim_report_expr("Incidence_vaxprotect1", ~ (S_vaxprotect1_to_E_vaxprotect1) * (S_vaxprotect1))
      %>% add_sim_report_expr("Incidence_vaxdose2", ~ (S_vaxdose2_to_E_vaxdose2) * (S_vaxdose2))
      %>% add_sim_report_expr("Incidence_vaxprotect2", ~ (S_vaxprotect2_to_E_vaxprotect2) * (S_vaxprotect2))
      %>% add_sim_report_expr("Incidence", sum(foi_vec * S_vec))
      %>% add_conv("^Incidence")
      %>% add_lag_diff("^(X|D)total$")
      %>% update_condense_map(c(
        Stotal = "S",
        Etotal = "E",
        Itotal = "I",
        Htotal = 'H',
        ICU = 'ICU',
        Rtotal = "R",
        lag_1_diff_Xtotal = 'hosp',
        Xtotal = "X",
        lag_1_diff_Dtotal = 'death',
        Dtotal = "D",
        S_unvax_to_E_unvax = "foi_unvax",
        S_vaxdose1_to_E_vaxdose1 = "foi_vaxdose1",
        S_vaxprotect1_to_E_vaxprotect1 = "foi_vaxprotect1",
        S_vaxdose2_to_E_vaxdose2 = "foi_vaxdose2",
        S_vaxprotect2_to_E_vaxprotect2 = "foi_vaxprotect2",
        Incidence_unvax = "incidence_unvax",
        Incidence_vaxdose1 = "incidence_vaxdose1",
        Incidence_vaxprotect1 = "incidence_vaxprotect1",
        Incidence_vaxdose2 = "incidence_vaxdose2",
        Incidence_vaxprotect2 = "incidence_vaxprotect2",
        Incidence = 'incidence',
        conv_Incidence_unvax = 'report_unvax',
        conv_Incidence_vaxdose1 = 'report_vaxdose1',
        conv_Incidence_vaxprotect1 = 'report_vaxprotect1',
        conv_Incidence_vaxdose2 = 'report_vaxdose2',
        conv_Incidence_vaxprotect2 = 'report_vaxprotect2',
        conv_Incidence = 'report'
    ))
    )
  }
  model = update_tmb_indices(model)
  if (spec_ver_gt('0.1.0')) {
    model = update_initial_state(model, silent = TRUE)
  }
  model$classic_macpan_model = TRUE
  model
}

#' Converts a classic ageified model to a flexmodel
#'
#' @inheritDotParams flexmodel
#' @family flexmodels
#' @family canned_models
#' @export
make_ageified_model <- function(...,min_age = 0, max_age = 100, do_ageing = FALSE){
  args = list(...)
  McMasterPandemic::unpack(args)

  stopifnot(has_age(params))
  state = unclass(state)
  model = flexmodel(params = unlist(params), state = setNames(state, make.names(names(state))), start_date = start_date, end_date=end_date)

  epi_cat = c(attr(state, "epi_cat"))
  age_cat = c(attr(state, "age_cat"))
  age_cat = gsub("-", ".", age_cat)
  age_cat = gsub("\\+", ".", age_cat)

  beta0 = vec(paste0("beta0", 1:length(age_cat)))
  C = vec(c("(Ca)", "(Cp)", "(1 - iso_m) * (Cm)", "(1 - iso_s) * (Cs)"))
  N = struc_block(vec(inverse(paste0("N", 1:length(age_cat)))), row_times=1, col_times=4)

  pmat = as.struc(matrix(paste0('pmat', 1:(length(age_cat)*length(age_cat))), length(age_cat), length(age_cat)))

  Istate = c('Ia', 'Ip', 'Im', 'Is')
  Ivec = expand_names(Istate, age_cat)
  Imat = as.struc(matrix(Ivec, ncol = 4, byrow=TRUE))
  foi_vec = rowSums((beta0 %*% t(C)) * (pmat %*% (Imat * N)))

  model = (model
           # Flow within age categories,
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
           %>% rep_rate("E", "Ia",      ~(           alpha ) * (sigma))
           %>% rep_rate("E", "Ip",      ~(1-alpha) * (sigma))
           %>% rep_rate("Ip", "Im",     ~(              mu ) * (gamma_p))
           %>% rep_rate("Ip", "Is",     ~(1-mu) * (gamma_p))

           #Add variable rates within a single age category
           %>% vec_rate("S", "E", foi_vec)

           #Remove outflow to accumulators
           %>% add_outflow(".+", "^(S|E|I|H|ICU|D|R)")
  )

  #Add rates between age categories (if requested)
  if(do_ageing){
    ageing_rate = length(age_cat)/(365*(max_age-min_age))
    model = update_params(model, c(ageing_rate = ageing_rate))
    ageing_cat = epi_cat[epi_cat!="D"]
    for(i in ageing_cat){
      ageified_states = paste(i, age_cat, sep="_")
      states_from = ageified_states[1:(length(ageified_states)-1)]
      states_to = ageified_states[2:length(ageified_states)]
      model = model%>%rep_rate(states_from, states_to, ~(ageing_rate))
    }
  }


  model$classic_macpan_model = TRUE
  return(model)
}

#' Converts a classic testified model to a flexmodel
#'
#' @inheritDotParams flexmodel
#' @family flexmodels
#' @family canned_models
#' @export
make_testified_model <-function(...){
  args = list(...)
  McMasterPandemic::unpack(args)

  stopifnot(McMasterPandemic:::has_testing(params=params))

  state = make_state(params=params)
  state_names = names(state)
  state_names[state_names == "N"]="Neg"
  state_names[state_names == "P"]="Pos"
  state = setNames(state, state_names)
  state = unclass(state)

  model = flexmodel(params=unlist(params), state = state, start_date=start_date, end_date=end_date)


  epi_cat = c(attr(state, "epi_cat"))
  expand_set = McMasterPandemic:::exclude_states(epi_cat, non_expanded_states)

  posvec = paste("P", expand_set, sep="_")
  negvec = vector(mode="character", length = length(expand_set))
  for(i in 1:length(negvec)){
    negvec[i]=complement(posvec[i])
  }
  posvec = vec(posvec)
  negvec = vec(negvec)

  wtsvec = vector(mode = "character", length = length(expand_set))
  for(i in 1:length(expand_set)){
    if(any(expand_set[i]==asymp_cat)) wtsvec[i] = "W_asymp"
    else wtsvec[i] = "1"
  }
  wtsvec = vec(wtsvec)

  testing_intensity = vec(rep("testing_intensity", length(posvec)))
  u_to_n_flow = testing_intensity*wtsvec*negvec
  u_to_p_flow = testing_intensity*wtsvec*posvec

  Ia_set = paste("Ia", test_extensions, sep="_")
  Im_set = paste("Im", test_extensions, sep="_")
  Is_set = paste("Is", test_extensions, sep="_")
  Is_set_rd = paste("Is", test_extensions[test_extensions!="u"], sep="_")
  Ip_set = paste("Ip", test_extensions, sep="_")
  ICUs_set = paste("ICUs", test_extensions, sep="_")
  ICUs_set_rd = paste("ICUs", test_extensions[test_extensions!="u"], sep="_")
  ICUd_set = paste("ICUd", test_extensions, sep="_")
  ICUd_set_rd = paste("ICUd", test_extensions[test_extensions!="u"], sep="_")
  H2_set = paste("H2", test_extensions, sep="_")
  H_set = paste("H", test_extensions, sep="_")
  H_set_rd = paste("H", test_extensions[test_extensions!="u"], sep="_")
  R_set = paste("R", test_extensions, sep="_")
  E_set = paste("E", test_extensions, sep="_")
  S_set = paste("S", test_extensions, sep="_")

  u_set = paste(expand_set, "u", sep="_")
  p_set = paste(expand_set, "p", sep="_")
  p_set_rd = paste(c("H", "ICUs", "ICUd"), "p", sep="_")
  n_set = paste(expand_set, "n", sep="_")
  n_set_rd = paste(c("H", "ICUs", "ICUd"), "n", sep="_")
  t_set = paste(expand_set, "t", sep="_")

  Is_u_outflow = vec("(1 - nonhosp_mort) * (gamma_s) * (    phi1)", "(1 - nonhosp_mort) * (gamma_s) * (1 - phi1) * (1 - phi2)", "(1 - nonhosp_mort) * (gamma_s) * (1 - phi1) * (    phi2)")
  Is_u_outflow_p = Is_u_outflow * vec(rep("P_Is", 3))
  Is_u_outflow_n = Is_u_outflow * vec(rep("1-P_Is", 3))


  Istate = vec(c(Is_set, Ip_set, Im_set, Is_set))
  trans_rates =
    vec(
      rep('Ca', 4),
      rep('Cp', 4),
      rep('(1 - iso_m) * (Cm)', 4),
      rep('(1 - iso_s) * (Cs)', 4)) *
    struc('(beta0) * (1/N)')
  trans_rate = t(trans_rates)%*%Istate

  model = (model
           # Flow from expanded states to expanded states with constant rate
           %>% rep_rate(Ia_set, R_set, ~                      (gamma_a))
           %>% rep_rate(Im_set, R_set, ~                      (gamma_m))
           %>% rep_rate(Is_set_rd, H_set_rd, ~ (1 - nonhosp_mort) * (gamma_s) * (    phi1))
           %>% rep_rate(Is_set_rd, ICUs_set_rd, ~ (1 - nonhosp_mort) * (gamma_s) * (1 - phi1) * (1 - phi2))
           %>% rep_rate(Is_set_rd, ICUd_set_rd, ~ (1 - nonhosp_mort) * (gamma_s) * (1 - phi1) * (    phi2))
           %>% rep_rate(ICUs_set, H2_set, ~                                  (    psi1))
           %>% rep_rate(H2_set, R_set,    ~                                  (    psi3))
           %>% rep_rate(H_set,  R_set,    ~ (rho))
           %>% rep_rate(E_set, Ia_set,      ~(           alpha ) * (sigma))
           %>% rep_rate(E_set, Ip_set,      ~(1-alpha) * (sigma))
           %>% rep_rate(Ip_set, Im_set,     ~(              mu ) * (gamma_p))
           %>% rep_rate(Ip_set, Is_set,     ~(1-mu) * (gamma_p))

           # Flow from expanded states to non-expanded states with constant rate
           %>% rep_rate(Is_set,   "D",    ~ (    nonhosp_mort) * (gamma_s))
           %>% rep_rate(Is_set,   "X",    ~ (1 - nonhosp_mort) * (gamma_s) * (    phi1))
           %>% rep_rate(ICUd_set, "D",    ~                                  (    psi2))

           # Force of infection flow
           %>% rep_rate(S_set, E_set, as.character(trans_rate))

           # Flow related to testing
           %>% vec_rate(u_set, n_set, u_to_n_flow)
           %>% vec_rate(u_set, p_set, u_to_p_flow)
           %>% rep_rate(n_set, u_set, ~ omega)
           %>% rep_rate(p_set, t_set, ~ omega)

           # Special flows related to hospital visits
           %>% vec_rate(rep("Is_u", 3), p_set_rd, Is_u_outflow_p)
           %>% vec_rate(rep("Is_u", 3), n_set_rd, Is_u_outflow_n)

           # Pos and Neg Accumulators
           %>% vec_rate(u_set, rep("Pos", length(u_set)), u_to_p_flow)
           %>% vec_rate(u_set, rep("Neg", length(u_set)), u_to_n_flow)
           %>% add_outflow(".+", "^(S|E|I|H|ICU|D|R)")
  )

  model$classic_macpan_model = TRUE
  return(model)
}

##' SIR Model
##'
##' @inheritDotParams flexmodel
##' @family flexmodels
##' @family canned_models
##' @export
make_sir_model = function(...) {
  l... = list(...)
  stopifnot(all(names(l...$state) %in% c("S", "I", "R")))
  stopifnot(all(c("S", "I", "R") %in% names(l...$state)))
  stopifnot(all(names(l...$params) %in% c("beta", "gamma", "N", "nb_disp_S", "nb_disp_I", "nb_disp_R")))
  stopifnot(all(c("beta", "gamma") %in% names(l...$params)))
  model = flexmodel(...)
  model$do_make_state = FALSE
  if("N" %in% names(l...$params)) {
    stopifnot(l...$params[["N"]] == sum(l...$state))
  } else {
    model = add_state_param_sum(model, "N", alt_group(names(l...$state), exact = TRUE))
  }
  model = (model
    %>% add_rate("S", "I", ~ (beta) * (I) * (1/N))
    %>% add_rate("I", "R", ~ (gamma))
    %>% add_outflow
    %>% update_condense_map
    %>% update_tmb_indices
  )
  model
}

#' Hello World
#'
#' Return a hello-world SIR model
#'
#' @seealso \code{\link{make_sir_model}}
#' @family flexmodels
#' @family canned_models
#' @export
make_hello_world_model = function() {
  state = c(S = 20000, I = 100, R = 0)
  (
    flexmodel(
      params = c(
        gamma = 0.06,
        beta = 0.15,
        N = sum(state)
      ),
      state = state,
      start_date = "2016-07-08",
      end_date = "2016-12-31",
      do_hazard = FALSE,
      do_make_state = FALSE
    )
    %>% add_rate("S", "I", ~ (1/N) * (beta) * (I))
    %>% add_rate("I", "R", ~ (gamma))
  )
}
