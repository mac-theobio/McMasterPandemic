# Example Rate Matrix Structure
ratemat_struc = list(
  # recovery
  list(from = "Ia", to = "R",
       formula = ~ (gamma_a)),
  # hospitalizations
  list(from = "Is", to = "ICUs",
       formula = ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s)),
  # force of infection
  list(from = "S",  to = "E",
       formula = ~
         (Ia) * (beta0) * (1/N) * (Ca) +
         (Ip) * (beta0) * (1/N) * (Cp) +
         (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
         (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
)

library(dplyr)

opt_list = 'beta0' # vector of parameters to optimize
params = McMasterPandemic::read_params('ICU1.csv')
params_to_opt = params[opt_list]
params_to_fix = params[!(names(params) %in% opt_list)]
state = McMasterPandemic::make_state(params = params)
M = McMasterPandemic::make_ratemat(state, params)
state_params = c(state, params)
pfun = McMasterPandemic:::pfun

#' Rate Structure
#'
#' Define how the rate of flow from one compartment to another
#' depends on the parameters and state variables.
#'
#' @param from Name of state where flow is happening from
#' @param to Name of state where flow is happening to
#' @param formula Model formula defining dependence of the rate on
#' parameters and state variables
#' @return List with the information in the arguments.
rate = function(from, to, formula) {
  stopifnot(is.character(from),
            is.character(to),
            length(from) == 1L,
            length(to) == 1L,
            inherits(formula, 'formula'))

  # regex pattern for finding any variable
  # (i.e. any parameter or state variable)
  # variable_regex looks like this '(beta0|Ca|...|zeta|S|E|Ia|...|V)'
  variable_regex = paste0(
    '(',
    paste0(c(names(params), names(state)), collapse = '|'),
    ')', sep = '')
  get_variables = function(x) {
    r = regexpr(variable_regex, x)
    regmatches(x, r)
  }
  # this only works because complements (1 - x) and inverses (1 / x) are
  # so similar in structure
  find_operators = function(x, operator) {
    grepl(paste0('\\( *1 *', operator, ' *', variable_regex, collapse = ''), x)
  }
  factor_table = function(x) {
    data.frame(
      var = unlist(lapply(x, get_variables)),
      compl = unlist(lapply(x, find_operators, '-')),
      invrs = unlist(lapply(x, find_operators, '/')))
  }
  product_list = function(x) {
    # x is a list defining rate matrix structure. return an altered x with
    # factor list appended to the structure of each non-zero rate matrix
    # element
    x$factors = (x$formula
                 %>% as.character %>% getElement(2L)
                 %>% strsplit(split = '\\+') %>% getElement(1L)
                 %>% strsplit(split = '\\*')
                 %>% lapply(factor_table) %>% bind_rows(.id = 'product_index')
                 #%>% mutate(value = c(state, params)[var])
                 #%>% mutate(factor = ifelse(compl, 1 - value, ifelse(invrs, 1 / value, value)))
    )
    x$ratemat_indices =
      do.call(pfun, c(x[c('from', 'to')], list(mat = M)))
    x
  }

  # TODO: test for formula structure
  # TODO: test that from and to are available in the state vector
  structure(
    product_list(list(from = from, to = to, formula = formula)),
    class = 'rate-struct')
}

#' Rate Matrix Structure
#'
#' @param ... objects of class rate-struct, created by the rate function
#' @return list of rate-struct objects
mk_ratemat_struct = function(...) {

  mk_indices = function(x) {
    apply(outer(x, names(state_params), '=='), 1, which)
  }
  update_w_indices = function(x) {
    # line up the indices into the factor vector, with the lists of
    # factors for each product
    x$factors = left_join(x$factors, all_factors, by = c('var', 'compl', 'invrs'))
    x
  }

  struct = list(...)
  lapply(struct, inherits, 'rate-struct') %>% lapply(stopifnot)
  all_factors = (struct
    %>% lapply(`[[`, "factors")
    %>% bind_rows
    %>% select(var, compl, invrs)
    %>% distinct
    %>% mutate(state_param_index = mk_indices(var))
    %>% mutate(factor_index = seq_along(var))
  )
  structure(
    lapply(struct, update_w_indices),
    class = 'ratemat-struct',
    all_factors = all_factors)
}

`ratemat_env<-` = function(struct, value) {
  set_formula_env = function(x) {
    environment(x$formula) = value
    x
  }
  lapply(struct, set_formula_env)
}

save_ratemat_envs = function(struct) {
  save_formula_env = function(x) {
    x$ratemat_env = environment(x$formula)
    x
  }
  lapply(struct, save_formula_env)
}

restore_ratemat_envs = function(struct) {
  restore_formula_env = function(x) {
    if((!is.null(x$ratemat_env)) & (is.environment(x$ratemat_env))) {
      environment(x$formula) = x$ratemat_env
    } else {
      warning("no saved environment available to restore")
    }
    x
  }
}

update_ratemat2 = function(ratemat_struct, model_env = .GlobalEnv) {
  ratemat_env(ratemat_struct) = model_env
  update_rate = function()

  restore
}

ratemat_struct = mk_ratemat_struct(
  rate("E", "Ia", ~ (alpha) * (sigma)),
  rate("E", "Ip", ~ (1 - alpha) * (sigma)),
  rate("Ia", "R", ~ (gamma_a)),
  rate("Ip", "Im", ~ (mu) * (gamma_p)),
  rate("Ip", "Is", ~ (1 - mu) * (gamma_p)),
  rate("Im", "R", ~ (gamma_m)),
  rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s)),
  rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s)),
  rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s)),
  rate("Is", "D", ~ (nonhosp_mort) * (gamma_s)),
  rate("ICUs", "H2", ~ (psi1)), ## ICU to post-ICU acute care
  rate("ICUd", "D", ~ (psi2)), ## ICU to death
  rate("H2", "R", ~ (psi3)), ## post-ICU to discharge
  ## H now means 'acute care' only; all H survive & are discharged
  #list(from = "H", to = "D", formula = ~ 0),
  rate("H", "R", ~ (rho)), ## all acute-care survive
  rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s)), ## assuming that hosp admissions mean *all* (acute-care + ICU)
  # force of infection
  rate("S",  "E", ~
         (Ia) * (beta0) * (1/N) * (Ca) +
         (Ip) * (beta0) * (1/N) * (Cp) +
         (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
         (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
)



