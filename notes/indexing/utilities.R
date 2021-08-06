library(dplyr)

opt_list = 'beta0' # vector of parameters to optimize
params = McMasterPandemic::read_params('ICU1.csv')
params_to_opt = params[opt_list]
params_to_fix = params[!(names(params) %in% opt_list)]
state = McMasterPandemic::make_state(params = params)
M = McMasterPandemic::make_ratemat(state, params)
state_params = c(state, params)
pfun = McMasterPandemic:::pfun
do_step = McMasterPandemic::do_step

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
      var   = unlist(lapply(x, get_variables)),
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
    x$factors = left_join(x$factors, all_factors,
                          by = c('var', 'compl', 'invrs'))
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

update_ratemat2 = function(M, params, state, ratemat_struct) {
  for(i in seq_along(ratemat_struct)) {
    M[ratemat_struct[[i]]$ratemat_indices] = eval(
      ratemat_struct[[i]]$formula[[2]],
      as.list(c(state, params)))
  }
  M
}

do_step2 = function(state, M, params, ratemat_struct) {
  M = update_ratemat2(M, params, state, ratemat_struct)
  flows = M * state
  p_states = c("S", "E", "Ia", "Ip", "Im", "Is", "H",
               "H2", "ICUs", "ICUd", "D", "R")
  outflow = rowSums(flows[, p_states])
  inflow = colSums(flows)
  state - outflow + inflow
}
