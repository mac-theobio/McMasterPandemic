#' Initialize Compartmental Model
#'
#' TODO: `recompute` argument for functions below
#' TODO: `overwrite` argument for add functions that force `recompute`
#' TODO: create versions of `rate/add_rate` functions that allow structure.
#'       maybe `rates/add_rates`??
#'
#' @param params a \code{param_pansim} object
#' @param state a \code{state_pansim} object
#' @return object representing a compartmental model
#' @export
init_model <- function(params, state) {
  model = list(
    state = state,
    params = params,
    ratemat = make_ratemat(state, params),
    rates = list(),
    tmb_indices = list()  # TODO: clarify index structure here once we converge
  )

  if(spec_version_greater_than('0.0.1')) {
    model$parallel_accumulators = character()
  }

  return(model)
}

spec_version = function() {
  # https://canmod.net/misc/flex_specs
  parse_version(getOption('MP_flex_spec_version'))
}

spec_check = function(introduced_version, msg_if_version_too_old) {
  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")
  if(parse_version(current_version) < parse_version(introduced_version)) {
  stop("\n\n", msg_if_version_too_old, "\n",
       "The specification version currently being used is ",
       getOption('MP_flex_spec_version'), '\n',
       "See ", getOption('MP_flex_spec_doc_site'),
       " for more information on specification versions.",
       call. = FALSE)
  }
}

spec_version_equal_to = function(version) {
  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")
  parse_version(current_version) == parse_version(version)
}

spec_version_greater_than = function(version) {
  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")
  parse_version(current_version) > parse_version(version)
}

spec_version_less_than = function(version) {
  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")
  parse_version(current_version) < parse_version(version)
}

spec_version_between = function(version_left, version_right) {
  # https://canmod.net/misc/flex_specs
  current_version = getOption("MP_flex_spec_version")
  spec_version_greater_than(version_left) &
    spec_version_less_than(version_right)
}

#' Rate Structure
#'
#' Define how the rate of flow from one compartment to another
#' depends on the parameters and state variables.
#'
#' @param model compartmental model
#' @param from Name of state where flow is happening from
#' @param to Name of state where flow is happening to
#' @param formula Model formula defining dependence of the rate on
#' parameters and state variables
#' @return another compartmental model with an additional non-zero rate matrix
#' element specified
#' @export
add_rate = function(model, from, to, formula) {
  added_rate = (
    rate(from, to, formula, model$state, model$params, model$ratemat)
    %>% list
    %>% setNames(paste(from, to, sep = '_to_'))
  )
  model$rates = c(model$rates, added_rate)
  return(model)
}

#' @export
rate = function(from, to, formula, state, params, ratemat) {
  # TODO: test for formula structure
  # TODO: test that from and to are available in the state vector
  M = ratemat
  stopifnot(is.character(from),
            is.character(to),
            length(from) == 1L,
            length(to) == 1L,
            inherits(formula, 'formula'))

  # regex pattern for finding variables
  # (e.g. any parameter or state variable)
  # variable_regex looks like this '(beta0|Ca|...|zeta|S|E|Ia|...|V)'
  variable_regex = function(...) {
    character_class =
      (list(...)
       %>% lapply(names)
       %>% unlist
       %>% paste0(collapse = "|")
      )
    paste0('(', character_class, ')', sep = '')
  }
  get_variables = function(x) {
    r = regexpr(variable_regex(params, state), x)
    regmatches(x, r)
  }
  # this only works because complements (1 - x) and inverses (1 / x) are
  # so similar in structure
  find_operators = function(x, operator) {
    grepl(
      paste0(
        '\\( *1 *',
        operator,
        ' *',
        variable_regex(params, state),
        collapse = ''),
      x)
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
                 %>% lapply(factor_table) %>% bind_rows(.id = 'prod_indx')
    )
    x$ratemat_indices =
      do.call(McMasterPandemic:::pfun, c(x[c('from', 'to')], list(mat = M)))
    x$factors$var_indx = (
      x$factors$var
      %>% outer(names(c(state, params)), '==')
      %>% apply(1, which)
    )
    if(spec_version_greater_than('0.0.1')) {
      x$state_dependent = any(x$factors$var_indx <= length(state))
    }
    x
  }
  structure(
    product_list(list(from = from, to = to, formula = formula)),
    class = 'rate-struct')
}

#' Add Parallel Accumulators
#'
#' Add parallel accumulators to a compartmental model.
#'
#' @param model TODO
#' @param state_patterns regular expressions for identifying states as
#' parallel accumulators
#' @return another compartmental model with parallel accumulators specified
#' @export
add_parallel_accumulators = function(model, state_patterns) {
  model$parallel_accumulators = parallel_accumulators(model, state_patterns)
  return(model)
}

#' @export
parallel_accumulators = function(model, state_patterns) {
  spec_check(
    introduced_version = '0.0.2',
    msg_if_version_too_old =
      "Parallel accumulators are not introduced until spec version 0.0.2.")
  (
    state_patterns
    %>% lapply(function(x) {grep(x, colnames(model$ratemat), value = TRUE)})
    %>% unlist
  )
}

#' Add TMB Indices
#'
#' Add, to a compartmental model, indices used to access appropriate values
#' during simulation and calibration using TMB
#'
#' @param model compartmental model
#' @param another compartmental model with indices for TMB
#' @export
add_tmb_indices = function(model) {
  model$tmb_indices = tmb_indices(model)
  return(model)
}

#' @export
tmb_indices = function(model) {
  state_dependent_rates = function(model) {
    model$rates[sapply(model$rates, '[[', 'state_dependent')]
  }

  ratemat_indices = function(rates, state_params) {
    sp = state_params
    ratemat_indices = sapply(rates, `[[`, 'ratemat_indices')
    spi = {lapply(rates, function(y) {y$factors$var_indx}) %>% unlist}
    count = sapply(rates, function(y) {
      nrow(y$factors)
    })
    modifier = lapply(unname(rates), '[[', 'factors') %>%
      bind_rows(.id = 'rate_indx') %>%
      mutate(add = as.logical(c(0, diff(as.numeric(prod_indx))))) %>%
      mutate(modifier = 4 * add + 2 * invrs + compl) %>%
      `$`("modifier")
    names(spi) = colnames(ratemat_indices) = names(count) = NULL
    list(
      from = ratemat_indices[1,],
      to = ratemat_indices[2,],
      count = count,
      spi = spi,
      modifier = modifier
    )
  }

  sp = c(model$state, model$params)
  indices = list(make_ratemat_indices = ratemat_indices(model$rates, sp))

  if(spec_version_greater_than('0.0.1')) {
    indices$update_ratemat_indices =
      ratemat_indices(state_dependent_rates(model), sp)
    indices$par_accum_indices =
      which(colnames(model$ratemat) %in% model$parallel_accumulators)
  }
  return(indices)
}
