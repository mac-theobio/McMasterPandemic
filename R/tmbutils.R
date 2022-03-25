# constructing names and strings ----------------------

#' Paste with Underscore Separator
#' @export
`%_%` = function(x, y) paste(x, y, sep = "_")

#' Paste with Blank Separator
#'
#' Like Python string `+`
#' @export
`%+%` = function(x, y) paste(x, y, sep = "")

wrap_exact = function(x) {
  "^" %+% x %+% "$"
}

##' Regex Alternation Group
##'
##' @export
alt_group = function(x, exact = FALSE, negate = FALSE) {
  x = "(" %+% paste0(x, collapse = "|") %+% ")"
  if (negate) {
    x = "(?!(?:" %+% x %+% ")$).*"
  }
  if (exact) {
    x = "^" %+% x %+% "$"
  }
  x
}

names_or_values = function(x) {
  if (!is.character(x)) {
    x = names(x)
  }
  x
}

##' @export
all_in = names_or_values

##' @export
any_var = function(x) {
  (x
   %>% names_or_values
   %>% alt_group(exact = TRUE)
  )
}

##' @export
all_except = function(x) {
  (x
   %>% names_or_values
   %>% alt_group(exact = TRUE, negate = TRUE)
  )
}

# constructing vectors ---------------------

#' @export
layered_zero_state = function(...) {
  state_nms = (list(...)
   %>% lapply(as.character)
   %>% Reduce(f = expand_names)
  )
  setNames(rep(0, length(state_nms)), state_nms)
}

# null-safe coercion ----------------------------

# Used in tmb_fun -- needs work, but this is a
# good approach generally

#' @export
null_to_char0 = function(x) {
  if(is.null(x)) return(character(0L))
  as.character(x)
}

#' @export
null_to_int0 = function(x) {
  if(is.null(x)) return(integer(0L))
  as.integer(x)
}

#' @export
null_to_num0 = function(x) {
  if(is.null(x)) return(numeric(0L))
  as.numeric(x)
}

#' @export
null_to_log0 = function(x) {
  if(is.null(x)) return(logical(0L))
  as.logical(x)
}

#' @export
null_to_charNA = function(x) {
  if(is.null(x)) return(as.character(NA))
  as.character(x)
}

#' @export
null_to_intNA = function(x) {
  if(is.null(x)) return(as.integer(NA))
  as.integer(x)
}

#' @export
null_to_numNA = function(x) {
  if(is.null(x)) return(as.numeric(NA))
  as.numeric(x)
}

#' @export
null_to_logNA = function(x) {
  if(is.null(x)) return(as.logical(NA))
  as.logical(x)
}

#' @export
null_to_0 = function(x) {
  if(is.null(x)) return(0L)
  as.integer(x)
}

#' @export
int0_to_0 = function(x) {
  if(length(x) == 0L) return(0L)
  as.integer(x)
}

# utilities for rate functions ----------------------------

## regex pattern for finding variables
## (e.g. any parameter or state variable)
## variable_regex looks like this '(beta0|Ca|...|zeta|S|E|Ia|...|V)'
variable_regex <- function(...) {
  return(getOption("MP_name_search_regex"))
  # only works because there are no reserved words allowed yet
  return('[A-z]([A-z][0-9])+')

  # this strategy failed (e.g. Isum matches Is)
  character_class <-
    (list(...)
     %>% lapply(names)
     %>% unlist()
     %>% paste0(collapse = "|")
    )
  paste0("(", character_class, ")", sep = "")
}
get_variables <- function(x) {
  r = gregexpr(variable_regex(), x)
  unlist(regmatches(x, r))
}
## FIXME: this only works because complements (1 - x) and
## inverses (1 / x) are so similar in structure
find_operators <- function(x, operator) {
  grepl(
    paste0("\\( *1 *", operator, " *",
           #variable_regex(params, state),
           variable_regex(),
           collapse = ""
    ), x
  )
}
factor_table <- function(x) {
  data.frame(
    var = unlist(lapply(x, get_variables)),
    compl = unlist(lapply(x, find_operators, "-")),
    invrs = unlist(lapply(x, find_operators, "/"))
  )
}

find_pos_grep <- function(tags, x) {
  pos_list = (tags
              %>% sprintf(fmt = "^%s(_|$)")
              %>% lapply(grep, x)
              %>% unlist(use.names = FALSE)
  )
}

# FIXME: this is not being used anywhere
check_in_rate_fun = function(pf) {
  stopifnot(inherits(get("model", envir = pf), 'flexmodel'))
  f = get("formula", envir = pf)
  stopifnot(
    inherits(f, "formula") |
      inherits(f, "struc") |
      inherits(f, "character")
  )
}

cross = function(from, to, mat) {
  expand.grid(
    from_pos = unique(find_pos_grep(from, rownames(mat))),
    to_pos = unique(find_pos_grep(to, colnames(mat)))
  )
}

pwise = function(from, to, mat) {
  if (length(from) == 1L) {from = rep(from, length(to))}
  if (length(to) == 1L) {to = rep(to, length(from))}
  from_pos = find_pos_grep(from, rownames(mat))
  to_pos = find_pos_grep(to, colnames(mat))
  if(length(from_pos) != length(to_pos)) {
    stop("\nargument 'from' matches to ", length(from_pos), " row indices.",
         "\nargument 'to' matches to ", length(to_pos), " column indices.",
         "\nbut these numbers must match for valid pairwise indexing of ",
         "the rate matrix")
  }
  cbind(from_pos, to_pos)
}

#' @export
block = function(from, to, mat) {
  stop('Blockwise rate specification is not implemented')
}

#' @export
check_from_to = function(from, to, state_nms) {
  bads = c(
    from[!from %in% state_nms],
    to[!to %in% state_nms]
  )
  if (length(bads) > 0L) {
    stop(
      "the following states involved in flows are not in the model:\n",
      paste0(bads, collapse = ', ')
    )
  }
}

#' Parse a Flexmodel Formula
#'
#' @param x one-sided formula, character vector, or struc object describing
#' the dependence of rate matrix elements on parameter and/or state variables
#' @return depends on the spec version (TODO: add detail once we converge)
#' @export
parse_formula = function(x) {
  (x
   %>% rateform_as_char
   %>% strsplit(split = "\\+") %>% getElement(1L)
   %>% strsplit(split = "\\*")
  )


  #pf = function(y) {
  #  (y
  #   %>% strsplit(split = "\\+") %>% getElement(1L)
  #   %>% strsplit(split = "\\*")
  #  )
  #}
  #return(pf(y))


  #o = lapply(y, pf)
  #if(spec_ver_lt('0.1.0')) {
  #  return(o[[1L]])
  #}
  #return(o)
}

#' @rdname parse_formula
#' @export
rateform_as_char = function(x) {
  y = as.character(x)
  if(inherits(x, "formula")) y = y[[2L]]
  y
}

#' @param x character vector of names to look for in \code{vec} or
#' \code{names(vec)}
#' @param vec character vector or object with names
find_vec_indices <- function(x, vec) {
  if(!is.character(vec)) vec = names(vec)
  missing_variables = x[!x %in% vec]
  if(length(missing_variables) > 0L) {
    stop("the following variables were used but not found in the model:\n",
         paste0(missing_variables, collapse = "\n"))
  }
  (x
    %>% as.character()
    %>% outer(vec, "==")
    %>% apply(1, which)
  )
}

combine_rates = function(rates) {
  rate = list()
  rate$from = rates[[1]]$from
  rate$to = rates[[1]]$to
  rate$formula = (rates
    %>% lapply(getElement, "formula")
    %>% lapply(rateform_as_char)
    %>% lapply(unlist)
    %>% c(list(sep = " + "))
    %>% do.call(what = "paste")
  )
  factors = lapply(rates, getElement, "factors")
  prod_indx = (factors
   %>% lapply(getElement, "prod_indx")
   %>% lapply(as.integer)
  )
  max_prod_indx = c(0, (rates[-length(rates)]
    %>% lapply(getElement, "factors")
    %>% lapply(getElement, "prod_indx")
    %>% lapply(as.integer)
    %>% lapply(max)
    %>% unlist(use.names = FALSE)
    %>% cumsum
  ))
  rate$factors = bind_rows(factors)
  rate$factors$prod_indx = (`+`
    %>% mapply(max_prod_indx, prod_indx, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    %>% unlist
  )
  rate$factors$prod_indx = as.character(rate$factors$prod_indx)
  rate$ratemat_indices = rates[[1]]$ratemat_indices
  lapply_any = function(name) {
    (rates
      %>% lapply(getElement, name)
      %>% unlist
      %>% any
    )
  }
  rate$state_dependent = lapply_any("state_dependent")
  rate$time_varying = lapply_any("time_varying")
  rate$sum_dependent = lapply_any("sum_dependent")
  structure(rate, class = "rate-struct")
}

# reduce rates that act on a single from-to state pair
# to a single total rate
reduce_rates = function(rates) {
  repeated_transitions = (rates
   %>% names
   %>% table
   %>% `>`(1L)
   %>% which
   %>% names
  )
  if (length(repeated_transitions) == 0L) return(rates)
  if (getOption("MP_warn_repeated_rates")) {
    warning("More than one rate was specified for at least one from-to state pair")
  }
  combined_rates = (repeated_transitions
    %>% lapply(function(name) rates[names(rates) %in% name])
    %>% lapply(combine_rates)
    %>% setNames(repeated_transitions)
  )
  rates[names(rates) %in% repeated_transitions] = NULL
  rates = c(rates, combined_rates)
  rates
}

# getting information about rates and sums ----------------------

#' @export
get_rate_info = function(model, what) lapply(model$rates, '[[', what)

#' @export
get_factr_info = function(model, what) lapply(model$factrs, '[[', what)

#' @export
get_rate_from = function(model) get_rate_info(model, 'from')

#' @export
get_rate_to = function(model) get_rate_info(model, 'to')

#' @export
get_rate_ratemat_indices = function(model) get_rate_info(model, 'ratemat_indices')

#' @export
get_rate_formula = function(model) get_rate_info(model, 'formula')

#' @export
get_rate_factors = function(model) get_rate_info(model, 'factors')

#' @export
get_rate_state_dependent = function(model) get_rate_info(model, 'state_dependent')

#' @export
get_rate_time_varying = function(model) get_rate_info(model, 'time_varying')

#' @export
get_rate_sum_dependent = function(model) get_rate_info(model, 'sum_dependent')

#' @export
get_factr_formula = function(model) get_factr_info(model, 'formula')

#' @export
get_n_products = function(model) {
  (model
   %>% get_rate_info('factors')
   %>% lapply('[[', 'prod_indx')
   %>% lapply(unique)
   %>% lapply(length)
  )
}

#' @export
get_n_variables = function(model) {
  (model
   %>% get_rate_info('factors')
   %>% lapply('[[', 'var')
   %>% lapply(length)
  )
}

#' @export
get_n_factors = function(model) {
  (model
   %>% get_rate_info('factors')
   %>% lapply(nrow)
  )
}

#' @export
get_sum_info = function(model, what) lapply(model$sums, '[[', what)

#' @export
get_sum_summands = function(model) get_sum_info(model, 'summands')

#' @export
get_sum_indices = function(model) get_sum_info(model, 'sum_indices')

#' @export
get_sum_initial_value = function(model) get_sum_info(model, 'initial_value')

#' @export
get_factr_initial_value = function(model) get_factr_info(model, 'initial_value')

#' @export
get_rate_vars = function(model) {
  (model
   %>% get_rate_factors
   %>% lapply(getElement, "var")
  )
}

#' @export
get_rates_with_vars = function(model, var_pattern) {
  ii = (model
    %>% get_rate_vars
    %>% lapply(grepl, pattern = var_pattern)
    %>% lapply(any)
    %>% sapply(isTRUE)
  )
  get_rates(model)[ii]
}

##' @param include_formula include a column for the expanded rate formula
##' @export
rate_summary = function(model, include_formula = FALSE) {
  summary = data.frame(
    from = get_rate_info(model, "from") %>% unlist,
    to = get_rate_info(model, "to") %>% unlist,
    n_fctrs = get_n_factors(model) %>% unlist,
    n_prdcts = get_n_products(model) %>% unlist,
    n_vrbls = get_n_variables(model) %>% unlist,
    state_dependent = get_rate_info(model, "state_dependent") %>% unlist,
    time_varying = get_rate_info(model, "time_varying") %>% unlist,
    sum_dependent = get_rate_info(model, "sum_dependent") %>% unlist
  )
  if(include_formula) {
    summary$formula = (model
      %>% get_rate_info('formula')
      %>% lapply(rateform_as_char)
      %>% unlist
    )
  }
  summary
}

#' @export
get_rates = function(model) {
  model$rates
}

##' @param x parameter vector or flexmodel
##' @export
has_time_varying <- function(x) {
  spec_check(
    feature = "Time-varying parameters",
    introduced_version = "0.0.3"
  )
  if(inherits(x, 'flexmodel')) {
    # more reliable test than looking at parameter attributes
    return(length(x$timevar$piece_wise$schedule$Value) > 0L)
  } else {
    return("tv_param_indices" %in% names(attributes(x)))
  }
}

##' @export
time_varying_rates <- function(model) {
  model$rates[which_time_varying_rates(model)]
}

##' @export
which_time_varying_rates <- function(model) {
  sd  <- get_rate_info(model, "state_dependent") %>% unlist
  tv  <- get_rate_info(model, "time_varying") %>% unlist
  smd <- get_rate_info(model, "sum_dependent") %>% unlist
  which(sd | tv | smd)
}

##' @export
state_dependent_rates <- function(model) {
  i = get_rate_info(model, "state_dependent") %>% unlist
  model$rates[i]
}

##' @export
sum_dependent_rates = function(model) {
  i = get_rate_info(model, "sum_dependent") %>% unlist
  model$rates[i]
}

#' @export
compute_rates = function(model) {
  (model
   %>% get_rate_info("formula")
   %>% eval_formulas(get_var_list(model))
   %>% setNames(names(model$rates))
  )
}

#' @export
compute_factrs = function(model) {
  (model
   %>% get_factr_info("formula")
   %>% eval_formulas(get_var_list(model))
   %>% setNames(names(model$factrs))
  )
}

#' @export
eval_formulas = function(formula_list, var_list) {
  (formula_list
   %>% lapply(function(x) ifelse(inherits(x, 'formula'), as.character(x[2]), x))
   %>% sapply(as.character)
   %>% struc
   %>% struc_eval(var_list)
   %>% c()
  )
}

#' @export
get_var_list = function(model) {
  c(
    as.list(model$params),
    as.list(model$state),
    as.list(model$sum_vector),
    as.list(model$factr_vector)
  )
}

compute_rate_from_indices = function(model, i) {
  unpack(get_indices_per_rate(model, i))
  from = model$tmb_indices$make_ratemat_indices$from[i]
  to = model$tmb_indices$make_ratemat_indices$to[i]
  count = model$tmb_indices$make_ratemat_indices$count[i]
  sp = c(model$state, model$param, model$sum_vector)
  factr_algo(spi, sp, modifier)
}

compute_factr_from_indices = function(model, i) {
  stop("not implemented")
}

factr_algo = function(spi, sp, modifier) {
  result = 0
  prod = 1.0
  for (j in seq_along(spi)) {
    x = sp[spi[j]]
    if (modifier[j] > 3) {
      result = result + prod
      prod = 1
    }
    if (modifier[j] %in% c(1, 3, 5, 7)) {
      x = 1-x
    } else if (modifier[j] %in% c(2, 3, 6, 7)){
      if (x != 0) {
        x = 1/x
      }
    }
    prod = prod * x
  }
  result = result + prod
  return(result)
}

get_indices_per_rate = function(model, i) {
  unpack(model$tmb_indices$make_ratemat_indices)
  start = sum(c(1, count)[1:i])
  end = start + count[i] - 1L
  list(spi = spi[start:end],
       modifier = modifier[start:end],
       start = start, end = end)
}

#' @export
lookup_pairwise = function(from, to, M) {
  i = pwise(from, to, M)
  data.frame(
    from = rownames(M)[i[,"from_pos"]],
    to = colnames(M)[i[,"to_pos"]]
  )
}


#' Rate Matrix Loopup Table
#'
#' @param state state_pansim object
#' @param ratemat rate matrix
#' @export
rate_matrix_lookup = function(ratemat) {
  ratemat = as(ratemat, "dgTMatrix")
  (data.frame(
    from_pos = ratemat@i + 1,
    to_pos = ratemat@j + 1)
    %>% mutate(
      from_state = ratemat@Dimnames[[1]][from_pos],
      to_state = ratemat@Dimnames[[2]][to_pos]
    )
  )
}





# computing indices for tmb ----------------------

##' @export
sum_indices = function(sums, state, params) {
  sum_index_list = lapply(sums, "[[", "sum_indices")
  sumidx = c(seq_len(length(sums))) + length(state) + length(params)
  sumcount = sapply(sum_index_list, length, USE.NAMES = FALSE) %>% unname
  summandidx = c(unlist(sum_index_list, use.names = FALSE))
  clean_if_empty = function(x) {
    if((length(x) == 0L) | is.null(x)) return(integer(0L))
    x
  }
  list(sumidx = clean_if_empty(sumidx),
       sumcount = clean_if_empty(sumcount),
       summandidx = clean_if_empty(summandidx))
}

##' @export
factr_indices = function(factrs, state_param_sums) {
  if (length(factrs) == 0L) {
    indices = list(
      spi_factr = integer(0L),
      count = integer(0L),
      spi = integer(0L),
      modifier = integer(0L)
    )
    return(indices)
  }
  sp = state_param_sums
  spi_factr = (factrs
    %>% seq_along
    %>% `+`(length(sp))
  )
  spi <- get_var_indx(factrs)
  count <- get_var_counts(factrs)
  modifier <- get_var_modifiers(factrs)
  names(spi) = names(count)  = names(count) = NULL
  indices <- list(
    spi_factr = spi_factr,
    count = count,
    spi = spi,
    modifier = modifier
  )
  return(indices)
}

##' @export
ratemat_indices <- function(rates, state_params) {
  sp <- state_params
  ratemat_indices <- sapply(rates, `[[`, "ratemat_indices")
  spi <- get_var_indx(rates)
  count <- get_var_counts(rates)
  modifier <- get_var_modifiers(rates)
  names(spi) <- colnames(ratemat_indices) <- names(count) <- NULL
  indices <- list(
    from = ratemat_indices[1, ],
    to = ratemat_indices[2, ],
    count = count,
    spi = spi,
    modifier = modifier
  )
  return(indices)
}

get_var_indx = function(l) {
  {
    lapply(l, function(y) {
      y$factors$var_indx
    }) %>% unlist()
  }
}

get_var_counts = function(l) {
  sapply(l, function(y) {
    nrow(y$factors)
  })
}

get_var_modifiers = function(l) {
  # an 'item' is either a rate or (common) factr
  (l
   %>% unname
   %>% lapply("[[", "factors")
   %>% bind_rows(.id = "item_indx")
   %>% mutate(new_item = as.logical(c(0, diff(as.numeric(item_indx)))))
   %>% mutate(new_prod = as.logical(c(0, diff(as.numeric(prod_indx)))))
   %>% mutate(add = new_prod & (!new_item))
   %>% mutate(modifier = 4 * add + 2 * invrs + compl)
   %>% `$`("modifier")
  )
}

##' @family flexmodels
##' @export
state_mapping_indices = function(
  state,
  eigen_drop_pattern,
  infected_drop_pattern,
  susceptible_pattern) {
  spec_check(introduced_version = "0.1.1",
             feature = "Disease free state updates")

  if(length(susceptible_pattern) == 0L) {
    susceptible_idx = integer(0L)
  } else {
    susceptible_idx = grep(susceptible_pattern, names(state), perl = TRUE)
  }

  eigen = eigen_drop_pattern
  infected = infected_drop_pattern
  all = names(state)
  if((length(eigen) == 0L) | (length(infected) == 0L)) {
    index_set = list(
      x = list(vector('list', 3L)),
      i = list(vector('list', 3L)),
      j = list(vector('list', 3L)))
  } else {
    index_set = make_nested_indices(all, nlist(eigen, infected), invert = TRUE)
  }

  c(
    (index_set
     %>% getElement('i')
     %>% unlist(recursive = FALSE)
     %>% setNames(c("all_to_eigen_idx",
                    "all_to_infected_idx",
                    "eigen_to_infected_idx"))
    ),
    (index_set
     %>% getElement('j')
     %>% unlist(recursive = FALSE)
     %>% setNames(c("all_drop_eigen_idx",
                    "all_drop_infected_idx",
                    "eigen_drop_infected_idx"))
    ),
    nlist(susceptible_idx)
  )
}

##' @export
disease_free_indices = function(model) {
  unpack(model)

  df_state_par_idx = (disease_free
    %>% lapply(`[`, 'params_to_use')
    %>% unlist(use.names = FALSE)
    %>% find_vec_indices(params)
  )
  df_state_count = (disease_free
    %>% lapply(`[[`, 'states_to_update')
    %>% lapply(length)
    %>% unlist
  )
  df_state_idx = (disease_free
    %>% lapply(`[`, 'states_to_update')
    %>% unlist(use.names = FALSE)
    %>% find_vec_indices(state)
  )
  nlist(df_state_par_idx, df_state_count, df_state_idx)
}

##' @export
initial_population_indices = function(model) {
  unpack(model)
  infected_idx = which(initial_population$infected == names(params))
  total_idx = which(initial_population$total == names(params))
  nlist(total_idx, infected_idx)
}

##' @export
initialization_mapping_indices = function(model) {
  unpack(model)
  state_mapping_indices(
    state,
    initialization_mapping$eigen,
    initialization_mapping$infected,
    initialization_mapping$susceptible)
}

##' @export
linearized_param_indices = function(model) {
  unpack(model)
  lin_param_vals = (linearized_params
                   %>% lapply(`[`, 'update_value')
                   %>% unlist(use.names = FALSE)
  )
  lin_param_count = (linearized_params
                    %>% lapply(`[[`, 'params_to_update')
                    %>% lapply(length)
                    %>% unlist(use.names = FALSE)
  )
  lin_param_idx = (linearized_params
                  %>% lapply(`[`, 'params_to_update')
                  %>% unlist(use.names = FALSE)
                  %>% find_vec_indices(params)
  )

  nlist(lin_param_vals, lin_param_count, lin_param_idx)
}

##' @export
outflow_indices = function(outflow, ratemat) {
  if(length(outflow) == 0L) {
    return(list(
      row_count = integer(),
      col_count = integer(),
      rows = integer(),
      cols = integer()
    ))
  }
  indices = lapply(outflow, function(o) {
    list(
      state = grep(o$from, rownames(ratemat), perl = TRUE),
      flow = grep(o$to, colnames(ratemat), perl = TRUE)
    )
  })
  n = length(indices)
  indices$row_count = seq(n)
  indices$col_count = seq(n)
  indices$rows = vector()
  indices$cols = vector()

  for(i in 1:n) {
    rows = indices[[i]]$state
    cols = indices[[i]]$flow

    indices$row_count[i] = length(rows)
    indices$col_count[i] = length(cols)

    indices$rows = c(indices$rows, rows)
    indices$cols = c(indices$cols, cols)
  }

  indices

  # TODO: add warning if any state_indices are equal to other state_indices
  #       that are already added to some other call to outflow in the model.
  #       in general state_indices should be mutually exclusive across each
  #       outflow.
  #nlist(state_indices, flow_state_indices)
  #all = names(model$state)
  #lapply(model$outflow, McMasterPandemic:::make_nested_indices, x = all)
}

#' Nested Indices
#'
#' Create and return a nested set of character vectors from a sequence
#' of regular expressions, and return indices into each vector for recovering
#' other shorter vectors that are lower in the hierarchy.
#'
#' @section Motivating Example
#'
#' There are three nested state vectors
#' a_states -- all states
#' p_states -- excludes accumulators
#' e_states -- includes only infected states
#'
#' There are three index vectors
#' i_ap -- indexes a_states to yield p_states
#' i_ae -- indexes a_states to yield e_states
#' i_pe -- indexes p_states to yield e_states
#'
#' There are also three inverse index vectors
#' j_ap -- indexes a_states to yield setdiff(a_states, p_states)
#' j_ae -- indexes a_states to yield setdiff(a_states, e_states)
#' j_pe -- indexes p_states to yield setdiff(p_states, e_states)
#'
#' In general -- assume that y_states are nested in x_states
#' i_xy -- indexes x_states to yield y_states
#' j_xy -- indexes x_states to yeild setdiff(x_states, y_states)
#' i_yx -- doesn't exist because not all x_states are also y_states
#'
#' @param x character vector
#' @param patterns named list or vector of regular expressions to be applied
#' sequentially to create a nested set of character vectors
#' @return List with three elements.
#' \describe{
#'   \item{\code{x}}{Nested character vectors.}
#'   \item{\code{i}}{
#'     List of lists, with inner list giving the indices into character vectors
#'     that are higher in the hierarchy. For example, one may read this
#'     expression, \code{x$e == x$p\\[i$e$p\\]}, as "e equals p at the index
#'     that takes p to e".
#'   }
#'   \item{\code{j}}{Set difference version of \code{i}.}
#' }
make_nested_indices = function(x, patterns, invert = FALSE) {

  stopifnot(all(sapply(patterns, is.character)))
  stopifnot(all(sapply(patterns, length) == 1L))
  stopifnot(!any(is.null(names(patterns))))

  nms = c(deparse(substitute(x)), names(patterns))

  x_list = list(x)
  i_list = j_list = list()

  # p indexes list of patterns
  # v indexes list of vectors (i.e. list of nested subsets of x)
  for(p in seq_along(patterns)) {
    x_list[[p+1]] = grep(patterns[[p]], x_list[[p]], value = TRUE, perl = TRUE, invert = invert)
    #i_list[[to]][[from]]
    i_list[[p+1]] = j_list[[p+1]] = list()

    for(v in seq_len(p)) {
      indicators = x_list[[v]] %in% x_list[[p+1]]
      i_list[[p+1]][[v]] = which( indicators)
      j_list[[p+1]][[v]] = which(!indicators)

    }
    names(i_list[[p+1]]) = names(j_list[[p+1]]) = nms[1:p]
  }

  drop_null_elements = function(e) e[!sapply(e, is.null)]
  list(
    x = setNames(x_list, nms) %>% drop_null_elements,
    i = setNames(i_list, nms) %>% drop_null_elements,
    j = setNames(j_list, nms) %>% drop_null_elements)
}

lin_state_timevar_params = function(schedule) {
  stop('function lin_state_timevar_params is under construction')
}


# retrieving information from tmb objective function --------------

#' Extract Parameter Vector to Pass to a TMB Function
#'
#' Get all of the parameters from a flexmodel in a vector that is
#' ready to be passed to a TMB AD objective function, gradient function,
#' hessian function, simulate function, report function.
#'
#' Currently this includes \code{params} and if \code{spec_ver_gt('0.1.0')}
#' it also includes time-varying multipliers in
#' \code{model$timevar$piece_wise$schedule$last_tv_mult}.
#'
#' @param model flexmodel
#' @export
tmb_params = function(model) {
  # TODO: when we start using the TMB map argument to only pass parameters
  # that are allowed to change, we will need to account for this here
  full_param_vec = c(unlist(model$params))
  if(spec_ver_gt('0.1.0')) {
    if(has_time_varying(model)) {
      last_tv_mult = model$timevar$piece_wise$schedule$last_tv_mult
      if (any(is.na(last_tv_mult)))
        stop('missing time-varying multipliers need to be replaced for use with TMB')
      full_param_vec = c(full_param_vec, last_tv_mult)
    }
  }
  full_param_vec
}

#' @export
simulation_dates = function(model) {
  seq.Date(
    as.Date(model$start_date),
    as.Date(model$end_date),
    by = 1)
}

#' @param model flexmodel
#' @param sim_params parameter vector to pass to a TMB objective function
#' @export
changing_ratemat_elements = function(model, sim_params = NULL) {
  if(is.null(sim_params)) sim_params = tmb_params(model)
  tmb_fun(model)$simulate(sim_params)$concatenated_ratemat_nonzeros
}

#' @export
simulate_changing_ratemat_elements = function(model, sim_params = NULL, add_dates = FALSE) {
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

#' @export
concatenated_state_vector = function(model, sim_params = NULL) {
  if(is.null(sim_params)) sim_params = tmb_params(model)
  tmb_fun(model)$simulate(sim_params)$concatenated_state_vector
}

#' @export
structure_state_vector = function(x, iters, state_nms) {
  matrix(
    x,
    nrow = iters + 1L,
    ncol = length(state_nms),
    dimnames = list(1:(iters + 1), state_nms),
    byrow = TRUE
  )
}

#' @export
simulate_state_vector = function(model, sim_params = NULL, add_dates = FALSE) {
  state_sims = (model
   %>% concatenated_state_vector(sim_params)
   %>% structure_state_vector(model$iters, names(model$state))
   %>% as.data.frame
  )
  if(add_dates) {
    state_sims = (model
      %>% simulation_dates
      %>% data.frame
      %>% setNames("Date")
      %>% cbind(state_sims)
    )
  }
  state_sims
}

#' @export
initial_state_vector = function(model, sim_params = NULL) {
  # FIXME: minor performance optimization could be made
  # here by restricting the number of iterations temporarily
  # so that the full simulation is not run just to get the
  # first value -- i think all that is needed is this (but
  # need to test):
  # model$iters = 1
  (model
   %>% concatenated_state_vector(sim_params)
   %>% head(length(model$state))
   %>% setNames(names(model$state))
  )
}

#' @export
final_state_vector = function(model, sim_params = NULL) {
  (model
   %>% concatenated_state_vector(sim_params)
   %>% tail(length(model$state))
   %>% setNames(names(model$state))
  )
}

#' @export
penultimate_state_vector = function(model, sim_params = NULL) {
  (model
   %>% concatenated_state_vector(sim_params)
   %>% tail(2 * length(model$state))
   %>% head(length(model$state))
   %>% setNames(names(model$state))
  )
}

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

#' @export
initial_ratemat = function(model, sim_params = NULL) {
  if(is.null(sim_params)) sim_params = tmb_params(model)
  tmb_fun(model)$simulate(sim_params)$ratemat
}

# benchmarking and comparison in tests --------------------

#' Compare TMB-based and classic MacPan Simulations
#'
#' Exceptions to drop-in replacement:
#' 1. don't require that the attributes are in the same order
#' 2. don't require that the r version returns everything that the tmb
#'    version does (e.g. flexmodel)
#' 3. don't require that the row.names are identical (is this ok?
#'    the r version counts iterations with skips, but is this informative?)
#' 4. don't require that the call is identical (obvious i guess, but
#'    being exhaustive)
#' 5. don't look for parameter called S0 in classic simulation results
#'    (FIXME: this one might not be necessary anymore)
#' 6. don't require the ordering of params_timevar rows to be the same
#'    (FIXME: this should probably be handled in run_sim itself)
#'
#' @param classic_sim result of `run_sim` without using TMB
#' @param tmb_sim result of `run_sim` using TMB
#' @param tolerance numerical tolerance
#' @param compare_attr compare attributes or just the simulations themselves
#' @export
compare_sims = function(classic_sim, tmb_sim, tolerance = NULL, compare_attr = TRUE) {
  if (is.null(tolerance)) {
    if (require(testthat)) {
      tolerance = testthat_tolerance()
    } else {
      tolerance = .Machine$double.eps^0.5
    }
  }
  if(compare_attr) {
    params_to_keep = which(names(attr(tmb_sim, 'params')) != "S0")
    attr(tmb_sim, 'params')[params_to_keep] = attr(tmb_sim, 'params')[params_to_keep]
    attr(tmb_sim, 'row.names') = attr(classic_sim, 'row.names') = NULL
    attr(tmb_sim, 'call') = attr(classic_sim, 'call') = NULL

    # FIXME: this should probably be handled in run_sim itself
    params_timevar = attr(tmb_sim, "params_timevar")
    if(!is.null(params_timevar)) {
      attr(tmb_sim, "params_timevar") = arrange(params_timevar, Date, Symbol)
      attr(attr(tmb_sim, "params_timevar"), "row.names") = NULL
      attr(attr(classic_sim, "params_timevar"), "row.names") = NULL
    }

    for(a in names(attributes(classic_sim))) {
      expect_equal(attr(tmb_sim, a), attr(classic_sim, a), tolerance = tolerance)
    }
  }
  attributes(classic_sim) <- attributes(tmb_sim) <- NULL
  expect_equal(tmb_sim, classic_sim, tolerance = tolerance)
  TRUE
}


#' Time Two Expressions
#'
#' Evaluate two expressions in the environment in which
#' \code{time_wrap} is called, and record the system
#' time associated with each expression.
#'
#' @param expr1 first expression
#' @param expr2 second expression
#' @param units time units with which to measure the system time
#' of each expression (see \code{\link{difftime}})
#' @return list with three components: (1) time to evaluate
#' the first expression, (2) time to evaluate the second expression,
#' and (3) the time for the first expression divided by the time
#' for the second
#'
#' @export
time_wrap = function(expr1, expr2, units = 'secs') {
  e = parent.frame()

  strt1 = Sys.time()
  eval(expr1, e)
  nd1 = Sys.time()
  dff1 = difftime(nd1, strt1, units = units)

  strt2 = Sys.time()
  eval(expr2, e)
  nd2 = Sys.time()
  dff2 = difftime(nd2, strt2, units = units)

  list(dff1, dff2, as.numeric(dff2)/as.numeric(dff1))
}

#' @export
compare_sims_cbind = function(classic_sim, tmb_sim, index) {
  stopifnot(length(index) == 1L)
  x = cbind(unlist(classic_sim[,index]),
        unlist(tmb_sim[,index]))
  colnames(x) = c('classic', 'tmb') %_% names(classic_sim[index])
  x
}

#' @export
compare_sims_plot = function(...) {
  x = compare_sims_cbind(...)
  ss = sub('^classic_', '', colnames(x))
  plot(x[,1], type = 'l', main = ss[[1L]][1])
  lines(x[,2], col = 'red', lty = 2)
}

#' @export
compare_sims_diffmat = function(classic_sim, tmb_sim, state, op = `-`) {
  op(
    as.matrix(classic_sim[,names(state)]),
    as.matrix(tmb_sim[,names(state)]))
}

#' @export
compare_sims_unlist = function(classic_sim, tmb_sim, tolerance = testthat_tolerance()) {
  all.equal(unlist(tmb_sim), unlist(classic_sim), tolerance)
  TRUE
}

# spec version and global option control ----------------------

##' Set and Reset Spec Version
##'
##' Set the spec version and optionally compile a
##' specific c++ file associated with this version.
##'
##' The user need to have write permissions to
##' \code{cpp_path}, because object and shared
##' object files will be created in this
##' directory.
##'
##' \code{reset_spec_version} returns specs
##' being assumed and c++ file being used to
##' factory-fresh settings.
##'
##' @param v character string with spec version
##' @param cpp_path string containing path
##' to a directory containing a file called
##' \code{macpan.cpp}, which is to be used
##' to construct the objective function --
##' the default is \code{NULL}, which will
##' make use of the packaged source file
##' @param use_version_directories is cpp_path
##' organized by sub-directories with names
##' given by the spec version number? (as
##' they are in \code{inst/tmb})
##' \code{McMasterPandemic.cpp}
##'
##' @rdname set_spec_version
##' @importFrom tools file_path_sans_ext
##' @export
set_spec_version = function(v, cpp_path = NULL, use_version_directories = TRUE) {
  spec_version <- as.character(unlist(v))[1]
  print(spec_version)
  options(MP_flex_spec_version = spec_version)

  if(is.null(cpp_path)) {
    options(MP_flex_spec_dll = "McMasterPandemic")
  } else {
    sub_dir = ifelse(use_version_directories, spec_version, '')
    cpp <- file.path(cpp_path, sub_dir, "macpan.cpp")
    dll <- file_path_sans_ext(cpp)
    options(MP_flex_spec_dll = basename(dll))
    compile(cpp)
    dyn.load(dynlib(dll))
  }
}

#' @rdname set_spec_version
#' @export
reset_spec_version = function() {
  flex_version <- readLines(
    system.file("tmb/recommended_spec_version",
                package = "McMasterPandemic"))
  options(MP_flex_spec_version = flex_version)
  options(MP_flex_spec_dll = "McMasterPandemic")
}

##' TMB Mode
##'
##' Set options so that R engine runs in a manner
##' that is comparable with the TMB engine
##'
##' @export
tmb_mode = function() {
  options(
    MP_use_state_rounding = FALSE,
    MP_vax_make_state_with_hazard = FALSE)
}

##' R Mode
##'
##' Set options so that TMB engine runs in a manner
##' that is comparable with the R engine
##'
##' @export
r_mode = function() {
  options(
    MP_tmb_models_match_r = TRUE
  )
}

##' Compare R and TMB Engines
##'
##' Set options so that the TMB and R engines
##' can be compared.
##'
##' This is useful for testing.
##'
##' @export
r_tmb_comparable = function() {
  r_mode()
  tmb_mode()
}

##' Get McMasterPandemic Options
##'
##' @export
get_macpan_options = function() {
  (options()
   %>% names
   %>% grep(pattern = "^MP_", value = TRUE)
   %>% sapply(getOption, simplify = FALSE)
  )
}


# converting between classic and tmb params_timevar -------
#   this is annoying:
#   macpan_ontario orders by Symbol then Date
#   macpan.cpp orders by breaks then spi (which is roughly Date then Symbol)
#' @export
timevar_order_indices = function(schedule) {
  # FIXME: need to do something different if rel_orig/rel_prev are both used?

  # arrange classic method
  schedule = arrange(schedule, Symbol, Date)
  schedule$classic_row_num = seq_len(nrow(schedule))

  # arrange tmb method
  schedule = arrange(schedule, breaks, tv_spi)
  schedule$classic_order = order(schedule$classic_row_num)
  schedule$tmb_row_num = seq_len(nrow(schedule))

  # arrange back to classic
  schedule = arrange(schedule, Symbol, Date)
  schedule$tmb_order = order(schedule$tmb_row_num)

  select(schedule, classic_order, tmb_order)
}
