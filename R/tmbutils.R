# constructing names and strings ----------------------

#' Paste with Underscore Separator
#' @export
`%_%` = function(x, y) paste(x, y, sep = "_")

#' Paste with Blank Separator
#'
#' Like Python string `+`
#' @export
`%+%` = function(x, y) paste(x, y, sep = "")

##' Regex Alternation Group
alt_group = function(x) {
  "(" %+% paste0(x, collapse = "|") %+% ")"
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

#' @rdname parse_formula
#' @export
rateform_as_char = function(x) {
  y = as.character(x)
  if(inherits(x, "formula")) y = y[[2L]]
  y
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

#' @param x character vector of names to look for in \code{vec} or
#' \code{names(vec)}
#' @param vec character vector or object with names
find_vec_indices <- function(x, vec) {
  if(!is.character(vec)) vec = names(vec)
  (x
    %>% as.character()
    %>% outer(vec, "==")
    %>% apply(1, which)
  )
}



# getting information about rates and sums ----------------------

#' @export
get_rate_info = function(model, what) lapply(model$rates, '[[', what)

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
      # TODO: this probably won't work for structured formulae
      %>% lapply(getElement, 2L)
      %>% as.character
    )
  }
  summary
}

#' @export
get_rates = function(model) {
  from = flex_from(model)
  to = flex_to(model)
  mapply(function(from, to) {
    model$ratemat[from, to]
  }, from = from, to = to,
  SIMPLIFY = TRUE)
}

##' @param x parameter vector or flexmodel
##' @export
has_time_varying <- function(x) {
  spec_check(
    feature = "Time-varying parameters",
    introduced_version = "0.0.3"
  )
  if(inherits(x, 'flexmodel')) x = x$params
  "tv_param_indices" %in% names(attributes(x))
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

compute_rates = function(model) {
  (model
   %>% get_rate_info("formula")
   %>% lapply(function(x) ifelse(inherits(x, 'formula'), as.character(x[2]), x))
   %>% sapply(as.character)
   %>% struc
   %>% struc_eval(c(as.list(model$params), as.list(model$state), as.list(model$sum_vector)))
   %>% c()
   %>% setNames(names(model$rates))
  )
}

compute_rate_from_indices = function(model, i) {
  unpack(get_indices_per_rate(model, i))
  from = model$tmb_indices$make_ratemat_indices$from[i]
  to = model$tmb_indices$make_ratemat_indices$to[i]
  count = model$tmb_indices$make_ratemat_indices$count[i]
  sp = c(model$state, model$param, model$sum_vector)
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
check_rates = function(model, eps = 1e-5) {
  (data.frame(get_rates(test_model), compute_rates(test_model))
   %>% setNames(c('get', 'compute'))
   %>% mutate(diff = abs(get - compute))
   %>% mutate(bads = (diff > eps) | is.nan(compute) | is.na(get))
   %>% filter(bads)
  )
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
ratemat_indices <- function(rates, state_params) {
  sp <- state_params
  ratemat_indices <- sapply(rates, `[[`, "ratemat_indices")
  spi <- {
    lapply(rates, function(y) {
      y$factors$var_indx
    }) %>% unlist()
  }
  count <- sapply(rates, function(y) {
    nrow(y$factors)
  })
  modifier <- (rates
               %>% unname
               %>% lapply("[[", "factors")
               %>% bind_rows(.id = "rate_indx")
               %>% mutate(new_rate = as.logical(c(0, diff(as.numeric(rate_indx)))))
               %>% mutate(new_prod = as.logical(c(0, diff(as.numeric(prod_indx)))))
               %>% mutate(add = new_prod & (!new_rate))
               %>% mutate(modifier = 4 * add + 2 * invrs + compl)
               %>% `$`("modifier")
  )
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

  df_state_par_idx = (disease_free$state$simple
    %>% lapply(`[`, 'params_to_use')
    %>% unlist(use.names = FALSE)
    %>% find_vec_indices(params)
  )
  df_state_count = (disease_free$state$simple
    %>% lapply(`[[`, 'states_to_update')
    %>% lapply(length)
    %>% unlist
  )
  df_state_idx = (disease_free$state$simple
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
    initialization_mappings$eigen,
    initialization_mappings$infected,
    initialization_mappings$susceptible)
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
      state = grep(o$state_patterns, rownames(ratemat)),
      flow = grep(o$flow_state_patterns, colnames(ratemat))
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
  full_param_vec = model$params
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

#' @param model flexmodel
#' @param sim_params parameter vector to pass to a TMB objective function
#' @export
changing_ratemat_elements = function(model, sim_params = NULL) {
  if(is.null(sim_params)) sim_params = tmb_params(model)
  tmb_fun(model)$simulate(sim_params)$concatenated_ratemat_nonzeros
}

#' @export
simulate_changing_ratemat_elements = function(model, sim_params = NULL) {
  updateidx = model$tmb_indices$updateidx
  ratemat_elements = (model
    %>% changing_ratemat_elements(sim_params)
    %>% matrix(nrow = model$iters + 1L,
               ncol = length(updateidx),
               byrow = TRUE)
  )

  colnames(ratemat_elements) = names(updateidx)
  as.data.frame(ratemat_elements)
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
simulate_state_vector = function(model, sim_params = NULL) {
  (model
   %>% concatenated_state_vector(sim_params)
   %>% structure_state_vector(model$iters, names(model$state))
   %>% as.data.frame
  )
}

#' @export
initial_state_vector = function(model, sim_params = NULL) {
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

# benchmarking and comparison in tests

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
#'
#' @param classic_sim result of `run_sim` without using TMB
#' @param tmb_sim result of `run_sim` using TMB
#' @param tolerance numerical tolerance
#' @importFrom testthat testthat_tolerance
#' @export
compare_sims = function(classic_sim, tmb_sim, tolerance = testthat_tolerance()) {
  params_to_keep = which(names(attr(tmb_sim, 'params')) != "S0")
  attr(tmb_sim, 'params')[params_to_keep] = attr(tmb_sim, 'params')[params_to_keep]
  attr(tmb_sim, 'row.names') = attr(classic_sim, 'row.names') = NULL
  attr(tmb_sim, 'call') = attr(classic_sim, 'call') = NULL
  for(a in names(attributes(classic_sim))) {
    expect_equal(attr(tmb_sim, a), attr(classic_sim, a), tolerance = tolerance)
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

##' Set Spec Version
##'
##' Set the spec version and optionally compile a
##' specific c++ file associated with this version.
##'
##' The user need to have write permissions to
##' \code{cpp_path}, because object and shared
##' object files will be created in this
##' directory.
##'
##' @param v character string with spec version
##' @param cpp_path string containing path
##' to a directory containing a file called
##' \code{macpan.cpp}, which is to be used
##' to construct the objective function --
##' the default is \code{NULL}, which will
##' make use of the packaged source file
##' \code{McMasterPandemic.cpp}
##'
##' @importFrom tools file_path_sans_ext
##' @export
set_spec_version = function(v, cpp_path = NULL) {
  spec_version <- as.character(unlist(v))[1]
  print(spec_version)
  options(MP_flex_spec_version = spec_version)

  if(is.null(cpp_path)) {
    options(MP_flex_spec_dll = "McMasterPandemic")
  } else {
    cpp <- file.path(cpp_path, spec_version, "macpan.cpp")
    dll <- file_path_sans_ext(cpp)
    options(MP_flex_spec_dll = basename(dll))
    compile(cpp)
    dyn.load(dynlib(dll))
  }
}
