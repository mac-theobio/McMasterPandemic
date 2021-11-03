# constructing names and strings ----------------------

#' Paste with Underscore Separator
#' @export
`%_%` = function(x, y) paste(x, y, sep = "_")

#' Paste with Blank Separator
#'
#' Like Python string `+`
#' @export
`%+%` = function(x, y) paste(x, y, sep = "")

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

#' Parse a Flexmodel Formula
#'
#' @param x one-sided formula, character vector, or struc object describing
#' the dependence of rate matrix elements on parameter and/or state variables
#' @return depends on the spec version (TODO: add detail once we converge)
#' @export
parse_formula = function(x) {
  y = as.character(x)
  if(inherits(x, "formula")) y = y[[2L]]
  pf = function(y) {
    (y
     %>% strsplit(split = "\\+") %>% getElement(1L)
     %>% strsplit(split = "\\*")
    )
  }
  return(pf(y))
  o = lapply(y, pf)
  if(spec_ver_lt('0.1.0')) {
    return(o[[1L]])
  }
  return(o)
}

find_vec_indices <- function(x, vec) {
  (
    x
    %>% as.character()
    %>% outer(names(vec), "==")
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


##' @export
rate_summary = function(model) {
  data.frame(
    from = get_rate_info(model, "from") %>% unlist,
    to = get_rate_info(model, "to") %>% unlist,
    n_fctrs = get_n_factors(model) %>% unlist,
    n_prdcts = get_n_products(model) %>% unlist,
    n_vrbls = get_n_variables(model) %>% unlist,
    state_dependent = get_rate_info(model, "state_dependent") %>% unlist,
    time_varying = get_rate_info(model, "time_varying") %>% unlist,
    sum_dependent = get_rate_info(model, "sum_dependent") %>% unlist
  )
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

##' @param x parameter vector
##' @export
has_time_varying <- function(x) {
  spec_check(
    feature = "Time-varying parameters",
    introduced_version = "0.0.3"
  )
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


# retrieving information from tmb objective function --------------

#' @export
simulate_changing_ratemat_elements = function(model) {
  updateidx = model$tmb_indices$updateidx
  ratemat_elements =
    matrix(tmb_fun(model)$simulate()$concatenated_ratemat_nonzeros,
           nrow = model$iters + 1L,
           ncol = length(updateidx),
           byrow = TRUE)
  colnames(ratemat_elements) = names(updateidx)
  as.data.frame(ratemat_elements)
}

#' @export
simulate_state_vector = function(model) {
  matrix(tmb_fun(model)$simulate()$concatenated_state_vector,
         nrow = model$iters + 1L,
         ncol = length(model$state),
         byrow = TRUE) %>% as.data.frame %>% setNames(names(model$state))
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
#'
#' @param classic_sim result of `run_sim` without using TMB
#' @param tmb_sim result of `run_sim` using TMB
#' @export
compare_sims = function(classic_sim, tmb_sim) {
  attr(tmb_sim, 'row.names') = attr(classic_sim, 'row.names') = NULL
  attr(tmb_sim, 'call') = attr(classic_sim, 'call') = NULL
  for(a in names(attributes(classic_sim))) {
    expect_equal(attr(tmb_sim, a), attr(classic_sim, a))
  }

  attributes(classic_sim) <- attributes(tmb_sim) <- NULL
  expect_equal(tmb_sim, classic_sim)
  TRUE
}


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
