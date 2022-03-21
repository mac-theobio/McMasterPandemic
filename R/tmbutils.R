# objects ----------------------

valid_loss_functions = c('negative_binomial')
valid_prior_families = c('flat', 'normal')
valid_trans = c('log', 'log10', 'logit', 'cloglog', 'inverse')

init_initialization_mapping = list(
  eigen = character(0L),
  infected = character(0L),
  susceptible = character(0L)
)

init_initial_population = list(
  total = character(0L),
  infected = character(0L)
)

init_factr_vector = numeric(0L)

init_condensation = list(
  include = list(),
  # ordered list of condensation steps
  steps = list()
)

init_condensation_map = character(0L)

init_observed = list(
  data = data.frame(
    date = Date(),
    var = character(),
    value = numeric()
  ),
  loss_params = data.frame(
    Parameter = character(),
    Distribtion = character(),
    Variable = character()
  )
)

init_tmb_indices = list(
  make_ratemat_indices = list(
    from = integer(0L),
    to = integer(0L),
    count = integer(0L),
    spi = integer(0L),
    modifier = integer(0L)
  ),
  par_accum_indices = integer(0L),
  updateidx = integer(0L),
  sum_indices = list(
    sumidx = integer(0L),
    sumcount = integer(0L),
    summandidx = integer(0L)
  ),
  factr_indices = list(
    spi_factr = integer(0L),
    count = integer(0L),
    spi = integer(0L),
    modifier = integer(0L)
  ),
  condense_indices = list(
    sri_output = integer(0L),
    sr_count = integer(0L),
    sri = integer(0L),
    sr_modifier = integer(0L)
  ),
  opt_params = list(
    index_table = data.frame(
      param_nms = character(0),
      param_trans = character(0),
      prior_distr = character(0),
      prior_trans = character(0),
      count_hyperparams = integer(0),
      init_trans_params = numeric(0),
      prior_trans_id = integer(0),
      param_trans_id = integer(0),
      prior_distr_id = integer(0),
      opt_param_id = integer(0)
    ),
    hyperparameters = numeric(0L),
    index_tv_table = data.frame(
      param_nms = character(0),
      param_trans = character(0),
      prior_distr = character(0),
      prior_trans = character(0),
      count_hyperparams = integer(0),
      init_trans_params = numeric(0),
      prior_trans_id = integer(0),
      param_trans_id = integer(0),
      prior_distr_id = integer(0),
      opt_tv_mult_id = integer(0)
    ),
    hyperparameters_tv = numeric(0L)
  )
)




# test functions ------------------------------------------------

is_len1_char = function(x) (length(x) == 1L) & is.character(x)

assert_len1_char = function(x) {
  nm = deparse(substitute(x))
  if(!is_len1_char(x)) {
    stop(nm, ' must be a length-1 character')
  }
}

is_len1_int = function(x) {
  (length(x) == 1L) & isTRUE(all.equal(x, as.integer(x)))
}

assert_len1_int = function(x) {
  nm = deparse(substitute(x))
  if(!is_len1_int(x)) {
    stop(nm, ' must be a length-1 integer')
  }
}

#' @export
exists_opt_params = function(model) {
  (model
   $  tmb_indices
   $  ad_fun_map
   %>% unlist
   %>% is.na
   %>% all
   %>% `!`
  )
}

# constructing names and strings ----------------------

#' Paste with Underscore Separator
#' @export
`%_%` = function(x, y) paste(x, y, sep = "_")

#' Paste with Blank Separator
#'
#' Like Python string `+`
#'
#' @export
`%+%` = function(x, y) paste(x, y, sep = "")

#' @export
is_empty = function(x) {
  x = as.character(x)
  is.na(x) | (nchar(x) == 0L) | is.nan(x)
}

#' @export
omit_empty = function(x) x[!is_empty(x)]

#' Get Substrings by Indices and Separators
#'
#' For example \code{index_sep('a_b_c', 2, '_')} equals \code{'b'}.
#' For example \code{index_sep('a_b_c', c(1, 3), '_')} equals \code{'a_c'}.
#' For example \code{index_sep('a_b_c', -2, '_')} equals \code{'a_c'}.
#' For example \code{index_sep('a_b_c', 4, '_')} equals \code{''}.
#' For example \code{index_sep('a', 1, '_')} equals \code{'a'}.
#' For example \code{index_sep('a', 2, '_')} equals \code{''}.
#' For example \code{index_sep(c('a_b', 'c'), 2, '_')} equals \code{c('b', '')}.
#' For example \code{index_sep('a_b_c', c(3, 1), '_')} equals \code{'c_a'}.
#'
#' @param x character vector
#' @param i integer vector without sign mixing
#' @param sep length-one character vector
#' @param complement if \code{TRUE} the indices in \code{i} that are not matched are returned
#'
#' @export
index_sep = function(x, i, sep = "_") {
  complement = FALSE
  if (any(i < 0L)) {
    if (!all(i < 0L)) stop("cannot mix positive and negative indices")
    complement = TRUE
    i = -1 * i
  }
  stopifnot(length(sep) == 1L)
  if (complement) {
    n_separated_items = nchar(x) - nchar(gsub(sep, '', x)) + 1
    if (length(n_separated_items) > 1L) {
      stop('cannot use complement method with multiple inputs')
    }
    i = setdiff(seq_len(n_separated_items), i)
  }
  (x
   %>% as.character
   %>% strsplit(sep)
   %>% lapply(function(x) {
     ifelse(
       length(x) == 0L,
       '',
       paste0(omit_empty(x[i]), collapse = sep)
     )
   })
   %>% unlist
   %>% unname
  )
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

##' @export
initial_sim_report_names = function(model) {
  c(
    names(model$state),
    names(which_time_varying_rates(model)),
    names(model$sum_vector),
    names(model$factr_vector)
  )
}

##' @export
intermediate_sim_report_names = function(model) {
  c(
    initial_sim_report_names(model),
    names(model$sim_report_exprs)
  )
}

##' @export
final_sim_report_names = function(model) {
  c(
    intermediate_sim_report_names(model),
    unlist(lapply(model$lag_diff, getElement, "output_names")),
    unlist(lapply(model$conv, getElement, "output_names"))
  )
}

# flexmodel to latex (experimental) -------------------

##' @export
make_latex_symbols = function(nms) {
  greek_letters = c(
    "alpha",
    "beta",
    "gamma", "Gamma",
    "delta", "Delta",
    "epsilon", "varepsilon",
    "zeta",
    "eta",
    "theta", "vartheta", "Theta",
    "iota",
    "kappa",
    "lambda", "Lambda",
    "mu",
    "nu",
    "xi", "Xi",
    "pi", "Pi",
    "rho", "varrho",
    "sigma", "Sigma",
    "tau",
    "upsilon", "Upsilon",
    "phi", "varphi", "Phi",
    "chi",
    "psi", "Psi",
    "omega", "Omega"
  )

  # 1. split into two parts, the second of which becomes a subscript
  #    -- splitting mechanisms checked in this order
  #    a) first underscore
  #    b) letters followed by one or more numbers
  #    c) case change that is not in greek_letters or items in user-defined
  #       word_patterns
  #    e) do not split
  #
  # 2. figure out if the first part is a greek letter, and if so escape for use
  #    with latex
  #

  pats = list(
    lu = '^([a-z]+)([A-Z]+)(.*)$',
    ul = '^([A-Z]+)([a-z]+)(.*)$',
    nm = '^([A-z]+)([0-9]+)(.*)$')
  nms_first = nms
  nms_second = rep('', length(nms))
  for(i in seq_along(pats)) {
    l = grepl(pats[[i]], nms, perl = TRUE)
    nms_first[l] = sub(pats[[i]], '\\1', nms[l], perl = TRUE)
    nms_second[l] = sub(pats[[i]], '\\2\\3', nms[l], perl = TRUE)
  }
  l = grepl('^(\\w+)_', nms)
  nms_first[l] = index_sep(nms[l], 1)
  nms_second[l] = index_sep(nms[l], -1)

  is_greek = nms_first %in% greek_letters
  nms_first[is_greek] = "\\" %+% nms_first[is_greek]

  is_symb_word = (nchar(nms_first) > 1L) & !is_greek
  is_subscr_word = nchar(nms_second) > 1L
  nms_first[is_symb_word] = "\\text{" %+% nms_first[is_symb_word] %+% "}"
  nms_second[is_subscr_word] = "\\text{" %+% nms_second[is_subscr_word] %+% "}"

  latex_vars = nms_first
  empty_second = is_empty(nms_second)
  latex_vars[!empty_second] = nms_first[!empty_second] %_% nms_second[!empty_second]
  latex_vars
}


#' @export
make_latex_rates = function(model) {
  (get_rate_factors(model)
   %>% bind_rows(.id = "rate")
   %>% separate(rate, c("from", "to"), '_to_') # fragile
   %>% mutate(latex_var = make_latex_symbols(var))
   %>% mutate(latex_var = ifelse(compl, "1-" %+% latex_var, latex_var))
   %>% mutate(latex_var = ifelse(invrs, "1/" %+% latex_var, latex_var))
   %>% mutate(latex_var = "(" %+% latex_var %+% ")")
   %>% group_by(latex_var, from, to)
   %>% mutate(multiplicity = n())
   %>% ungroup()
   %>% group_by(from, to)
   %>% mutate(is_common_fact = (max(prod_indx) != 1L) & (multiplicity == max(prod_indx)))
   %>% ungroup()
   %>% group_by(prod_indx, from, to)
   %>% summarise(
     latex = paste0(latex_var[!is_common_fact], collapse = ""),
     latex_common_fact = paste0(latex_var[is_common_fact], collapse = "")
   )
   %>% ungroup
   %>% group_by(from, to)
   %>% summarise(
     left = ifelse(is_empty(latex_common_fact[1]), "", "\\left("),
     right = ifelse(is_empty(latex_common_fact[1]), "", "\\right)"),
     latex = latex_common_fact[1] %+% left %+% paste0(latex, collapse = "+") %+% right
   )
   %>% ungroup
   %>% select(-left, -right)
  )
}

#' @export
make_latex_flows = function(model) {
  rates = rate_summary(model, include_formula = TRUE, include_latex = TRUE)
  inflow = (rates
            %>% group_by(to)
            %>% summarise(flow = '+' %+% paste0(latex, collapse = "+"))
            %>% ungroup
            %>% rename(state = to)
  )
  outflow = (rates
             %>% group_by(from)
             %>% summarise(flow = '-\\left(' %+% paste0(latex, collapse = "+") %+% '\\right)')
             %>% ungroup
             %>% rename(state = from)
  )
  flows = (bind_rows(inflow, outflow)
           %>% group_by(state)
           %>% summarise(flow = paste0(flow, collapse = ''))
           %>% ungroup
  )
  latex_states = make_latex_symbols(flows$state)
  latex_states %+% "(t+1) = " %+% latex_states %+% "(t)" %+% flows$flow
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

#' @export
const_named_vector = function(nms, cnst) {
  setNames(rep(cnst[[1]], length(nms)), nms)
}

#' Merge One Vector into Another by Name
#'
#' If an item in \code{u} has the same name as an item
#' in \code{v} then replace the value in \code{v} with that
#' in \code{u}, otherwise create a new element in \code{v}.
#'
#' @param v named vector or list
#' @param u named vector of list
#' @export
merge_named_vectors = function(v, u) {
  if (is.null(names(v)) | is.null(names(u)) ) {
    stop("v and u must be named vectors")
  }
  for (nm in names(u)) {
    v[nm] = u[nm]
  }
  return(v)
}

#' @export
merge_named_vec_attr = function(v, u, a) {
  attr_v = attributes(v)
  attr_u = attributes(u)
  attr_v[[a]] = merge_named_vectors(attr_v[[a]], attr_u[[a]])
  attributes(v) = attr_v
  return(v)
}

ff = function(x) {
  if (is.atomic(x)) {
    return(x)
  }
  n = length(x) # number of hyperparameters
  d = lapply(x, dim) %>% unique
  if (length(d) != 1L) {
    stop('inconsistent hyperparameter dimensions')
  }
  d = unlist(d)
  if (!is.null(d)) {
    stop('cannot handle matrix-valued hyperparameters yet')
  }
  l = lapply(x, length) %>% unique
  if (length(l) != 1L) {
    stop('inconsistent hyperparameter lengths')
  }
  l = unlist(l)
  ii = outer(seq(from = 1, by = l, len = n), seq_len(l) - 1, `+`) %>% c
  unlist(x)[ii]
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
rate_summary = function(model, include_formula = FALSE, include_latex = FALSE) {
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

  if(include_latex) {
    summary = inner_join(
      summary,
      make_latex_rates(model),
      by = c("from", "to")
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



# utilities for parsing opt_params formulas ------------------

#' Parse Fitted Parameter Formula
#'
#' \code{trans1_parameter ~ trans2_prior(hyperparameter1, ..., laplace = FALSE)}
#'
#' If \code{parameter}
#'
#' @param x opt_params formula
#' @param params parameter vector
#' @param regex_param_nm should the parameter name be treated as a regular
#' expression used to produce a vector of parameter names
#' @export
parse_opt_param = function(x, params, params_timevar = NULL) {

  # TODO: allow for transformations associated with the regularization function
  # TODO: allow for vector-valued hyperparameters to be used with time variation
  stopifnot(is.call(x))
  stopifnot(as.character(x)[[1]] == '~')
  stopifnot(length(x) == 3L)

  # string parameter expressions are treated as regex patterns
  regex_param_nm = is.character(x[[2]])

  # non-regex patterns are split on plus
  if (!regex_param_nm) {
    param_nm = (x[2]
      %>% as.character
      %>% strsplit("+", fixed = TRUE)
      %>% lapply(trimws)
      %>% unlist
    )
  } else {
    param_nm = x[[2]]
  }

  # separate parameter name from transformation
  trans_nm = sapply(param_nm, index_sep,  1L, USE.NAMES = FALSE)
  for (i in seq_along(trans_nm)) {
    if (trans_nm[i] %in% valid_trans) {
      param_nm[i] = sapply(param_nm[i], index_sep, -1L, USE.NAMES = FALSE)
    } else {
      trans_nm[i] = ""
    }
  }

  # expand param_nm into a vector of parameter names that match
  if (regex_param_nm) {
    param_nm = grep(param_nm, names(params), value = TRUE, perl = TRUE)
  }

  # lengths of the parameter vector, p, and the time series, l
  p = length(param_nm)
  l = 0
  if (!is.null(params_timevar)) {
    l = sapply(param_nm, length_tv_mult_series)
    if (length(unique(l)) != 1L) {
      stop("all time varying parameters must have the same series length ",
           "if they are going to get the same prior")
    }
    l = l[1]
  }

  # cases:
  #  p = 1, l = 0: baseline prior on one parameter (scalar hyperparameters)
  #  p = 1, l > 0: tv_mult prior or priors on one parameter (scalar or l-vector hyperparameters)
  #  p > 1, l = 0: baseline priors on several parameters (scalar or p-vector hyperparameters)
  #  p > 1, l > 0: tv_mult priors on one parameter (scalar or l-by-p matrix hyperparameters)

  stopifnot(all(param_nm %in% names(params)))

  # re = "^" %+% alt_group(valid_trans %_% '') %+% "?" %+% alt_group(names(params)) %+% "$"
  # if (!grepl(re, param_nm)) {
  #   stop('unrecognized transformation or parameter is not in the model')
  # }
  # trans_nm = sub(re, '\\1', param_nm, perl = TRUE) %>% sub(pattern = '_$', replacement = '', perl = TRUE)
  # param_nm = sub(re, '\\2', param_nm, perl = TRUE)

  # extract prior and starting value information
  param_spec = x[[3]]
  prior_family = ''
  if (is.call(param_spec)) {
    rhs_func = as.character(param_spec[[1]])
    prior_trans = index_sep(rhs_func, 1)
    if (prior_trans %in% valid_trans) {
      prior_family = index_sep(rhs_func, -1)
    } else {
      prior_trans = ''
      prior_family = rhs_func
    }
    if (prior_family %in% valid_prior_families) {
      reg_params = vector("list", p * l)
      dim(reg_params) = c(l, p)
      for (i in l) {
        for(j in p) {
          reg_params[[i, j]] = try(eval(param_spec[[i+1]]))
        }
      }
      #reg_params = vector("list", length(param_spec) - 1)
      for (i in seq_along(reg_params)) {
        reg_params[[i]] = try(eval(param_spec[[i+1]]))
        reg_params[[i]] = set_hyperparam_dims(reg_params[[i]], l, p, param_nm)
      }
    } else {
      # TODO: worry about eval environments??
      reg_params = try(eval(param_spec))
      if (!is.numeric(reg_params)) {
        stop('unrecognized parameter specification')
      }
      reg_params = set_hyperparam_dims(reg_params, l, p, param_nm)
    }
  } else {
    if (!is.numeric(param_spec)) {
      stop('unrecognized parameter specification')
    }
    reg_params = set_hyperparam_dims(param_spec, l, p, param_nm)
    prior_family = ''
    prior_trans = ''
  }
  nlist(param_nm, trans_nm, prior_family, prior_trans, reg_params)

}

length_tv_mult_series = function(x) {
  sum(is.na(params_timevar$Value) & (params_timevar$Symbol == x))
}

set_hyperparam_dims = function(x, l, p, param_nm) {
  if (length(x) == 1L) {
    if (l > 0 & p > 1) {
      x = matrix(x, l, p)
    } else if (l > 0) {
      x = rep(x, l)
    } else if (p > 1) {
      x = rep(x, p)
    }
  } else {
    if (l > 0 & p > 1) {
      if (!all(dim(x) == c(l, p))) {
        stop("incompatible hyperparameter dimensions")
      }
    } else if (l > 0) {
      if (length(x) != l) {
        stop("incompatible numbers of hyperparameters")
      }
    } else if (p > 1) {
      if (length(x) != p) {
        stop("incompatible numbers of hyperparameters")
      }
    } else {
      stop("incompatible numbers of hyperparameters")
    }
  }
  if (is.matrix(x)) {
    colnames(x) = param_nm
  }
  x
}

parse_and_resolve_opt_form = function(x, params) {
  pf = parse_opt_form(x)
  pf$param = resolve_param(pf$param, params)
  pf
}

tmb_opt_form = function(pf, params, params_timevar = NULL) {
  if (is.null(params_timevar)) {
    if (!all(pf$param$param_nms %in% names(params))) {
      stop('parameters declared for optimization are missing from the model')
    }
  } else {
    if (!all(pf$param$param_nms %in% filter(params_timevar, is.na(Value))$Symbol)) {
      stop('time varying multipliers declared for optimization are not scheduled to vary in simulation time')
    }
  }
  fd = get_form_dim(pf, params, params_timevar)
  hyperparams_vec = (pf$prior$reg_param
    %>% lapply(hyperparam_to_vec, l = fd$l, p = fd$p)
    %>% unlist
    %>% matrix(nrow = fd$h, ncol = max(fd$l, 1) * fd$p, byrow = TRUE)
    %>% c
  )
  init_trans_params = c(pf$prior$reg_param[[1]])
  param_nms = rep(pf$param$param_nms, each = max(fd$l, 1))
  param_trans = rep(pf$param$trans, each = max(fd$l, 1))
  prior_distr = pf$prior$distr
  d = data.frame(param_nms, param_trans, prior_distr = pf$prior$distr, prior_trans = pf$prior$trans, count_hyperparams = fd$h, init_trans_params)
  d$prior_trans_id = find_vec_indices(d$prior_trans, c('', valid_trans))
  d$param_trans_id = find_vec_indices(d$param_trans, c('', valid_trans))
  d$prior_distr_id = find_vec_indices(d$prior_distr, valid_prior_families)
  if (is.null(params_timevar)) {
    d$opt_param_id = find_vec_indices(d$param_nms, params)
  } else {
    lookup_tv_table = (params_timevar
      %>% mutate(v = ifelse(is.na(Value), Symbol, ''))
    )
    lookup_tv_vec = lookup_tv_table$v
    d$tv_breaks = filter(lookup_tv_table, is.na(Value) & (Symbol == param_nms[1]))$breaks
    d$opt_tv_mult_id =
      unlist(lapply(
        unique(param_nms),
        find_vec_indices,
        lookup_tv_vec
      ))
  }
  nlist(d, hyperparams_vec)
}
get_hyperparam_list = function(reg_params,  n_hyperparams, n_params, n_breaks) {
  h = n_hyperparams
  p = n_params
  l = n_breaks
  ll = vector("list", max(1, l) * p)
  dim(ll) = c(max(1, l), p)
  for (k in seq_len(h)) {
    hyper_mat = matrix(reg_params[[k]], max(1, l), p)
    for (i in seq_len(max(1, l))) {
      for (j in seq_len(p)) {
        ll[[i, j]] = c(ll[[i, j]], hyper_mat[i, j])
      }
    }
  }
  ll
}
get_form_dim = function(x, params, params_timevar = NULL) {
  param_nm = x$param$param_nms
  # lengths of the parameter vector, p, and the time series, l
  h = x$prior$reg_params %>% length
  p = length(param_nm)
  l = 0
  if (!is.null(params_timevar)) {
    l = sapply(param_nm, length_tv_mult_series)
    if (length(unique(l)) != 1L) {
      stop("all time varying parameters must have the same series length ",
           "if they are going to get the same prior")
    }
    l = l[1]
  }
  l = unname(l)
  lapply(x$prior$reg_params, valid_dim_hyperparam, l, p)
  nlist(h, p, l)
}
dim_hyperparam = function(x) {
  if (is.null(dim(x))) {
    y = c(length(x), 1)
  } else {
    y = dim(x)
  }
  if (length(y) != 2) stop('only matrices, vectors, and scalar hyperparameters are allowed')
  y
}
valid_dim_hyperparam = function(x, l, p) {
  dh = dim_hyperparam(x)
  if (!any(
    isTRUE(all.equal(dh, c(l, p))),
    isTRUE(all.equal(dh, c(1, 1))),
    isTRUE(all.equal(dh, c(p, 1))) & (l == 0L),
    isTRUE(all.equal(dh, c(l, 1))) & (p == 1L)
  )) {
    stop('inconsistent hyperparameter dimensions')
  }
  NULL
}
hyperparam_to_vec = function(hyperparam, l, p) {
  hyperparam = (hyperparam
                %>% matrix(max(1, l), p)
                %>% t
                %>% c
  )
  hyperparam
}
resolve_param = function(x, params) {
  if (inherits(x, 'param_regex')) {
    param_nms = grep(x$regex, names(params), perl = TRUE, value = TRUE)
    trans = rep(x$trans, length(param_nms))
  } else {
    param_nms = x
    if (inherits(param_nms, 'param_vector')){
      param_nms = unclass(param_nms)
    }
    if (!is.character(param_nms)) {
      stop('unable to resolve parameter optimization specification')
    }
    trans = index_sep(param_nms, 1)
    which_have_trans = trans %in% valid_trans
    trans = ifelse(which_have_trans, trans, '')

    param_nms = ifelse(which_have_trans, Vectorize(index_sep)(param_nms, -1), param_nms)
  }
  nlist(param_nms, trans)
}
#' Recursive function to parse a formula containing sums of names
parse_name_sum = function(x) {
  y = character(0L)
  if (is.call(x)) {
    if (as.character(x[[1]]) != '+') {
      stop('"+" is the only operator/function allowed when summing, but "',
           as.character(x[[1]]), '" was used')
    }
    y = c(parse_name_sum(x[[2]]), parse_name_sum(x[[3]]))
    return(y)
  } else if (is.name(x)) {
    return(as.character(x))
  } else {
    stop('parameters must be expressed as valid R names (e.g. not character strings) ',
         'when summing parameters on the left-hand-side. if you are trying to ',
         'programmatically create formulas, you should instead programmatically ',
         'create the results of parsing the formula')
  }
}
#' Recursive function to parse a formula containing parameter
#' optimization specifications
#'
#' trans1_param ~ trans2_prior(hyperparameters, ...)
#'
#' trans1_param can be:
#'   (1) a symbol/name (or sum of symbols)
#'   (2) an expression that evaluates (in the environment of the formula)
#'       to a character vector
#'   (3) a literal length-1 character vector
#'
#' if trans1_param is a symbol:
#' trans1 = name of a valid parameter transformation, giving the scale
#'          on which the optimizer treats the parameter(s)
#' param = name of a parameter in the parameter vector
#'
#' if trans1_param is an expression
#'
parse_opt_form = function(x, e = NULL) {
  if (is.call(x)) {
    func = parse_opt_form(x[[1]], e)
    if (func == '~') {
      if (!is.null(e)) {
        # note that this message will display wrongly if someone
        # calls this function directly with an environment,
        # but this should be ok given that we are not going to
        # export this function
        stop('one may not use more than one tilde in optimization formulas')
      }
      # capture the environment of the formula so that it can
      # be recursively passed down, and used when/if eval is called
      e = environment(x)
      if (is.character(x[[2]])) {
        trans = index_sep(x[[2]], 1)
        if (trans %in% valid_trans) {
          regex = index_sep(x[[2]], -1)
        } else {
          regex = x[[2]]
          trans = ''
        }
        param = structure(nlist(regex, trans), class = 'param_regex')
      } else {
        param = parse_opt_form(x[[2]], e)
        if (!inherits(param, 'param_sum')) {
          param = structure(param, class = 'param_vector')
        }
      }
      prior = parse_opt_form(x[[3]], e)
      if (is.numeric(prior)) {
        stop('need to specify a prior distribution or initial value for the optimizer')
      }
      if (inherits(prior, 'param_sum')) {
        stop('summing prior distributions is not allowed')
      }
      # TODO: check to make sure that param_trans and prior_trans
      # are identical, and throw an informative error message
      # describing the missing feature of putting a prior over a
      # different scale than the parameter in the objective function
      return(nlist(param, prior))
    } else if (func == '+') {
      return(structure(parse_name_sum(x), class = 'param_vector'))
    } else {
      trans = index_sep(func, 1)
      if (trans %in% valid_trans) {
        distr = index_sep(func, -1)
      } else {
        distr = func
        trans = ''
      }
      if (distr %in% valid_prior_families) {
        reg_params = lapply(x[-1], parse_opt_form, e)
        return(structure(nlist(distr, trans, reg_params), class = 'prior'))
      } else {
        # if the function of the call is not a valid prior distribution,
        # assume that we just want to evaluate the call and do so in the
        # environment of the formula -- is this a good idea? -- could there
        # be issues if we decide to allow two tildes?
        return(parse_opt_form(eval(x, e), e))
      }
    }
  } else if (is.name(x)) {
    return(parse_opt_form(as.character(x), e))
  } else if (is.character(x)) {
    return(x)
  } else if (is.numeric(x)) {
    return(x)
  } else {
    stop('not a valid type')
  }
}

# computing indices and data for tmb ----------------------

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

#' @export
sim_report_expr_indices = function(exprs, init_sim_report_nms) {
  if (length(exprs) == 0L) {
    indices = list(
      sri_output = integer(0L),
      sr_count = integer(0L),
      sri = integer(0L),
      sr_modifier = integer(0L)
    )
    return(indices)
  }
  sri_output = seq_along(exprs) + length(init_sim_report_nms)
  sri <- get_var_indx(exprs)
  sr_count <- get_var_counts(exprs)
  sr_modifier <- get_var_modifiers(exprs)
  names(sri) = names(sr_count)  = names(sr_modifier) = NULL
  indices <- list(
    sri_output = sri_output,
    sr_count = sr_count,
    sri = sri,
    sr_modifier = sr_modifier
  )
  return(indices)
}

#' @export
lag_diff_indices = function(model) {
  indices_for_one_pattern = function(pattern_input) {
    nms = intermediate_sim_report_names(model)
    indices = grep(pattern_input$var_pattern, nms, perl = TRUE)
    data.frame(
      sri = indices,
      delay_n = rep(pattern_input$delay_n, length(indices))
    )
  }
  (model$lag_diff
    %>% lapply(indices_for_one_pattern)
    %>% Reduce(f = rbind)
    %>% as.list
  )
}

#' @export
conv_indices = function(model) {
  indices_for_one_pattern = function(pattern_input) {
    nms = intermediate_sim_report_names(model)
    indices = grep(pattern_input$var_pattern, nms, perl = TRUE)
    conv_par_indices = lapply(pattern_input$conv_pars, find_vec_indices, model$params)
    init_c_delay_cv = model$params[pattern_input$conv_pars$c_delay_cv]
    init_c_delay_mean = model$params[pattern_input$conv_pars$c_delay_mean]
    init_c_prop = model$params[pattern_input$conv_pars$c_prop]
    qmax = length(make_delay_kernel(init_c_prop, init_c_delay_mean, init_c_delay_cv)) + 1L
    #qmax = qgamma(
    #  0.95,
    #  1/init_c_delay_cv^2,
    #  init_c_delay_mean * init_c_delay_cv^2
    #)
    data.frame(
      sri = indices,
      c_prop_idx = rep(conv_par_indices$c_prop, length(indices)),
      c_delay_cv_idx = rep(conv_par_indices$c_delay_cv, length(indices)),
      c_delay_mean_idx = rep(conv_par_indices$c_delay_mean, length(indices)),
      qmax = rep(qmax, length(indices))
    )
  }
  (model$conv
    %>% lapply(indices_for_one_pattern)
    %>% Reduce(f = rbind)
    %>% as.list
  )
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
#' @section Motivating Example:
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

##' @export
tmb_observed_data = function(model) {

  initial_table = (model$observed$loss_params
   %>% mutate(loss_id = 1) # only negative binomial (id=1) currently
   %>% mutate(spi_loss_param = find_vec_indices(
     Parameter %_% Variable,
     c(model$state, model$params)
   ))
   %>% mutate(variable_id = find_vec_indices(
     Variable,
     model$condensation_map[final_sim_report_names(model)]
   ))
   %>% select(variable_id, loss_id, spi_loss_param)
  )
  variables_by_distributions = (initial_table
    %>% group_by(variable_id, loss_id)
    %>% summarise(loss_param_count = n())
    %>% ungroup
  )
  parameters = (initial_table
    %>% select(spi_loss_param)
  )
  comparisons = (model$observed$data
   %>% na.omit
   %>% rename(observed = value)
   %>% mutate(time_step = find_vec_indices(
     as.character(date),
     as.character(simulation_dates(model))
   ))
   %>% mutate(history_col_id = find_vec_indices(
     var,
     model$condensation_map[final_sim_report_names(model)]
   ))
   %>% select(time_step, history_col_id, observed)
   %>% arrange(time_step, history_col_id)
  )
  c(
    as.list(variables_by_distributions),
    as.list(parameters),
    as.list(comparisons)
  )
}

#' @export
tmb_opt_params = function(model) {
  indices = init_tmb_indices$opt_params

  if (length(model$opt_params) > 0L) {

    opt_tables = (model$opt_params
      %>% lapply(tmb_opt_form, model$params)
    )
    indices$index_table = (opt_tables
      %>% lapply(getElement, 'd')
      %>% do.call(what = 'rbind')
    )
    indices$hyperparameters = (opt_tables
      %>% lapply(getElement, 'hyperparams_vec')
      %>% unlist
    )
  }
  if (length(model$opt_tv_params) > 0L) {

    opt_tv_tables = (model$opt_tv_params
      %>% lapply(tmb_opt_form, model$params, model$timevar$piece_wise$schedule)
    )
    indices$index_tv_table = (opt_tv_tables
      %>% lapply(getElement, 'd')
      %>% do.call(what = 'rbind')
    )
    indices$hyperparameters_tv = (opt_tv_tables
      %>% lapply(getElement, 'hyperparams_vec')
      %>% unlist
    )
  }
  indices
}

# tmb_opt_params = function(model) {
#   return(tmb_opt_indices(model))
#
#   param_nms = lapply(model$opt_params, getElement, "param_nm")
#   params_per_formula = sapply(param_nms, length)
#   param_spi = lapply(param_nms, find_vec_indices, c(model$state, model$params))
#   trans_nms = lapply(model$opt_params, getElement, "trans_nm")
#
#   prior_families = lapply(model$opt_params, getElement, "prior_family")
#   reg_params = lapply(model$opt_params, getElement, "reg_params")
#
#   prior_families = rep(unlist(prior_families), params_per_formula)
#
#   count_reg_params = (model$opt_params
#     %>% lapply(getElement, "reg_params")
#     %>% sapply(length)
#   )
#   reg_params = (model$opt_params
#     %>% lapply(getElement, "reg_params")
#     %>% unlist
#   )
#   trans_id = find_vec_indices(trans_nms, c('', valid_trans))
#   prior_family_id = find_vec_indices(prior_families, c('', valid_prior_families))
#   nlist(param_spi, trans_id, count_reg_params, reg_params, prior_family_id)
# }


# retrieving information from tmb objective function --------------

#' Extract Parameter Vector to Pass to a TMB Function
#'
#' Get all of the parameters from a flexmodel in a vector that is
#' ready to be passed to a TMB AD objective function, gradient function,
#' hessian function, simulate function, report function.
#'
#' Currently this includes \code{params} and if \code{spec_ver_gt('0.1.0')}
#' it also includes time-varying multipliers in
#' \code{model$timevar$piece_wise$schedule$last_tv_mult}. If
#' \code{spec_ver_gt('0.1.2')}, only the parameters to be optimized are
#' returned.
#'
#'
#' @param model flexmodel
#' @export
tmb_params = function(model) {
  full_param_vec = c(unlist(model$params))
  if (spec_ver_gt('0.1.0')) {
    if(has_time_varying(model)) {
      last_tv_mult = model$timevar$piece_wise$schedule$last_tv_mult
      if (any(is.na(last_tv_mult)))
        stop('missing time-varying multipliers need to be replaced for use with TMB')
      full_param_vec = c(full_param_vec, last_tv_mult)
    }
  }
  if (spec_ver_gt('0.1.2')) {
    if (exists_opt_params(model)) {
      full_param_vec = setNames(
        full_param_vec[tmb_map_indices(model)],
        tmb_param_names(model)
      )
    }
  }
  full_param_vec
}

#' Parameter Vectors Transformed for TMB
#'
#' @params model flexmodel
#' @param vec_type type of vector to return:
#' (1) tmb_fun_arg = same length as the argument for TMB objective function,
#' (2) params = same length as the params element in \code{model}
#' (3) tv_mult = same length as the \code{model$timevar$piece_wise$schedule$init_tv_mult}
#'
#' @export
tmb_params_trans = function(model, vec_type = c('tmb_fun_arg', 'params', 'tv_mult')) {
  spec_check(
    introduced_version = '0.2.0',
    feature = 'parameter transformations'
  )
  vec_type = match.arg(vec_type)
  model = update_tmb_indices(model)
  if (vec_type != 'tmb_fun_arg') {
    opt_params = model$tmb_indices$opt_params
    table = switch(vec_type
      , params = opt_params$index_table
      , tv_mult = opt_params$index_tv_table
    )
    id = switch(vec_type
      , params = table$opt_param_id
      , tv_mult = table$opt_tv_mult_id
    )
    vec = switch(vec_type
      , params = model$params
      , tv_mult = model$timevar$piece_wise$schedule$init_tv_mult
    )
    vec[id] = table$init_trans_params
    return(vec)
  } else {
    init_trans_params = (model
       $  tmb_indices
       $  opt_params
       $  index_table
      %>% with(setNames(init_trans_params, param_nms))
    )
    init_trans_tv_mult = (model
       $  tmb_indices
       $  opt_params
       $  index_tv_table
      %>% getElement('init_trans_params')
    )
    return(c(init_trans_params, init_trans_tv_mult))
  }
}

#' @export
tmb_map_indices = function(model) {
  spec_check(
    introduced_version = '0.2.0',
    feature = 'flex models can use map argument of MakeADFun'
  )
  (model
    $  tmb_indices
    $  ad_fun_map
   %>% unlist
   %>% is.na
   %>% `!`
   %>% which
   %>% unname
  )
}

#' @export
tmb_param_names = function(model) {
  spec_check(
    introduced_version = '0.2.0',
    feature = 'flex models can use map argument of MakeADFun'
  )
  (model
     $  tmb_indices
     $  ad_fun_map
    %>% unlist
    %>% as.character
    %>% omit_empty
  )
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
  if (getOption("MP_force_full_outflow")) {
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

#' @export
simulation_history = function(model, add_dates = TRUE, sim_params = NULL) {
  if (is.null(sim_params)) sim_params = tmb_params(model)
  sim_hist = (tmb_fun(model)
    $  simulate(sim_params)
    $  simulation_history
   %>% as.data.frame
   %>% setNames(final_sim_report_names(model))
  )
  if (add_dates) {
    sim_hist = (model
      %>% simulation_dates
      %>% data.frame
      %>% setNames("Date")
      %>% cbind(sim_hist)
    )
  }
  sim_hist
}

#' Condensed set of Simulated Variables
#'
#' @export
simulation_condensed = function(model, add_dates = TRUE, sim_params = NULL) {
  cond_map = model$condensation_map
  cond_nms = names(cond_map)
  if (add_dates) {
    cond_nms = c("Date", cond_nms)
    cond_map = c(Date = "Date", cond_map)
  }
  sims = simulation_history(model, add_dates = TRUE)[cond_nms]
  names(sims) = c(cond_map)
  sims
}

#' @importFrom tidyr pivot_longer
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

#' Simulated Variables to Compare with Observed Data
#'
#' @export
simulation_fitted = function(model) {
  obsvars = unique(model$observed$data$var)
  simulation_condense(model)[obsvars]
}

# benchmarking and comparison in tests -----------------------

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
#' @importFrom testthat testthat_tolerance
#' @export
compare_sims = function(classic_sim, tmb_sim, tolerance = testthat_tolerance(), compare_attr = TRUE) {
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

#' Compare Automatic and Numerical Differentiation
#'
#' Compare differences in gradients of loss functions computed using
#' automatic and numerical differentiation
#'
#' @param model flexmodel
#' @param tolerance numerical tolerance to pass to \code{all.equal} --
#' if \code{NA}, then a list is returned that can be passed to
#' \code{all.equal} using \code{do.call}
#' @param ... additional arguments to pass to \code{numDeriv::grad}
#'
#' @return either (1) the return value of \code{all.equal} comparing
#' gradients using TMB auto differentiation to \code{numDeriv::grad} if
#' \code{is.numeric(tolerance)},
#' or (2) a list is returned that can be passed to
#' \code{all.equal} using \code{do.call} if \code{is.na(tolerance)}
#'
#' @importFrom numDeriv grad
#' @export
compare_grads = function(model, tolerance = 1e-5,  ...) {

  args = list(...)
  if (!'method.args' %in% names(args)) {
    args$method.args = list(
        eps = 1e-4, d = 0.1,
        zero.tol = sqrt(.Machine$double.eps / 7e-7), r = 4, v = 2,
        show.details = FALSE
      )
  }
  dd = tmb_fun(model)
  current <- numDeriv::grad(dd$fn, dd$par, ...)
  target <- dd$gr(dd$par)
  attributes(target) <- attributes(current) <- NULL
  if (is.na(tolerance)) return(nlist(target, current))
  all.equal(target, current, tolerance)
}

# spec versioning -----------------------

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
  print(spec_version)
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
