# objects ----------------------

valid_loss_params = list(
  negative_binomial = 'nb_disp',
  poisson = character(0L),
  log_normal = 'normal_sd',
  normal = 'normal_sd',
  beta = c('beta_shape1', 'beta_shape2')
)
valid_loss_functions = names(valid_loss_params)
valid_prior_families = c('flat', 'normal')
valid_trans = c('log', 'log10', 'logit', 'cloglog', 'inverse')
valid_tv_types = c(
  'abs',
  'rel_orig', 'rel_prev',
  'rel_orig_logit', 'rel_prev_logit'
)

filter_loss_params = function(loss_params, param_nms, dist_nm){
  filter(
    loss_params,
    Parameter == param_nms,
    Distribution == dist_nm
  )$Variable
}

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
init_pow_vector = numeric(0L)
init_pow = setNames(
  as.data.frame(matrix(nrow = 0L, ncol = 5)),
  c("pow_nms", "pow_arg1_nms", "pow_arg2_nms", "pow_const_nms", "initial_value")
)

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
    Distribution = character(),
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
  pow_indices = list(
    powidx = integer(0L),
    powarg1idx = integer(0L),
    powarg2idx = integer(0L),
    powconstidx = integer(0L)
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
      tv_breaks = integer(0),
      opt_tv_mult_id = integer(0)
    ),
    hyperparameters_tv = numeric(0L)
  )
)




# test functions ------------------------------------------------

check_len_bounds = function(min_length, max_length) {
  stopifnot(is.integer(min_length) & (is.integer(max_length) | is.infinite(max_length)) & min_length >= 0L & max_length >= min_length)
}

check_len = function(x) {
  len = length(x$default)
  check_len_bounds(x$min_len, len)
  check_len_bounds(len, x$max_len)
}

bundle_def = function(..., class) {
  structure(nlist(...), class = class)
}

def_vectors = c('logical', 'integer', 'double', 'character', 'list')
def_numerics = c('integer', 'double')
def_named_containers = c('nlist', 'dataframe')

make_def_function = function(def_name) {
  if (def_name %in% def_named_containers) {
    f = function(...) {
      l = list(...)
      stopifnot(!is.null(names(l)))
      stopifnot(!any(duplicated(names(l))))
      structure(l, class = 'def' %_% def_name)
    }
    return(f)
  }
  f = function() {
    l = as.list(environment())
    if (!l$optional & l$extra) {
      stop('extra components must be optional')
    }
    if (isFALSE(l$allow_dups)) {
      if (any(duplicated(l$default))) {
        stop('must either allow dups or supply an appropriate default if minimum length is greater than one')
      }
    }
    structure(l, class = 'def' %_% def_name)
  }
  if (def_name %in% def_vectors) {
    formals(f) = c(formals(f), alist(min_len = 0L, max_len = 1L, default = vector(def_name, min_len)))
    if (def_name != 'logical') {
      formals(f) = c(formals(f), alist(allow_dups = TRUE))
    }
  }
  if (def_name %in% def_numerics) {
    formals(f) = c(formals(f), alist(min_val = -Inf, max_val = Inf))
  }
  if (def_name == 'character') {
    formals(f) = c(formals(f), alist(pattern = '^.*$'))
  } else if (def_name == 'formula') {
    formals(f) = c(formals(f), alist(n_sides = 'any')) # 'one', 'two'
  } else if (def_name == 'list') {
    formals(f) = c(formals(f), alist(component_def = make_def_function('character')()))
  }
  formals(f) = c(formals(f), alist(optional = FALSE, extra = FALSE))
  f
}

def = sapply(c(
  'logical',
  'integer',
  'double',
  'character',
  'formula',
  'struc',
  'list',
  'nlist',
  'dataframe'
), make_def_function, simplify = FALSE)

def$model = function(name, model_spec_version) {
  def$nlist(
    name = def$character(min_len = 1L, max_len = Inf, default = name),
    model_spec_version = def$character(
      min_len = 1L, max_len = 1L,
      pattern = '[0-9]+\\.[0-9]+\\.[0-9]+',
      default = model_spec_version
    ),
    states = def$dataframe(
      state = def$character(min_len = 1L, max_len = Inf, allow_dups = FALSE),
      description = def$character(allow_dups = FALSE, optional = TRUE),
      layers = def$character(allow_dups = TRUE, option = TRUE, extra = TRUE)
    ),
    parameters = def$dataframe(
      parameter = def$character(allow_dups = FALSE),
      description = def$character(allow_dups = FALSE, optional = TRUE),
      layers = def$character(option = TRUE, extra = TRUE)
    ),
    rates = def$dataframe(
      from = def$character(),
      to = def$character(),
      expression = def$character(),
      outflow = def$logical()
    ),
    summaries = def$dataframe(
      summary = def$character(allow_dups = FALSE),
      type = def$character(pattern = "^(sum|intermediate|lag_n_diff|conv)$"),
      expression = def$character()
    ),
    structure = def$list(
      allow_dups = FALSE,
      component_def = def$struc()
    )
  )
}

check_def = function(x, ...) {
  UseMethod("check_def")
}

check_def.def_logical = function(x, ...) {
  check_len(x)
}

init_def = function(x, ...) {
  if (isTRUE(x$extra)) return(NULL)
  UseMethod("init_def")
}

init_def.def_logical = function(x, ...) {
  x$default
}

init_def.def_integer = function(x, ...) {
  x$default
}

init_def.def_double = function(x, ...) {
  x$default
}

init_def.def_character = function(x, ...) {
  x$default
}

init_def.def_formula = function(x, ...) {
  if (x$n_sides == 'two') return(. ~ 1)
  ~ 1
}

init_def.def_struc = function(x, ...) {
  struc('1')
}

init_def.def_list = function(x, ...) {
  if (x$min_len == 0L) return(list())
  rep(list(init_def(x$component_def)), x$min_len)
}

init_def.def_nlist = function(x, ...) {
  sapply(x, init_def, simplify = FALSE)
}

init_def.def_dataframe = function(x, ...) {
  l = sapply(x, init_def, simplify = FALSE)
  l[sapply(l, is.null)] = NULL
  as.data.frame(l)
}

is_len1_char = function(x) {
  isTRUE((length(x) == 1L) & is.character(x))
}

assert_len1_char = function(x) {
  nm = deparse(substitute(x))
  if(!is_len1_char(x)) {
    stop(nm, ' must be a length-1 character')
  }
}

is_len1_int = function(x) {
  # not good enough to just check for integer, because
  # a double may effectively be an integer
  (length(x) == 1L) & isTRUE(all.equal(x, as.integer(x)))
}

assert_len1_int = function(x) {
  nm = deparse(substitute(x))
  if(!is_len1_int(x)) {
    stop(nm, ' must be a length-1 integer')
  }
}

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

is_fitted_by_bbmle = function(model) {
  x = inherits(model, "flexmodel_calibrated")
  if (!isTRUE(x)) return(x)
  x & isTRUE(attr(class(model$opt_obj), "package") == "bbmle")
}

# constructing names and strings ----------------------

is_empty = function(x) {
  x = as.character(x)
  is.na(x) | (nchar(x) == 0L) | is.nan(x)
}

omit_empty = function(x) x[!is_empty(x)]

as_vector_no_attr = function(x) {
  nms = names(x)
  attributes(x) = NULL
  setNames(x, nms)
}

as_data_frame_no_row_names = function(x) {
  rownames(x) = NULL
  x
}

initial_sim_report_names = function(model) {
  c(
    names(model$state),
    names(which_time_varying_rates(model)),
    names(model$sum_vector),
    names(model$factr_vector),
    names(model$pow_vector)
  )
}

intermediate_sim_report_names = function(model) {
  c(
    initial_sim_report_names(model),
    names(model$sim_report_exprs)
  )
}

final_sim_report_names = function(model) {
  if (spec_ver_gt('0.2.0')) {
    lags = model$lag_diff_uneven
  } else {
    lags = model$lag_diff
  }
  c(
    intermediate_sim_report_names(model),
    unlist(lapply(lags, getElement, "output_names")),
    #unlist(lapply(model$lag_diff_uneven, getElement, "output_names")),
    unlist(lapply(model$conv, getElement, "output_names"))
  )
}

condensed_sim_report_names = function(model) {
  nms = final_sim_report_names(model)
  if (model$no_condensation) return(nms)
  condense_names(nms, model$condensation_map)
}

pad_lag_diffs = function(sims, lag_diff) {
  if(length(lag_diff) == 0L) return(sims)
  ff = function(x) {
    as.data.frame(x[c('delay_n', 'output_names')])
  }
  d = bind_rows(lapply(lag_diff, ff))
  for(i in 1:nrow(d)) {
    sims[1:as.integer(d[i,'delay_n']), d[i,'output_names']] = NA
  }
  sims
}

pad_convs = function(sims, conv, conv_indices) {
  if(length(conv) == 0L) return
  conv_nms = unlist(lapply(conv, getElement, 'output_names'))
  qmax = conv_indices$qmax
  for(i in 1:length(conv_nms)) {
    sims[1:(qmax[i]-2L), conv_nms[i]] = NA
  }
  sims
}

# constructing data frames ----------------

#' Vector to Data Frame
#'
#' Convert a named vector to a data frame with two columns, one for the
#' names and one for the values.
#'
#' @param x named vector
#' @param values_col name of the values column -- the default takes a guess
#' @param names_col name of the names column
#'
#' @export
v2d = function(x, values_col = NULL, names_col = "names") {
  stopifnot(!is.null(names(x)))
  if (is.null(values_col)) {
    values_col = make.names(deparse(substitute(x)))
  }
  x = setNames(data.frame(names(x), x), c(names_col, values_col))
  row.names(x) = NULL
  x
}

# flexmodel to latex (experimental) -------------------

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

unlist_params = function(x) {
  # -- would be nice to use get_attr and put_attr from utils.R, but
  #    the latter will set things back to pansim and this is not exactly
  #    what we want
  pattr = attributes(x)
  x = setNames(unlist(x), names(x))  # FIXME: will silently fail for nested lists
  attributes(x) = c(attributes(x), pattr)
  x
}

update_full_condensation_map = function(model) {
  srn = final_sim_report_names(model)
  model$condensation_map = setNames(srn, srn)
  model
}

condense_names = function(nms, nm_map) {
  nm_map[nms[nms %in% names(nm_map)]]
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

null_to_char0 = function(x) {
  if(is.null(x)) return(character(0L))
  as.character(x)
}

null_to_int0 = function(x) {
  if(is.null(x)) return(integer(0L))
  as.integer(x)
}

null_to_num0 = function(x) {
  if(is.null(x)) return(numeric(0L))
  as.numeric(x)
}

null_to_log0 = function(x) {
  if(is.null(x)) return(logical(0L))
  as.logical(x)
}

null_to_charNA = function(x) {
  if(is.null(x)) return(as.character(NA))
  as.character(x)
}

null_to_intNA = function(x) {
  if(is.null(x)) return(as.integer(NA))
  as.integer(x)
}

null_to_numNA = function(x) {
  if(is.null(x)) return(as.numeric(NA))
  as.numeric(x)
}

null_to_logNA = function(x) {
  if(is.null(x)) return(as.logical(NA))
  as.logical(x)
}

null_to_0 = function(x) {
  if(is.null(x)) return(0L)
  as.integer(x)
}

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

block = function(from, to, mat) {
  stop('Blockwise rate specification is not implemented')
}

## not currently working
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

# @param x character vector of names to look for in \code{vec} or
# \code{names(vec)}
# @param vec character vector or object with names
find_vec_indices <- function(x, vec) {
  if(!is.character(vec)) vec = names(vec)
  missing_variables = x[!x %in% vec]
  if(length(missing_variables) > 0L) {
    stop("the following variables were used but not found in the model:\n",
         paste0(missing_variables, collapse = "\n"),
         "\nyou might find the help page for ?avail_for_rate useful")
  }
  (x
    %>% as.character()
    %>% outer(vec, "==")
    %>% apply(1, which)
  )
}

# Nested Indices
#
# Create and return a nested set of character vectors from a sequence
# of regular expressions, and return indices into each vector for recovering
# other shorter vectors that are lower in the hierarchy.
#
# @section Motivating Example:
#
# There are three nested state vectors
# a_states -- all states
# p_states -- excludes accumulators
# e_states -- includes only infected states
#
# There are three index vectors
# i_ap -- indexes a_states to yield p_states
# i_ae -- indexes a_states to yield e_states
# i_pe -- indexes p_states to yield e_states
#
# There are also three inverse index vectors
# j_ap -- indexes a_states to yield setdiff(a_states, p_states)
# j_ae -- indexes a_states to yield setdiff(a_states, e_states)
# j_pe -- indexes p_states to yield setdiff(p_states, e_states)
#
# In general -- assume that y_states are nested in x_states
# i_xy -- indexes x_states to yield y_states
# j_xy -- indexes x_states to yeild setdiff(x_states, y_states)
# i_yx -- doesn't exist because not all x_states are also y_states
#
# @param x character vector
# @param patterns named list or vector of regular expressions to be applied
# sequentially to create a nested set of character vectors
# @return List with three elements.
# \describe{
#   \item{\code{x}}{Nested character vectors.}
#   \item{\code{i}}{
#     List of lists, with inner list giving the indices into character vectors
#     that are higher in the hierarchy. For example, one may read this
#     expression, \code{x$e == x$p\\[i$e$p\\]}, as "e equals p at the index
#     that takes p to e".
#   }
#   \item{\code{j}}{Set difference version of \code{i}.}
# }

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
  is_greater_than = function(x, y) x > y
  repeated_transitions = (rates
   %>% names
   %>% table
   %>% is_greater_than(1L)
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

# regenerate rates in a model. this is useful if parameters get inserted
# into a model and therefore throw off the indices that are saved in the
# rates component of a model -- the utility of this function is actually
# a sign that we have spaghetti code and therefore should create a more
# formal class structure
regen_model = function(model) {
  (model
   %>% regen_rates
   %>% regen_sim_report_expr
  )
}

regen_rates = function(model) {
  remodel = model
  remodel$rates = list()
  for(r in seq_along(model$rates)) {
    args = c(list(model = remodel), model$rates[[r]][c('from', 'to', 'formula')])
    remodel = do.call(add_rate, args)
  }
  remodel
}

##' Rate Summary
##'
##' Summarize the properties of the rates of transition amongst the
##' compartments in the model.
##'
##' @param model a \code{\link{flexmodel}} object
##' @param include_formula include a column for the expanded rate formula
##' @export
rate_summary = function(model, include_parse_info = TRUE, include_formula = FALSE, include_latex = FALSE) {
  summary = data.frame(
    from = get_rate_info(model, "from") %>% unlist,
    to = get_rate_info(model, "to") %>% unlist
  )

  if(include_parse_info) {
    summary$n_fctrs = get_n_factors(model) %>% unlist
    summary$n_prdcts = get_n_products(model) %>% unlist
    summary$n_vrbls = get_n_variables(model) %>% unlist
    summary$state_dependent = get_rate_info(model, "state_dependent") %>% unlist
    summary$time_varying = get_rate_info(model, "time_varying") %>% unlist
    summary$sum_dependent = get_rate_info(model, "sum_dependent") %>% unlist
  }

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

regen_sim_report_expr = function(model) {
  remodel = model
  remodel$sim_report_exprs = list()
  for(r in seq_along(model$sim_report_exprs)) {
    args = c(list(model = remodel), model$sim_report_exprs[[r]][c('expr_nm', 'formula')])
    remodel = do.call(add_sim_report_expr, args)
  }
  remodel
}


## @param x parameter vector or flexmodel
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

time_varying_rates <- function(model) {
  model$rates[which_time_varying_rates(model)]
}

which_time_varying_rates <- function(model) {
  sd  <- get_rate_info(model, "state_dependent") %>% unlist
  tv  <- get_rate_info(model, "time_varying") %>% unlist
  smd <- get_rate_info(model, "sum_dependent") %>% unlist
  which(sd | tv | smd)
}

state_dependent_rates <- function(model) {
  i = get_rate_info(model, "state_dependent") %>% unlist
  model$rates[i]
}

sum_dependent_rates = function(model) {
  i = get_rate_info(model, "sum_dependent") %>% unlist
  model$rates[i]
}

compute_rates = function(model) {
  (model
   %>% get_rate_info("formula")
   %>% eval_formulas(get_var_list(model))
   %>% setNames(names(model$rates))
  )
}

compute_factrs = function(model) {
  (model
   %>% get_factr_info("formula")
   %>% eval_formulas(get_var_list(model))
   %>% setNames(names(model$factrs))
  )
}

eval_formulas = function(formula_list, var_list) {
  (formula_list
   %>% lapply(function(x) ifelse(inherits(x, 'formula'), as.character(x[2]), x))
   %>% sapply(as.character)
   %>% struc
   %>% struc_eval(var_list)
   %>% c()
  )
}

get_var_list = function(model) {
  c(
    as.list(model$params),
    as.list(model$state),
    as.list(model$sum_vector),
    as.list(model$factr_vector),
    as.list(model$pow_vector)
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

lookup_pairwise = function(from, to, M) {
  i = pwise(from, to, M)
  data.frame(
    from = rownames(M)[i[,"from_pos"]],
    to = colnames(M)[i[,"to_pos"]]
  )
}


# Rate Matrix Lookup Table
#
# @param state state_pansim object
# @param ratemat rate matrix
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

# utilities for parsing opt_params and loss_params formulas and structures ------------------

remove_opt_param = function(model, focal_param) {
  opt_param_is_focal = (model$opt_params
    %>% lapply(getElement, 'param')
    %>% lapply(getElement, 'param_nms')
    %>% lapply(`==`, focal_param)
  )
  opt_params = model$opt_params[!unlist(lapply(opt_param_is_focal, all))]
  model$opt_params = mapply(function(x, y) {
      x$param$param_nms = x$param$param_nms[!y]
      x
    },
    opt_params,
    opt_param_is_focal,
    SIMPLIFY = FALSE
  )
  model
}

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
    l = sapply(param_nm, length_tv_mult_series, params_timevar)
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

length_tv_mult_series = function(x, params_timevar, tv_type) {
  sum(is.na(params_timevar$Value) & (params_timevar$Symbol == x) & (params_timevar$Type == tv_type))
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
  if (length(unique(c(pf$param$trans, pf$prior$trans))) != 1L) {
    example_bad = (pf$param$trans[1]
        %_% pf$param$param_nms[1]
        %~% (pf$prior$trans[1] %_% pf$prior$distr[1])
        %+% "(" %+% paste0(pf$prior$reg_params, collapse = ', ') %+% ")"
    )
    example_good = (pf$param$trans[1]
        %_% pf$param$param_nms[1]
        %~% (pf$param$trans[1] %_% pf$prior$distr[1])
        %+% "(" %+% paste0(pf$prior$reg_params, collapse = ', ') %+% ")"
    )
    stop(
      '\nthe transformation scale used when a parameter is an argument',
      '\nof an objective function, must be the same transformation used',
      '\nwhen it is an argument of a prior density function.',
      '\nfor example one may not optimize the transmission rate on the',
      '\nlog scale but take the normal prior density over the identity scale.',
      '\nthis restriction might be lifted in the future, which is the reason',
      '\nwhy the notation makes it possible to decouple these two types of',
      '\ntransformations.',
      '\n\nproblematic example inspired by your input:',
      '\n', example_bad,
      '\n\nvalid example inspired by your input:',
      '\n', example_good
    )
  }
  pf
}

parse_and_resolve_loss_form = function(x, hist_nms, param_nms) {
  pf = parse_loss_form(x)
  if (!all(unique(pf$Parameter) %in% param_nms)) {
    if (!all(is.na(pf$Parameter))) {
      stop(
        'parameters used in an error distribution are missing from the ',
        'model parameter vector'
      )
    } else if (length(pf$Parameter) != 1L) {
      stop(
        'error distribution without parameters is specified with parameters'
      )
    }
  }
  if (!all(unique(pf$Variable) %in% hist_nms)) {
    stop(
      'variables declared for observation error are missing from the ',
      'model simulation history'
    )
  }
  pf
}

tmb_opt_form = function(pf, params, params_timevar = NULL, tv_type = NULL) {
  if (is.null(params_timevar)) {
    stopifnot(is.null(tv_type))
    if (!all(pf$param$param_nms %in% names(params))) {
      stop('parameters declared for optimization are missing from the model')
    }
  } else {
    stopifnot(!is.null(tv_type))
    if (!all(pf$param$param_nms %in% filter(params_timevar, is.na(Value))$Symbol)) {
      stop('time varying multipliers declared for optimization are not scheduled to vary in simulation time')
    }
  }
  fd = get_form_dim(pf, params, params_timevar, tv_type)
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
    lookup_tv_vec = filter(lookup_tv_table, Type == tv_type)$v
    d$tv_breaks = filter(
      lookup_tv_table,
      is.na(Value) & (Symbol == param_nms[1]) & (Type == tv_type)
    )$breaks
    # d$tv_breaks = filter(lookup_tv_table, Symbol == param_nms[1])$breaks
    #(lookup_tv_table
    #  %>% mutate()
    #)
    d$opt_tv_mult_id = filter(
      lookup_tv_table,
      is.na(Value),
      Symbol %in% param_nms,
      Type == tv_type
    )$tv_mult_id
      # unlist(lapply(
      #   unique(param_nms),
      #   find_vec_indices,
      #   lookup_tv_vec
      # ))
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
get_form_dim = function(x, params, params_timevar = NULL, tv_type = NULL) {
  param_nm = x$param$param_nms
  # lengths of the parameter vector, p, and the time series, l
  h = x$prior$reg_params %>% length
  p = length(param_nm)
  l = 0
  if (!is.null(params_timevar)) {
    l = sapply(param_nm, length_tv_mult_series, params_timevar, tv_type)
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
  if (!all(param_nms %in% names(params))) {
    stop('cannot find opimization parameters in the params vector')
  }
  nlist(param_nms, trans)
}

# Recursive function to parse a formula containing sums of names
parse_name_sum = function(x) {
  stop("sums in opt_param formulas is under construction. ",
       "please specify the prior for each parameter separately.")
  y = character(0L)
  if (is.call(x)) {
    if (as.character(x[[1]]) != '+') {
      stop('you cannot apply functions or operators to terms in a sum, ',
           'but ', as.character(x[[1]]), ' was used')
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

parse_loss_form = function(x, e = NULL) {
  if (is.call(x)) {
    func = parse_loss_form(x[[1]], e)
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
        history_variable = parse_loss_form(x[[2]], e)
        # if (!inherits(history_variable, 'param_sum')) {
        #   history_variable = structure(history_variable, class = 'history_vector')
        # }
      }
      loss_dist = parse_loss_form(x[[3]], e)
      if (is.numeric(loss_dist)) {
        stop('need to specify a prior distribution or initial value for the optimizer')
      }
      loss_dist$Variable = history_variable
      return(loss_dist)
    } else if (func == '+') {
      stop('loss formulas cannot contain sums')
      # return(structure(parse_name_sum(x), class = 'param_vector'))
    } else {
      # at this point we should know that the function is a loss function
      loss_func = func
      if (loss_func %in% valid_loss_functions) {
        loss_params = lapply(x[-1], parse_loss_form, e)

        # check for the right number of loss parameters,
        # given the loss function
        if (length(loss_params) != length(valid_loss_params[[loss_func]])) {
          stop(
            length(loss_params), ' loss parameters were ',
            'specified, but exactly ', length(valid_loss_params[[loss_func]]),
            ' needs to be provided'
          )
        }

        # check that the loss parameters are specified as length-1
        # character vectors
        all_len1_char = all(unlist(lapply(loss_params, is_len1_char)))
        if (!all_len1_char) {
          stop(
            'loss parameters must be length-1 character vectors, ',
            'with names that will be added to the param vector'
          )
        }

        # check that the parameter names are syntactically valid
        name_re = wrap_exact(getOption("MP_name_search_regex"))
        all_valid_param_nms = all(
          unlist(lapply(loss_params, grepl, pattern = name_re))
        )
        if (!all_valid_param_nms) {
          stop('not all loss parameter names are syntactically valid')
        }

        unlist_loss_params = unlist(loss_params)
        if (is.null(unlist_loss_params)) unlist_loss_params = NA
        loss_params_data = data.frame(
          Parameter = unlist_loss_params,
          Distribution = loss_func
        )
        return(loss_params_data)
      } else {
        stop(
          '\n', loss_func,
          ' is not one of the following valid loss functions:\n',
          paste0(valid_loss_functions, collapse = ', ')
        )
      }
    }
  } else if (is.name(x)) {
    return(parse_loss_form(as.character(x), e))
  } else if (is.character(x)) {
    return(x)
  } else if (is.numeric(x)) {
    return(x)
  } else {
    stop('not a valid type')
  }
}

# date and time utilities -----------

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
  model = update_tmb_indices(model)
  if (inherits(model, "flexmodel_to_calibrate")) {
    p = tmb_params_trans(model)
    if (length(p) == 0L) {
      stop("observed data have been supplied, but optimization parameters have not been specified")
    }
    return(p)
  }
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
#' TODO: this should not be user-facing
#'
#' @param model flexmodel
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
    init_trans_params_ = with(model$tmb_indices$opt_params$index_table, {
      setNames(
        init_trans_params,
        param_trans %_% param_nms
      )
    })
    itvt = model$tmb_indices$opt_params$index_tv_table
    if (nrow(itvt) == 0L) return(c(init_trans_params_))
    init_trans_tv_mult_ = with(itvt, {
      setNames(
        init_trans_params,
        param_trans %_% param_nms %_% "t" %+% tv_breaks
      )
    })
    return(c(init_trans_params_, init_trans_tv_mult_))
  }
}

# FIXME: not used
tmb_params_init = function(model) {
  model = update_tmb_indices(model)
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
  c(init_trans_params, init_trans_tv_mult)
}

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

tmb_param_names = function(model) {
  stop("not working")
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

simulation_dates = function(model) {
  seq.Date(
    as.Date(model$start_date),
    as.Date(model$end_date),
    by = 1)
}

compute_num_iters = function(model) {
  (model$end_date
    %>% difftime(model$start_date, units = 'days')
    %>% as.integer
  )
}

update_params_calibrated = function(model) {
  # TODO: check if opt_par exists
  if (!inherits(model, 'flexmodel_calibrated')) {
    stop('this is not a flexmodel_calibrated object. ',
         'please first use nlminb_flexmodel or optim_flexmodel.')
  }
  obj_fun = tmb_fun(model)
  report = obj_fun$report(model$opt_par)
  params_calibrated = model$params
  params_calibrated[] = report$params
  params_calibrated_timevar = (model$timevar$piece_wise$schedule
    %>% select(Date, Symbol, Value, Type)
    %>% within(Value[is.na(Value)] <- report$tv_mult[is.na(Value)])
  )
  model$params = params_calibrated
  model = update_piece_wise(model, params_calibrated_timevar)
  model$opt_params = list()
  model$opt_tv_params = list()
  model = update_tmb_indices(model)
  model
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
#' @param na_is_zero should NAs be replaced with zeros?
#' @importFrom testthat testthat_tolerance
#' @export
compare_sims = function(classic_sim, tmb_sim, tolerance = NULL, compare_attr = TRUE, na_is_zero = FALSE) {
  if (is.null(tolerance)) {
    #if (require(testthat)) {
      tolerance = testthat_tolerance()
    #} else {
    #  tolerance = .Machine$double.eps^0.5
    #}
  }
  if (na_is_zero) {
    classic_sim[is.na(classic_sim)] = 0
    tmb_sim[is.na(tmb_sim)] = 0
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

#' @rdname compare_sims
#' @export
compare_sims_cbind = function(classic_sim, tmb_sim, index) {
  stopifnot(length(index) == 1L)
  x = cbind(unlist(classic_sim[,index]),
        unlist(tmb_sim[,index]))
  colnames(x) = c('classic', 'tmb') %_% names(classic_sim[index])
  x
}

#' @rdname compare_sims
#' @export
compare_sims_plot = function(...) {
  x = compare_sims_cbind(...)
  ss = sub('^classic_', '', colnames(x))
  plot(x[,1], type = 'l', main = ss[[1L]][1])
  lines(x[,2], col = 'red', lty = 2)
}

#' @rdname compare_sims
#' @export
compare_sims_diffmat = function(classic_sim, tmb_sim, state, op = `-`) {
  op(
    as.matrix(classic_sim[,names(state)]),
    as.matrix(tmb_sim[,names(state)]))
}

#' @rdname compare_sims
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
#' @param tmb_pars optional parameters to pass to the tmb functions
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
compare_grads = function(model, tolerance = 1e-5, tmb_pars = NULL, ...) {

  args = list(...)
  if (!'method.args' %in% names(args)) {
    args$method.args = list(
        eps = 1e-4, d = 0.1,
        zero.tol = sqrt(.Machine$double.eps / 7e-7), r = 4, v = 2,
        show.details = FALSE
      )
  }
  dd = tmb_fun(model)
  if (is.null(tmb_pars)) {
    tmb_pars = dd$par
  }
  current <- numDeriv::grad(dd$fn, tmb_pars, ...)
  target <- dd$gr(tmb_pars)
  attributes(target) <- attributes(current) <- NULL
  if (is.na(tolerance)) return(nlist(target, current))
  all.equal(target, current, tolerance)
}

#' @rdname compare_grads
#' @export
compare_hessians = function(model, tolerance = 1e-5, tmb_pars = NULL, ...) {
  args = list(...)
  if (!'method.args' %in% names(args)) {
    args$method.args = list(
      eps = 1e-4, d = 0.1,
      zero.tol = sqrt(.Machine$double.eps / 7e-7), r = 4, v = 2,
      show.details = FALSE
    )
  }
  dd = tmb_fun(model)
  if (is.null(tmb_pars)) {
    tmb_pars = dd$par
  }
  current <- numDeriv::hessian(dd$fn, tmb_pars, ...)
  target <- dd$he(tmb_pars)
  dms = dim(target)
  stopifnot(all.equal(dms, dim(current)))
  attributes(target) <- attributes(current) <- NULL
  dim(target) = dim(current) = dms
  if (is.na(tolerance)) return(nlist(target, current))
  all.equal(target, current, tolerance)
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
    MP_vax_make_state_with_hazard = FALSE,
    MP_force_dgTMatrix = TRUE)
}

##' R Mode
##'
##' Set options so that TMB engine runs in a manner
##' that is comparable with the R engine
##'
##' @export
r_mode = function() {
  options(
    MP_tmb_models_match_r = TRUE,
    MP_default_do_hazard = TRUE,
    MP_default_do_make_state = TRUE
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

##' Factory Fresh MacPan Options
##' @export
factory_fresh_macpan_options = function() {
  flex_version <- readLines(
        system.file("tmb/recommended_spec_version",
            package = "McMasterPandemic"
        )
    )[1]
  options(

        # -- spec version settings ---------------------------------------------

        # default spec version
        MP_flex_spec_version = flex_version,

        # url containing spec version descriptions
        MP_flex_spec_doc_site = "https://canmod.net/misc/flex_specs",

        # default object file for c++ implementation of the spec
        #   - the default means that src/McMasterPandemic.src will be used
        #   - a common alternative is "macpan", which will correspond to
        #     different spec versions in inst/tmb
        MP_flex_spec_dll = "McMasterPandemic",

        # https://stackoverflow.com/questions/8396577/check-if-character-value-is-a-valid-r-object-name
        # need to wrap this in wrap_exact for names that are not embedded in
        # character strings representing expressions
        MP_name_search_regex = "((([[:alpha:]]|[.][._[:alpha:]])[._[:alnum:]]*)|[.])",

        # -- default behaviour of classic macpan engine ------------------------

        # ??
        MP_badsum_action = "warning",

        # ??
        MP_badsum_tol = 1e-12,

        # default number of steps to take in the iterative eigenvector solver
        MP_rexp_steps_default = 100,

        # -- warnings associated with flexmodel structure ----------------------

        # warn if there are multiple rates for the same state transition
        MP_warn_repeated_rates = FALSE,

        # warn if there is no outflow in the model (rarely what users want)
        # -- setting to FALSE because of MP_auto_outflow
        MP_warn_no_outflow = FALSE,

        # warn if there are time-variation breakpoints outside of the
        # simulation bounds
        MP_warn_bad_breaks = TRUE,

        # -- flexmodel defaults ------------------------------------------------
        # -- see ?flexmodel for description of these arguments -----------------
        MP_default_do_hazard = FALSE,
        MP_default_do_hazard_lin = FALSE,
        MP_default_do_approx_hazard = FALSE,
        MP_default_do_approx_hazard_lin = FALSE,
        MP_default_do_make_state = FALSE,
        MP_default_do_sim_constraint = FALSE,
        MP_default_sim_lower_bound = 1e-12,

        # -- tmb_fun behaviour -------------------------------------------------

        # if the user does not specify any outflows, automatically add
        # an outflow for every inflow
        MP_auto_outflow = TRUE,

        # automatically update the tmb indices when (re)generating the
        # tmb ad fun
        MP_auto_tmb_index_update = TRUE,

        # do state condensation with c++ code as opposed to r
        MP_condense_cpp = TRUE,

        MP_silent_tmb_function = TRUE,

        # optimizer options ----------------------------------------------------
        MP_get_bbmle_init_from_nlminb = FALSE,

        # -- control how comparable r and tmb engines are ----------------------
        # -- see r_tmb_comparable ----------------------------------------------
        MP_use_state_rounding = TRUE,         # FALSE ~ comparable
        MP_vax_make_state_with_hazard = TRUE, # FALSE ~ comparable
        MP_tmb_models_match_r = FALSE,        # TRUE  ~ comparable
        MP_force_dgTMatrix = FALSE            # TRUE  ~ comparable
    )
}

# converting between classic and tmb params_timevar -------

#   this is annoying:
#   macpan_ontario orders by Symbol then Date
#   macpan.cpp orders by breaks then spi (which is roughly Date then Symbol)

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
