trans = 'testingtesting' # test to make sure parse_opt_form is evaluating in the right environment
distr = 1
test_formulas = list(
  "logit_^C(a|p)$" ~ log_flat(1), # params
  log_beta0 + log_mu ~ log_normal(matrix(1:4 2, 2), 1), # tv_mult
  logit_mu ~ normal(1:2, 2:1), # tv_mult
  logit_mu + Ca ~ logit_normal(1:2, 2:1), # params
  log_beta0 + logit_mu + Ca + Cp + N ~ normal(1:5, 6:10), # params
  c("log_beta0", "log_^C(a|p)$") ~ log_normal(1, 2), # fail -- can't find parameter
  c(trans, "log_mu") ~ log_normal(1, 2), # fail -- can't find parameter
  c("log_beta0") + mu ~ normal(0.1, 2), # fail -- invalid sum
  "log_beta0" + "mu" ~ normal(0.1, 2), # fail -- invalid sum
  logit_mu ~ normal(1:2, 2:1) + flat(10), # fail -- multiple priors
  logit_mu ~ flat(distr / 10) # param or tv_mult
)

# conflicted::conflict_prefer()


separate_param_trans = function(x, params) {
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

  if (length(trans_nm) == 1L) {
    trans_nm = rep(trans_nm, length(param_nm))
  }

  nlist(param_nm, trans_nm)
}
get_names_and_dims = function(x, params, params_timevar = NULL) {
  nms = separate_param_trans(x)
  unpack(nms)
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
  l = unname(l)
  nlist(param_nm, trans_nm, p, l)
}
separate_prior_trans = function(x) {
  rhs_func = x
  prior_trans = index_sep(rhs_func, 1)
  if (prior_trans %in% valid_trans) {
    prior_family = index_sep(rhs_func, -1)
  } else {
    prior_trans = ''
    prior_family = rhs_func
  }
  if (!prior_family %in% valid_prior_families) {
    prior_family = ''
    prior_trans = ''
  }
  nlist(prior_family, prior_trans)
}
get_hyperparam_list = function(x, n_params, n_breaks) {
  prior_spec = x[[3]]
  if (is.call(prior_spec)) {
    h = length(prior_spec) - 1L
    prior_list = separate_prior_trans(as.character(prior_spec[[1]]))
    prior_spec = prior_spec[-1]
  } else if (is.numeric(prior_spec)) {
    h = 1
    prior_spec = list(prior_spec)
    prior_list = list(prior_family = '', prior_trans = '')
  }
  p = n_params
  l = n_breaks
  ll = vector("list", max(1, l) * p)
  dim(ll) = c(max(1, l), p)
  for (k in seq_len(h)) {
    hyper_mat = matrix(eval(prior_spec[[k]]), max(1, l), p)
    for (i in seq_len(max(1, l))) {
      for (j in seq_len(p)) {
        ll[[i, j]] = c(ll[[i, j]], hyper_mat[i, j])
      }
    }
  }
  c(prior_list, nlist(ll))
}
