## mapping
## (1) dependence of elements of the rate matrix on parameters
##   (so we know which elements of the rate matrix to update when
##    a parameter changes)
## (2) how to construct elements of the rate matrix from
##    parameters (list of indices of multiplied parameters,
##    binary/boolean that specifies whether the parameter's complement
##   (1-x) is taken before multiplying
afun_ind <- function(from, to, val) {
  pos <- pfun(from, to, M)
  v <- substitute(val)
  pars <- vapply(all.vars(v),
                 FUN=match, table=params,
                 FUN.VALUE=numeric(1))
  ## we just need to identify complementarity
  ## do the Bad Thing for now (i.e., deparse + string operations
  ## (alternative is recursive formula processing, ugh.)
  s <- strsplit(deparse(v), "\\*")[[1]]
  complement <- grepl("1 +-",s)
  ratemat_vals <<- append(ratemat_vals, list(pos, pars, complement))
  for (p in pars) {
    ## not quite sure how to make this work ...
    ## ratemat_deps[[p]] <<- append(ratemate_deps[[p]], pos)
    ## ratemat_deps[[p]] <- append(ratemate_deps[[p]], list(pos))
    ## error in append_ratemat(...): object 'ratemate_deps' not found
    ## some weird scoping thing going on here with the `[[<-` operator?
    ratemate_deps <- ratemat_deps
    ratemat_deps[[p]] <- c(ratemate_deps[[p]], list(pos))
  }
  assign("ratemat_deps", ratemat_deps, parent.frame())
  ## add these to a global list
  return(NULL)
}

params <- c("x","y","z")
M <- matrix(0,5,5,dimnames=list(LETTERS[1:5], LETTERS[1:5]))
pfun <- McMasterPandemic:::pfun

## initialize data structures for dependence
ratemat_deps <- vector("list", length(params))
ratemat_vals <- list()

afun_ind("A", "B", x*y*(1-z))
afun_ind("B", "C", z*(1-y))

## construct dependence matrix for graphs?


## New Stuff ----

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

opt_list = 'beta0'
params = McMasterPandemic::read_params('ICU1.csv')
params_to_opt = params[opt_list]
params_to_fix = params[!(names(params) %in% opt_list)]
state = McMasterPandemic::make_state(params = params)
M = McMasterPandemic::make_ratemat(state, params)
state_params = c(state, params)
pfun = McMasterPandemic:::pfun

ratemat_struc = lapply(ratemat_struc, product_list)

product_list = function(x) {
  x$factors =
    x$formula %>%
    as.character %>%
    getElement(2L) %>%
    strsplit(split = '\\+') %>%
    getElement(1L) %>%
    strsplit(split = '\\*') %>%
    lapply(factor_table) %>%
    bind_rows(.id = 'product') %>%
    mutate(value = c(state, params)[var]) %>%
    mutate(factor = ifelse(compl, 1 - value, ifelse(invrs, 1 / value, value)))
  x$pos =
    do.call(McMasterPandemic:::pfun, c(x[c('from', 'to')], list(mat = M)))
  x
}

factor_table = function(x) {
  data.frame(
    var = unlist(lapply(x, get_variables)),
    compl = unlist(lapply(x, find_operators, '-')),
    invrs = unlist(lapply(x, find_operators, '/')))
}

factors =
  ratemat_struc %>%
  lapply(`[[`, "factors") %>%
  bind_rows %>%
  select(var, compl, invrs) %>%
  distinct %>%
  mutate(indices = apply(outer(var, names(state_params), '=='), 1, which))

state_params[factors$indices]

variables[variables %in% names(state)]
variables[variables %in% names(params_to_fix)]
variables[variables %in% names(params_to_opt)]

variable_regex = paste0('(', paste0(c(names(params), names(state)), collapse = '|'), ')', sep = '')
get_variables = function(x) {
  r = regexpr(variable_regex, x)
  regmatches(x, r)
}
find_operators = function(x, operator) {
  grepl(paste0('\\( *1 *', operator, ' *', variable_regex, collapse = ''), x)
}

ratemat_struc = list(
  list(from = "E", to = "Ia", formula = ~ (alpha) * (sigma)),
  list(from = "E", to = "Ip", formula = ~ (1 - alpha) * (sigma)),
  list(from = "Ia", to = "R", formula = ~ (gamma_a)),
  list(from = "Ip", to = "Im", formula = ~ (mu) * (gamma_p)),
  list(from = "Ip", to = "Is", formula = ~ (1 - mu) * (gamma_p)),
  list(from = "Im", to = "R", formula = ~ (gamma_m)),
  list(from = "Is", to = "H", formula = ~ (1 - nonhosp_mort) * (phi1) * (gamma_s)),
  list(from = "Is", to = "ICUs", formula = ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s)),
  list(from = "Is", to = "ICUd", formula = ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s)),
  list(from = "Is", to = "D", formula = ~ (nonhosp_mort) * (gamma_s)),
  list(from = "ICUs", to = "H2", formula = ~ (psi1)), ## ICU to post-ICU acute care
  list(from = "ICUd", to = "D", formula = ~ (psi2)), ## ICU to death
  list(from = "H2", to = "R", formula = ~ (psi3)), ## post-ICU to discharge
  ## H now means 'acute care' only; all H survive & are discharged
  #list(from = "H", to = "D", formula = ~ 0),
  list(from = "H", to = "R", formula = ~ (rho)), ## all acute-care survive
  list(from = "Is", to = "X", formula = ~ (1 - nonhosp_mort) * (phi1) * (gamma_s)), ## assuming that hosp admissions mean *all* (acute-care + ICU)
  # force of infection
  list(from = "S",  to = "E",
       formula = ~
         (Ia) * (beta0) * (1/N) * (Ca) +
         (Ip) * (beta0) * (1/N) * (Cp) +
         (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
         (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
) %>%
  lapply(product_list)

mk_indices = function(x) {
  apply(outer(x, names(state_params), '=='), 1, which)
}

factors =
  ratemat_struc %>%
  lapply(`[[`, "factors") %>%
  bind_rows %>%
  select(var, compl, invrs) %>%
  distinct %>%
  mutate(state_param_indices = mk_indices(var)) %>%
  mutate(factor_indices = seq_along(var))

update_w_indices = function(x) {
  x$factors = left_join(x$factors, factors, by = c('var', 'compl', 'invrs'))
  x
}
ratemat_struc = lapply(ratemat_struc, update_w_indices)
