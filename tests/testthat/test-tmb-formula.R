library(McMasterPandemic)
library(testthat)

params_icu = read_params("ICU1.csv")
params_sir = c(gamma = 0.06, beta = 0.15)

test_formulas = list(
  "logit_^C(a|p)$" ~ log_flat(1), # params
  log_beta0 + log_mu ~ log_normal(matrix(1:4, 2, 2), 1), # tv_mult
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


McMasterPandemic:::parse_opt_form(c('testingtesting', "log_mu") ~ log_normal(1, 2))


expect_error(
  McMasterPandemic:::parse_opt_form("log_beta0" + "mu" ~ normal(0.1, 2)),
  '^parameters must be expressed as valid R names'
)
expect_error(
  McMasterPandemic:::parse_opt_form(1 ~ 1 ~ 1),
  '^one may not use more than one tilde in optimization formulas'
)
expect_error(
  McMasterPandemic:::parse_opt_form(logit_mu ~ normal(1:2, 2:1) + flat(10)),
  '^you cannot apply functions or operators to terms in a sum'
)
expect_error(
  McMasterPandemic:::parse_and_resolve_opt_form(
    logit_mu + beta0 ~ normal(1:2, 2:1),
    params_icu
  ),
  '^\nthe transformation scale used when a parameter is an argument'
)
expect_error(
  McMasterPandemic:::parse_and_resolve_opt_form(
    notaparameter + mu ~ normal(1:2, 2:1),
    params_icu
  ),
  '^cannot find opimization parameters in the params vector'
)

x = 3
f = log_mu ~ log_normal(1, x)
McMasterPandemic:::parse_and_resolve_opt_form(
  f, params_icu
)
f = beta ~ normal(1, x)
mm = make_sir_model(
  params = c(beta = 1, gamma = 1), state = c(S = 1, I = 1, R = 1),
  start_date = '2000-01-01', end_date = '2000-01-01'
) %>% add_opt_params(f)
mm$opt_params


eps = 1e-10
x = exp(seq(from = 1e-12, to = 4e-10, length = 200)) - 1
y = x + eps * exp(-x / eps)
plot(x, y, type = 'l', ylim = c(0, max(y)))
abline(h = eps)
abline(h = 0)
