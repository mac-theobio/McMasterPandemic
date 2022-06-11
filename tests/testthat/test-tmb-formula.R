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
  '^sums in opt_param formulas is under construction'
)
expect_error(
  McMasterPandemic:::parse_opt_form(1 ~ 1 ~ 1),
  '^one may not use more than one tilde in optimization formulas'
)
expect_error(
  McMasterPandemic:::parse_opt_form(logit_mu ~ normal(1:2, 2:1) + flat(10)),
  '^sums in opt_param formulas is under construction'
)
expect_error(
  McMasterPandemic:::parse_and_resolve_opt_form(
    logit_mu + beta0 ~ normal(1:2, 2:1),
    params_icu
  ),
  '^sums in opt_param formulas is under construction'
)
expect_error(
  McMasterPandemic:::parse_and_resolve_opt_form(
    notaparameter + mu ~ normal(1:2, 2:1),
    params_icu
  ),
  '^sums in opt_param formulas is under construction'
)

