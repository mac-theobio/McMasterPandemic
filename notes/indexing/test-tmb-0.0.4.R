library(testthat)
library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)

spec_version = "0.0.4"
options(MP_flex_spec_version = spec_version)

test_files = "../../inst/tmb/"

cpp = file.path(test_files, spec_version, "macpan.cpp")
dll = file_path_sans_ext(cpp)

compile(cpp)
dyn.load(dynlib(dll))

params = read_params('ICU1.csv')
state = make_state(params = params)
M = McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

test_model = (init_model(
  params, state,
  start_date = "2021-09-09", end_date = "2021-10-09",
  timevar_piece_wise = data.frame(
    Date = c("2021-09-15", "2021-09-15", "2021-10-05"),
    Symbol = c("beta0", "mu", "beta0"),
    Value = c(0.5, 0.1, 0.05),
    Type = c('rel_orig', 'rel_prev', 'rel_orig')
  ))
  %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
  %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
  %>% add_rate("Ia", "R", ~ (gamma_a))
  %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
  %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
  %>% add_rate("Im", "R", ~ (gamma_m))
  %>% add_rate("Is", "H", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
  %>% add_rate("Is", "ICUs", ~ (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
  %>% add_rate("Is", "ICUd", ~ (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
  %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
  %>% add_rate("ICUs", "H2", ~ (psi1))
  %>% add_rate("ICUd", "D", ~ (psi2))
  %>% add_rate("H2", "R", ~ (psi3))
  %>% add_rate("H", "R", ~ (rho))
  %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
  %>% add_rate("S",  "E", ~
                 (Ia) * (beta0) * (1/N) * (Ca) +
                 (Ip) * (beta0) * (1/N) * (Cp) +
                 (Im) * (beta0) * (1/N) * (Cm) * (1-iso_m) +
                 (Is) * (beta0) * (1/N) * (Cs) * (1-iso_m))
  %>% add_parallel_accumulators(c("X", "N", "P", "V"))
  %>% add_tmb_indices()
)

#sp = c(test_model$state, test_model$params)
from = test_model$tmb_indices$make_ratemat_indices$from
to = test_model$tmb_indices$make_ratemat_indices$to
count = test_model$tmb_indices$make_ratemat_indices$count
spi = test_model$tmb_indices$make_ratemat_indices$spi
modifier = test_model$tmb_indices$make_ratemat_indices$modifier

updateidx = test_model$tmb_indices$updateidx
breaks = test_model$timevar$piece_wise$breaks
count_of_tv_at_breaks = test_model$timevar$piece_wise$count_of_tv_at_breaks
tv_spi = test_model$timevar$piece_wise$schedule$tv_spi
tv_val = test_model$timevar$piece_wise$schedule$tv_val
par_accum_indices = test_model$tmb_indices$par_accum_indices

dd <- MakeADFun(data = list(state = c(state), #c(test_model$state),
                            ratemat = M,
                            from = from,
                            to = to,
                            count = count,
                            spi = spi,
                            modifier = modifier,
                            updateidx = c(updateidx),
                            breaks = breaks,
                            count_of_tv_at_breaks = count_of_tv_at_breaks,
                            tv_spi = tv_spi,
                            tv_val = tv_val,
                            par_accum_indices = par_accum_indices,
                            numIterations = test_model$iters),
                parameters = list(params=c(test_model$params)),
                DLL=basename(dll))

tmb_sparse_ratemat = dd$report()$ratemat
concatenated_state_vector = dd$report()$concatenated_state_vector

print(concatenated_state_vector)
