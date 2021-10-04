library(McMasterPandemic)
library(TMB)
library(tools)
library(dplyr)
library(semver)

spec_version = "0.1.0"
options(MP_flex_spec_version = spec_version)

test_files = "../../inst/tmb/"

cpp = file.path(test_files, spec_version, "macpan.cpp")
dll = file_path_sans_ext(cpp)

compile(cpp)
dyn.load(dynlib(dll))

#grep("^MP", names(options()), value = TRUE)
#options(MP_flex_spec_version = "0.1.0")
#options(macpan_pfun_method = "grep")

params <- read_params("ICU1.csv")
state <- make_state(params = params)
M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

test_model <- (init_model(
  params, state,
  start_date = "2021-09-09", end_date = "2021-10-09",
  params_timevar = data.frame(
    Date = c("2021-09-15", "2021-09-15", "2021-10-05"),
    Symbol = c("beta0", "mu", "beta0"),
    Value = c(0.5, 0.1, 0.05),
    Type = c('rel_orig', 'rel_prev', 'rel_orig')
  )
)
%>% add_state_param_sum("Isum", "^I[a-z]")
%>% add_state_param_sum("ICUsum", c("ICUs", "ICUd"))
%>% add_rate("E", "Ia", ~ (alpha) * (sigma) * (1 / Isum))
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
%>% add_rate("H", "R", ~ (rho) * (1 / ICUsum))
%>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
%>% add_rate("S", "E", ~
               (Ia) * (beta0) * (1 / N) * (Ca) +
               (Ip) * (beta0) * (1 / N) * (Cp) +
               (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
               (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
%>% add_parallel_accumulators(c("X", "N", "P", "V"))
%>% add_tmb_indices()
)

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
tv_mult = test_model$timevar$piece_wise$schedule$Value
tv_orig = test_model$timevar$piece_wise$schedule$Type == "rel_orig"

do_hazard = TRUE

dd <- MakeADFun(data = list(state = c(state),
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
                            tv_mult = tv_mult,
                            tv_orig = tv_orig,
                            sumidx = test_model$tmb_indices$sum_indices$sumidx,
                            sumcount = test_model$tmb_indices$sum_indices$sumcount,
                            summandidx = test_model$tmb_indices$sum_indices$summandidx,
                            par_accum_indices = par_accum_indices,
                            do_hazard = do_hazard,
                            numIterations = test_model$iters),
                parameters = list(params=c(test_model$params)),
                DLL=basename(dll))

tmb_sparse_ratemat = dd$report()$ratemat
concatenated_state_vector = dd$report()$concatenated_state_vector

print(concatenated_state_vector)

