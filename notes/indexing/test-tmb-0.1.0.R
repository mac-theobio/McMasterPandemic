library(McMasterPandemic)
library(dplyr)

grep("^MP", names(options()), value = TRUE)
options(MP_flex_spec_version = "0.1.0")
options(macpan_pfun_method = "grep")

params <- read_params("ICU1.csv")
state <- make_state(params = params)
M <- McMasterPandemic::make_ratemat(state, params, sparse = TRUE)

test_model <- (init_model(
  params, state,
  start_date = "2021-09-09", end_date = "2021-10-09"
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
%>% add_rate("H", "R", ~ (rho) / (1 / ICUsum))
%>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
%>% add_rate("S", "E", ~
               (Ia) * (beta0) * (1 / N) * (Ca) +
               (Ip) * (beta0) * (1 / N) * (Cp) +
               (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
               (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
%>% add_parallel_accumulators(c("X", "N", "P", "V"))
%>% add_tmb_indices()
)


test_model$tmb_indices$sum_indices
