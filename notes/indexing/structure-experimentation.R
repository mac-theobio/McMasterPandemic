library(McMasterPandemic)
library(dplyr)

grep("^MP", names(options()), value = TRUE)
options(MP_flex_spec_version = "0.1.0")
options(macpan_pfun_method = "grep")

start_date <- "2020-02-01"
end_date <- "2020-09-01"

## initialize params
base_params <- read_params("PHAC.csv")
vax_params <- expand_params_vax(
  params = base_params,
  model_type = "twodose"
)
base_state <- make_state(params = base_params)
vax_state <- expand_state_vax(
  x = base_state,
  model_type = "twodose",
  unif = FALSE
)
vax_sim <- run_sim(
  params = vax_params,
  state = vax_state,
  start_date = start_date,
  end_date = end_date,
  condense_args = list(keep_all = TRUE)
)

# TODO:
#   - add multi_rate, which would be very similar to rep_rate
#     but have vectors for formulae



asymp_unvax_N = sum(struc('S_unvax', 'E_unvax', 'Ia_unvax', 'Ip_unvax', 'R_unvax'))
asymp_vaxprotect1_N = sum(struc('S_vaxprotect1', 'E_vaxprotect1', 'Ia_vaxprotect1', 'Ip_vaxprotect1', 'R_vaxprotect1'))



x <- params[["vax_prop_first_dose"]] * params[["vax_doses_per_day"]]/asymp_unvax_N
x[is.nan(x)] <- 0
vax_rate$dose1 <- x

x <- (1 - params[["vax_prop_first_dose"]]) * params[["vax_doses_per_day"]]/asymp_vaxprotect1_N
x[is.nan(x)] <- 0
vax_rate$dose2 <- x

// a, indices giving the states to update
// b, counts of other states required for each update
// c, indices giving the states required for each update
//
// a = [30, 31]
// b = [5, 5]
// c = [1, 2, 3, 4, 5, 13, 14, 15, 16, 17]
//
sp = c('asymp_unvax_N', 'S_unvax', 'E_unvax', 'Ia_unvax', 'Ip_unvax', 'R_unvax', 'asymp_vaxprotect1_N', 'S_vaxprotect1', 'E_vaxprotect1', 'Ia_vaxprotect1', 'Ip_vaxprotect1', 'R_vaxprotect1')
a = c(1, 7)
b = c(5, 5)
c = c(2, 3, 4, 5, 6, 8, 9, 10, 11, 12)
sp[a]
sp[c]

v = [A, B]  # two state variables
sv = A + B  # a third state variable that gets updated just before ratemat does
p1, p2 # are two parameters
r = p1/sv + p2/sv  # rate matrix update that depends on





rep_rate(to = "E", from = "Ia", formula = ~ (alpha) * (sigma),
         base_state, base_params,
         make_ratemat(base_state, base_params))

rep_rate("Ip", "Is", ~ (1 - mu) * (gamma_p),
         vax_state, vax_params,
         make_ratemat(vax_state, vax_params))

rep_rate(paste0("S_.*", vax_cat[2]), paste0("S_.*", vax_cat[3]),
         ~ (vax_response_rate),
         vax_state, vax_params,
         make_ratemat(vax_state, vax_params))

rep_rate(paste0("S_.*", vax_cat[2]), paste0("S_.*", vax_cat[3]),
         ~ (vax_response_rate),
         vax_state, vax_params,
         make_ratemat(vax_state, vax_params))

rate(paste0("Ip_.*", vax_cat[4]),
         paste0("Is_.*", vax_cat[4]),
         ~ (1 - vax_mu_dose1) * gamma_p,
         vax_state, vax_params,
         make_ratemat(vax_state, vax_params))


model <- (init_model
          %>% do.call(flex_args)
          %>% add_rate("E", "Ia", ~ (alpha) * (sigma))
          %>% add_rate("E", "Ip", ~ (1 - alpha) * (sigma))
          %>% add_rate("Ia", "R", ~ (gamma_a))
          %>% add_rate("Ip", "Im", ~ (mu) * (gamma_p))
          %>% add_rate("Ip", "Is", ~ (1 - mu) * (gamma_p))
          %>% add_rate("Im", "R", ~ (gamma_m))
          %>% add_rate("Is", "H", ~
                         (1 - nonhosp_mort) * (phi1) * (gamma_s))
          %>% add_rate("Is", "ICUs", ~
                         (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * (gamma_s))
          %>% add_rate("Is", "ICUd", ~
                         (1 - nonhosp_mort) * (1 - phi1) * (phi2) * (gamma_s))
          %>% add_rate("Is", "D", ~ (nonhosp_mort) * (gamma_s))
          %>% add_rate("ICUs", "H2", ~ (psi1))
          %>% add_rate("ICUd", "D", ~ (psi2))
          %>% add_rate("H2", "R", ~ (psi3))
          %>% add_rate("H", "R", ~ (rho))
          %>% add_rate("Is", "X", ~ (1 - nonhosp_mort) * (phi1) * (gamma_s))
          %>% add_rate("S", "E", ~
                         (Ia) * (beta0) * (1 / N) * (Ca) +
                         (Ip) * (beta0) * (1 / N) * (Cp) +
                         (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
                         (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
          %>% add_parallel_accumulators(c("X", "N", "P", "V"))
          %>% update_tmb_indices()
)

# force of infection under vaxination
#  -- note the need to use brackets so that the matrix product (ie the sums)
#     is the final step (this is because we don't see the + signs within
#     elements when we multiply)
A  =   struc("(1)", "(1 - VE1)", "(1 - VE2)")
B  = t(struc("(Ca)", "(Cp)", "(Cm)", "(Cs)"))
C  =   struc("(Ia)", "(Ip)", "(Im)", "(Is)")
AB = kronecker(A, B)
foi = AB %*% C * struc('(beta) * (1/N)')



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
  %>% update_tmb_indices()
)
test_model$tmb_indices$sum_indices
dd <- tmb_fun(test_model)
length(dd$par)
length(params)
dd$fn(dd$par)
sum(state)

v = "(alpha) * (sigma) * (1 / Isum)"

nms = names(c(state, params))
r = regexpr('[A-z]+', v)
regmatches(v, r)
