library(McMasterPandemic)
library(dplyr)

grep("^MP", names(options()), value = TRUE)
options(MP_flex_spec_version = "0.1.0")
options(macpan_pfun_method = "grep")

params <- read_params("ICU1.csv")
state <- make_state(params = params)
vax_params <- expand_params_vax(
  params = params,
  model_type = "twodose"
)
vax_state <- expand_state_vax(
  x = state,
  model_type = "twodose",
  unif = FALSE
)

# problem dimensions
(epi_states = c(attr(vax_state, "epi_cat")))
(asymp_cat = c("S", "E", "Ia", "Ip", "R"))
(vax_cat = c(attr(vax_state, "vax_cat")))
(dose_from = rep(asymp_cat, 2))
(dose_to = c(asymp_cat, rep("V", 5)))

# Specify structure of the force of infection calculation
Istate = (McMasterPandemic:::expand_names(
  c('Ia', 'Ip', 'Im', 'Is'),   # = Icats
  attr(vax_params, 'vax_cat')) # = vax_cats
  %>% struc
)
baseline_trans_rates =
  struc(
                   'Ca',
                   'Cp',
    '(1 - iso_m) * (Cm)',
    '(1 - iso_s) * (Cs)') *
  struc('(beta0) * (1/N)')
vax_trans_red = struc_block(struc(
  '1',
  '1',
  '(1 - vax_efficacy_dose1)',
  '(1 - vax_efficacy_dose1)',
  '(1 - vax_efficacy_dose2)'),
  row_times = 1, col_times = 5)

alpha   = c("alpha", "alpha", "vax_alpha_dose1", "vax_alpha_dose1", "vax_alpha_dose2")
mu      = c("mu",    "mu",    "vax_mu_dose1",    "vax_mu_dose1",    "vax_mu_dose2")
sigma   = struc("sigma")
gamma_p = struc("gamma_p")
E_to_Ia_rates  = struc(           alpha ) * sigma
E_to_Ip_rates  = struc(complement(alpha)) * sigma
Ip_to_Im_rates = struc(              mu ) * gamma_p
Ip_to_Is_rates = struc(complement(   mu)) * gamma_p

vax_params_test = vax_params
vax_state_test = vax_state
#vax_params_test$beta0 = 0
#vax_params_test$vax_prop_first_dose = 1
#vax_state_test[-1] = 0
#vax_state_test[1] = vax_params_test$N
#vax_state_test[] = vax_state[] + runif(length(vax_state))

test_model <- (init_model(
  vax_params, vax_state,
  start_date = "2021-09-09", end_date = "2021-10-09"
)

  # Flow within vaccination categories,
  # with constant rates across categories
  %>% rep_rate("Ia",   "R",    ~                      (gamma_a))
  %>% rep_rate("Im",   "R",    ~                      (gamma_m))
  %>% rep_rate("Is",   "D",    ~ (    nonhosp_mort) * (gamma_s))
  %>% rep_rate("Is",   "H",    ~ (1 - nonhosp_mort) * (gamma_s) * (    phi1))
  %>% rep_rate("Is",   "X",    ~ (1 - nonhosp_mort) * (gamma_s) * (    phi1))
  %>% rep_rate("Is",   "ICUs", ~ (1 - nonhosp_mort) * (gamma_s) * (1 - phi1) * (1 - phi2))
  %>% rep_rate("Is",   "ICUd", ~ (1 - nonhosp_mort) * (gamma_s) * (1 - phi1) * (    phi2))
  %>% rep_rate("ICUs", "H2",   ~                                  (    psi1))
  %>% rep_rate("ICUd", "D",    ~                                  (    psi2))
  %>% rep_rate("H2",   "R",    ~                                  (    psi3))
  %>% rep_rate("H",    "R",    ~ (rho))

  # Flow within vaccination categories,
  # with rates that depend on category
  # (see struc objects created above)
  %>% vec_rate("E", "Ia",  E_to_Ia_rates)
  %>% vec_rate("E", "Ip",  E_to_Ip_rates)
  %>% vec_rate("Ip", "Im", Ip_to_Im_rates)
  %>% vec_rate("Ip", "Is", Ip_to_Is_rates)

  # Vaccination Response Rates
  %>% add_rate("R_vaxdose1", "R_vaxprotect1",  ~ (vax_response_rate_R))
  %>% add_rate("R_vaxdose2", "R_vaxprotect2",  ~ (vax_response_rate_R))
  %>% add_rate("S_vaxdose1", "S_vaxprotect1",  ~ (vax_response_rate))
  %>% add_rate("S_vaxdose2", "S_vaxprotect2",  ~ (vax_response_rate))

  # Forces of Infection
  %>% vec_rate(
    "S" %_% vax_cat,
    "E" %_% vax_cat,
    kronecker(vax_trans_red, t(baseline_trans_rates)) %*% Istate
  )

  # Sums across vaccination categories
  %>% add_state_param_sum("asymp_unvax_N",       asymp_cat %_% "unvax")
  %>% add_state_param_sum("asymp_vaxprotect1_N", asymp_cat %_% "vaxprotect1")

  # Flow among vaccination categories
  # (see dose_* above for epi states that are involved)
  %>% rep_rate(
    dose_from %_% 'unvax',
    dose_to   %_% 'vaxdose1',
    ~ (    vax_prop_first_dose) * (vax_doses_per_day) * (1 / asymp_unvax_N))
  %>% rep_rate(
    dose_from %_% 'vaxprotect1',
    dose_to   %_% 'vaxdose2',
    ~ (1 - vax_prop_first_dose) * (vax_doses_per_day) * (1 / asymp_vaxprotect1_N))

  # Technical steps
  %>% add_parallel_accumulators(c('V' %_% vax_cat, 'X' %_% vax_cat))
  #%>% add_parallel_accumulators(c("X", "N", "P", "V"))
  %>% update_tmb_indices()
)

updt =
c(dose_from %_% 'unvax',    "S" %_% vax_cat, dose_from %_% 'vaxprotect1') %_% 'to' %_%
c(dose_to   %_% 'vaxdose1', "E" %_% vax_cat, dose_to   %_% 'vaxdose2')

all(updt %in% names(test_model$tmb_indices$updateidx))
all(names(test_model$tmb_indices$updateidx) %in% updt)

names(test_model$tmb_indices$make_ratemat_indices)
for(i in 1:104) {
  test_model$rates[i]
  if(any(get_indices_per_rate(test_model, i)$modifier == 5)) break
}

indx_comp = sapply(1:104, compute_rate_from_indices, model = test_model)
expr_comp = compute_rates(test_model)
bads = which(abs(indx_comp - expr_comp) > 1e-5)
test_model$rates[bads]
indx_comp[bads]
expr_comp[bads]
plot(indx_comp, expr_comp)
abline(a = 0, b = 1)
c(unlist(vax_state_test), unlist(vax_params_test))[93]

get_indices_per_rate(test_model, 81)

$from
[1] 15

$to
[1] 16

$count
[1] 90

$spi
[1] 72 71 85  3 73 71 85  4 88 74 71 85  5 89 75 71 85  6 72 71 85 17 73 71 85 18 88 74 71 85 19 89 75 71 85 20 72
[38] 71 85 31 73 71 85 32 88 74 71 85 33 89 75 71 85 34 72 71 85 45 73 71 85 46 88 74 71 85 47 89 75 71 85 48 72 71
[75] 85 59 73 71 85 60 88 74 71 85 61 89 75 71 85 62

$modifier
[1] 0 0 2 0 4 0 2 0 5 0 0 2 0 5 0 0 2 0 4 0 2 0 4 0 2 0 5 0 0 2 0 5 0 0 2 0 4 0 2 0 4 0 2 0 5 0 0 2 0 5 0 0 2 0 4 0
[57] 2 0 4 0 2 0 5 0 0 2 0 5 0 0 2 0 4 0 2 0 4 0 2 0 5 0 0 2 0 5 0 0 2 0

compute_rate_from_indices(test_model, bads[1])

check_rates(test_model)

tmb_strt = Sys.time()
tmb_sim <- run_sim(
  params = vax_params_test, state = vax_state,
  start_date = "2021-09-09", end_date = "2021-10-09",
  condense = FALSE,
  step_args = list(do_hazard = TRUE),
  flexmodel = test_model,
  use_flex = TRUE
)
tmb_nd = Sys.time()
r_strt = Sys.time()
r_sim = run_sim(
  params = vax_params_test, state = vax_state,
  start_date = "2021-09-09", end_date = "2021-10-09",
  condense = FALSE,
  step_args = list(do_hazard = TRUE)
)
r_nd = Sys.time()
as.numeric(r_nd - r_strt) / as.numeric(tmb_nd - tmb_strt)

plot(tmb_sim$S_unvax, type = 'l')
lines(r_sim$S_unvax, col = 'red')

#vax_trans_red
#unvax        1.0  1.0  1.0  1.0  1.0
#vaxdose1     1.0  1.0  1.0  1.0  1.0
#vaxprotect1  0.4  0.4  0.4  0.4  0.4
#vaxdose2     0.4  0.4  0.4  0.4  0.4
#vaxprotect2  0.1  0.1  0.1  0.1  0.1
as.matrix(vax_trans_red)
1 - vax_params$vax_efficacy_dose1
1 - vax_params$vax_efficacy_dose2

#beta_0 -- column vector
#Ia           Ip           Im           Is
#6.666667e-07 1.000000e-06 1.000000e-06 1.000000e-06
beta_0 = as.character(baseline_trans_rates)
with(vax_params, eval(parse(text = beta_0[1])))
with(vax_params, eval(parse(text = beta_0[2])))
with(vax_params, eval(parse(text = beta_0[3])))
with(vax_params, eval(parse(text = beta_0[4])))


#beta_0_fastmatrix
#[,1]  [,2]  [,3]  [,4]         [,5]  [,6]  [,7]  [,8]         [,9] [,10] [,11]
#[1,] 6.666667e-07 1e-06 1e-06 1e-06 6.666667e-07 1e-06 1e-06 1e-06 6.666667e-07 1e-06 1e-06
#[2,] 6.666667e-07 1e-06 1e-06 1e-06 6.666667e-07 1e-06 1e-06 1e-06 6.666667e-07 1e-06 1e-06
#[3,] 2.666667e-07 4e-07 4e-07 4e-07 2.666667e-07 4e-07 4e-07 4e-07 2.666667e-07 4e-07 4e-07
#[4,] 2.666667e-07 4e-07 4e-07 4e-07 2.666667e-07 4e-07 4e-07 4e-07 2.666667e-07 4e-07 4e-07
#[5,] 6.666667e-08 1e-07 1e-07 1e-07 6.666667e-08 1e-07 1e-07 1e-07 6.666667e-08 1e-07 1e-07
#[,12]        [,13] [,14] [,15] [,16]        [,17] [,18] [,19] [,20]
#[1,] 1e-06 6.666667e-07 1e-06 1e-06 1e-06 6.666667e-07 1e-06 1e-06 1e-06
#[2,] 1e-06 6.666667e-07 1e-06 1e-06 1e-06 6.666667e-07 1e-06 1e-06 1e-06
#[3,] 4e-07 2.666667e-07 4e-07 4e-07 4e-07 2.666667e-07 4e-07 4e-07 4e-07
#[4,] 4e-07 2.666667e-07 4e-07 4e-07 4e-07 2.666667e-07 4e-07 4e-07 4e-07
#[5,] 1e-07 6.666667e-08 1e-07 1e-07 1e-07 6.666667e-08 1e-07 1e-07 1e-07

kr = as.matrix(kronecker(vax_trans_red, t(baseline_trans_rates)))
with(vax_params, eval(parse(text = kr[3, 9])))

#Ia_unvax       Ip_unvax       Im_unvax       Is_unvax    Ia_vaxdose1    Ip_vaxdose1
#1.933333e-06   2.900000e-06   2.900000e-06   2.900000e-06   1.933333e-06   2.900000e-06
#Im_vaxdose1    Is_vaxdose1 Ia_vaxprotect1 Ip_vaxprotect1 Im_vaxprotect1 Is_vaxprotect1
#2.900000e-06   2.900000e-06   1.933333e-06   2.900000e-06   2.900000e-06   2.900000e-06
#Ia_vaxdose2    Ip_vaxdose2    Im_vaxdose2    Is_vaxdose2 Ia_vaxprotect2 Ip_vaxprotect2
#1.933333e-06   2.900000e-06   2.900000e-06   2.900000e-06   1.933333e-06   2.900000e-06
#Im_vaxprotect2 Is_vaxprotect2
#2.900000e-06   2.900000e-06

II = as.character(Istate)
with(as.list(vax_state), eval(parse(text = II[1])))
c(struc_eval(Istate, as.list(vax_state))) == vax_state[gsub("\\(|\\)", "", II)]


tv_dat <- data.frame(
  Date = c("2021-09-15", "2021-09-20", "2021-10-05"),
  Symbol = c("beta0", "beta0", "beta0"),
  Value = c(0.5, 0.1, 0.05),
  Type = c("rel_prev", "rel_orig", "rel_prev")
)

test_model <- (init_model(
  params, state,
  start_date = "2021-09-09", end_date = "2021-10-09",
  params_timevar = tv_dat)
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
  #%>% add_rate("S", "E", ~
  #               (Ia) * (beta0) * (1 / N) * (Ca) +
  #               (Ip) * (beta0) * (1 / N) * (Cp) +
  #               (Im) * (beta0) * (1 / N) * (Cm) * (1 - iso_m) +
  #               (Is) * (beta0) * (1 / N) * (Cs) * (1 - iso_s))
  %>% add_parallel_accumulators(c("X", "N", "P", "V"))
  %>% update_tmb_indices()
)

# 4 * add + 2 * invrs + compl

#   spi modifier
#1    3        0    000 unmodified Ia
#2   15        0    000 unmodified beta0
#3   29        2    010 inverse N
#4   16        0    000 unmodified Ca
#5    4        4    100 unmodified Ia at the start of a product
#6   15        0    000 unmodified beta0
#7   29        2    010 inverse N
#8   17        0    000 unmodified Cp
#9    5        4    100 unmodified Im at the start of a product
#10  15        0    000 unmodified beta0
#11  29        2    010 inverse N
#12  18        0    000 unmodified Cm
#13  32        1    001 complement iso_m
#14   6        4    010 unmodified Is at the start of a product
#15  15        0    000 unmodified beta0
#16  29        2    010 inverse N
#17  19        0    000 unmodified Cs
#18  33        1    001 complement iso_s


data.frame(spi = test_model$tmb_indices$make_ratemat_indices$spi[48 - (17:0)],
           modifier = test_model$tmb_indices$make_ratemat_indices$modifier[48 - (17:0)])
