library(ggplot2)
library(McMasterPandemic)
library(dplyr)

options(macpan_pfun_method = "grep")

do_variant = McMasterPandemic:::do_variant

reset_spec_version()
tmb_mode()
options(macpan_pfun_method = "grep")

start_date <- "2021-02-01"
end_date <- "2021-09-01"

base_params <- read_params("PHAC.csv")
vax_params <- expand_params_vax(
  params = base_params,
  vax_doses_per_day = 1e4,
  model_type = "twodose"
)

params <- expand_params_variant(
  vax_params,
  variant_prop = 0.5,
  variant_advantage = 1.5,
  variant_vax_efficacy_dose1 = 0.3,
  variant_vax_efficacy_dose2 = 0.8
)

params$vax_prop_first_dose = 0.7
params$vax_prop_second_dose = 0.2
params$vax_alpha_dose3 = 0.5
params$vax_mu_dose3 = 1
params$vax_efficacy_dose3 = 0.95
params$variant_vax_efficacy_dose3 = 0.9
params$chal_waning_unvax = 0.01
params$chal_waning_dose1 = 0.01
params$chal_waning_dose2 = 0.01
params$chal_waning_dose3 = 0.01

params$leakiness = 1 ## back to original model (in theory)

do_variant = McMasterPandemic:::do_variant(params)

epi_states = c(
  "S", "E", "Ia", "Ip", "Im", "Is",
  "H", "H2", "ICUs", "ICUd", "D", "R",
  "X", "V") # 14 base epidemiological categories
(asymp_cat = c("S", "E", "Ia", "Ip", "R")) # 5 asymptomatic categories
vax_cat = c(
  "unvax",
  "vaxdose1", "vaxprotect1",
  "vaxdose2", "vaxprotect2",
  "vaxdose3", "vaxprotect3"
)
attr(params, "vax_cat") = vax_cat
(accum = c("X", "V")) # two base parallel accumulator states
(non_accum = base::setdiff(epi_states, accum)) # 12 base non-parallel accumulator states
(non_accum_non_S = non_accum[-1]) # 11 base non-susceptible/non-accumulator states

state = layered_zero_state(epi_states, vax_cat)

# dosing transitions across vaccination layers
(dose_from = rep(asymp_cat, 2))
(dose_to = c(asymp_cat, rep("V", 5)))

# ---------------------------
# Specify structure of the force of infection calculation
# ---------------------------

Istate = (c('Ia', 'Ip', 'Im', 'Is')
          %>% expand_names(vax_cat)
          %>% vec
)
Istate_vax = (c('Ia', 'Ip', 'Im', 'Is')
  %>% expand_names(vax_cat[3:7])
  %>% vec
)
baseline_trans_rates =
  vec(
    'Ca',
    'Cp',
    '(1 - iso_m) * (Cm)',
    '(1 - iso_s) * (Cs)') *
  struc('(beta0) * (1/N)')

if(!do_variant) {
   ## transmission reduction due to vaccine
  vax_trans_red = struc_block(vec(
    '1',
    '1',
    '(1 - vax_efficacy_dose1)',
    '(1 - vax_efficacy_dose1)',
    '(1 - vax_efficacy_dose2)',
    '(1 - vax_efficacy_dose2)',
    '(1 - vax_efficacy_dose3)'),
    row_times = 1, col_times = 7)
  ## proportion of pop protected by the vaccine (when there is some non-leakiness)
  vax_protection = struc_block(vec(
     '0',
     '0',
     '(1 - leakiness) * (vax_efficacy_dose1)',
     '(1 - leakiness) * (vax_efficacy_dose1)',
     '(1 - leakiness) * (vax_efficacy_dose2)',
     '(1 - leakiness) * (vax_efficacy_dose2)',
     '(1 - leakiness) * (vax_efficacy_dose3)'),
     row_times = 1, col_times = 7)
} else {
   ## transmission reduction due to vaccine
  vax_trans_red = struc_block(vec(
    '(1 - variant_prop) + (variant_advantage) * (variant_prop)',
    '(1 - variant_prop) + (variant_advantage) * (variant_prop)',
    '(1 - vax_efficacy_dose1) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose1) * (variant_advantage) * (variant_prop)',
    '(1 - vax_efficacy_dose1) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose1) * (variant_advantage) * (variant_prop)',
    '(1 - vax_efficacy_dose2) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose2) * (variant_advantage) * (variant_prop)',
    '(1 - vax_efficacy_dose2) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose2) * (variant_advantage) * (variant_prop)',
    '(1 - vax_efficacy_dose3) * (1 - variant_prop) + (1 - variant_vax_efficacy_dose3) * (variant_advantage) * (variant_prop)'),
    row_times = 1, col_times = 7)
  ## proportion of pop protected by the vaccine (when there is some non-leakiness)
  ## FIXME: ASK STEVE HOW TO WRITE THIS WITH SCALAR MULTIPLICATION OF 1-LEAKINESS
  vax_protection_vec = struc('1 - leakiness') * vec(
     #'0',
     #'0',
     '(vax_efficacy_dose1) * (1 - variant_prop) + (variant_vax_efficacy_dose1) * (variant_advantage) * (variant_prop)',
     '(vax_efficacy_dose1) * (1 - variant_prop) + (variant_vax_efficacy_dose1) * (variant_advantage) * (variant_prop)',
     '(vax_efficacy_dose2) * (1 - variant_prop) + (variant_vax_efficacy_dose2) * (variant_advantage) * (variant_prop)',
     '(vax_efficacy_dose2) * (1 - variant_prop) + (variant_vax_efficacy_dose2) * (variant_advantage) * (variant_prop)',
     '(vax_efficacy_dose3) * (1 - variant_prop) + (variant_vax_efficacy_dose3) * (variant_advantage) * (variant_prop)'
  )
  vax_protection = struc_block(vax_protection_vec, row_times = 1, col_times = 5)
}

chal_waning = vec(
  'chal_waning' %_% c(
    "unvax", "unvax", "dose1", "dose1", "dose2", "dose2", "dose3"
  )
)

alpha   = c("alpha", "alpha", "vax_alpha_dose1", "vax_alpha_dose1", "vax_alpha_dose2", "vax_alpha_dose2", "vax_alpha_dose3")
mu      = c("mu",    "mu",    "vax_mu_dose1",    "vax_mu_dose1",    "vax_mu_dose2",    "vax_mu_dose2",    "vax_mu_dose3")
sigma   = struc("sigma")
gamma_p = struc("gamma_p")
E_to_Ia_rates  = vec(           alpha ) * sigma
E_to_Ip_rates  = vec(complement(alpha)) * sigma
Ip_to_Im_rates = vec(              mu ) * gamma_p
Ip_to_Is_rates = vec(complement(   mu)) * gamma_p

model = (init_model(
      params = expand_params_S0(unlist(params), 1-1e-5),
      state = state,
      start_date = start_date,
      end_date = end_date,
      do_hazard = TRUE,
      do_make_state = TRUE
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
   %>% add_rate("R_vaxdose3", "R_vaxprotect3",  ~ (vax_response_rate_R)) # new
   %>% add_rate("S_vaxdose1", "S_vaxprotect1",  ~ (vax_response_rate))
   %>% add_rate("S_vaxdose2", "S_vaxprotect2",  ~ (vax_response_rate))
   %>% add_rate("S_vaxdose3", "S_vaxprotect3",  ~ (vax_response_rate)) # new

   # Forces of Infection
   %>% vec_rate(
     "S" %_% vax_cat,
     "E" %_% vax_cat,
     kronecker(vax_trans_red, t(baseline_trans_rates)) %*% Istate
   )

   # Non-leakiness in vaccine protection (when leakiness =/= 1)
   # %>% vec_rate(
   #    "S" %_% vax_cat,
   #    "S" %_% vax_cat,
   #    kronecker(leakiness*vax_protection, t(baseline_trans_rates)) %*% Istate
   # )

   %>% vec_rate(
      "S" %_% vax_cat[3:7],
      "R" %_% vax_cat[3:7],
      kronecker(vax_protection, t(baseline_trans_rates)) %*% Istate_vax
   )

   # Sums across vaccination categories
   %>% add_state_param_sum("asymp_unvax_N",       asymp_cat %_% "unvax")
   %>% add_state_param_sum("asymp_vaxprotect1_N", asymp_cat %_% "vaxprotect1")
   %>% add_state_param_sum("asymp_vaxprotect2_N", asymp_cat %_% "vaxprotect2") # new
   %>% add_state_param_sum(
     "vax_prop_first_and_second_dose",
     c("vax_prop_first_dose", "vax_prop_second_dose")
   ) # new

   # Flow among vaccination categories
   # (see dose_* above for epi states that are involved)
   %>% rep_rate(
     dose_from %_% 'unvax',
     dose_to   %_% 'vaxdose1',
     ~ (vax_prop_first_dose) * (vax_doses_per_day) * (1 / asymp_unvax_N))
   %>% rep_rate(
     dose_from %_% 'vaxprotect1',
     dose_to   %_% 'vaxdose2',
     ~ (vax_prop_second_dose) * (vax_doses_per_day) * (1 / asymp_vaxprotect1_N)) # reparameterized
   %>% rep_rate(
     dose_from %_% 'vaxprotect2',
     dose_to   %_% 'vaxdose3',
     ~ (1 - vax_prop_first_and_second_dose) * (vax_doses_per_day) * (1 / asymp_vaxprotect2_N)) # new

   # Waning challenge immunity (from infection challenge)
   %>% vec_rate(
     "R" %_% vax_cat,
     "S" %_% vax_cat,
     chal_waning
   ) # new

   # definitely need a better syntax here
   %>% add_outflow(
     ".+",
     "^" %+% alt_group(non_accum) %_% alt_group(vax_cat))

   # Update parameters for use with the linearized model
   # -- confirmed correct (TODO: check if params are missing? variant-related?)
   %>% update_linearized_params('^N$', 1) # scale population to 1
   %>% update_linearized_params('^E0$', 1e-5)
   %>% update_linearized_params('^vax_doses_per_day$', 0)
   %>% update_linearized_params('^vax_response_rate$', 0)
   %>% update_linearized_params('^vax_response_rate_R$', 0)

   # Set the disease-free equilibrium of the linearized model
   %>% update_disease_free_state('S_unvax', 'S0')

   # Perturb the disease-free equilibrium of the linearized model
   %>% update_disease_free_state('E_unvax', 'E0')

   # Define outflows for the linearized model
   # -- confirmed that this is producing the correct indices
   %>% add_linearized_outflow("^S", "^S") # S_pos, S_pos
   %>% add_linearized_outflow(
     "^" %+% alt_group(non_accum_non_S) %_% alt_group(vax_cat), # notS_pos
     "^" %+% alt_group(non_accum)       %_% alt_group(vax_cat)) # p_states

   # Define state mappings used to put the initial state values in
   # the correct positions of the initial state vector
   %>% add_state_mappings(

     # regular expression to find states to drop before computing
     # the eigenvector of the linearized system
     # -- generated indices are correct
     eigen_drop_pattern = '^(X|V)',

     # regular expression to find states to drop from the eigenvector
     # before distributing individuals among infected compartments
     # -- generated indices are correct
     infected_drop_pattern = '^(S|D|R)',

     # regular expression to find states in the initial population
     # of susceptibles
     initial_susceptible_pattern = '^S_unvax$'
   )

   # Set the total number of individuals and the total number of
   # infected individuals in the initial state vector
   %>% initial_population(total = 'N', infected = 'E0')

   %>% update_tmb_indices
   %>% update_initial_state(silent = FALSE)
)

(simulate_state_vector(model, add_dates = TRUE)
  %>% tidyr::pivot_longer(-Date)
  %>% tidyr::separate(name, into = c("epi", "vax"), sep = "_")
  #%>% filter(vax == "vaxprotect3")
  #%>% filter(epi == "S")
  %>% ggplot
   +  facet_grid(epi ~ vax, scales = 'free')
   +  geom_line(aes(Date, value))
)

sims = run_sim(
  params = params,
  state = model$state,
  start_date = model$start_date,
  end_date = model$end_date,
  condense = TRUE,
  flexmodel = model
)

fitdat = (sims
  %>% mutate(report = round(report))
  %>% select(date, report)
  %>% tidyr::pivot_longer(-date, names_to = "var")
  %>% filter(!is.na(value))
)
# shows a little waning immunity bump after
# the main peak
ggplot(fitdat) + geom_line(aes(date, value))

opt_pars = list(params = list(beta0 = 1))

fitted_model <- calibrate(
  base_params  = params
  , data       = fitdat
  , opt_pars   = opt_pars
  , time_args  = list(params_timevar = NULL)
  , sim_args   = list(ndt = 1,
                      step_args = list(do_hazard = TRUE),
                      flexmodel = model)
  , use_DEoptim = FALSE
  , debug = TRUE
  , debug_plot = FALSE
)

# recovers beta0
fitted_model$mle2@coef
