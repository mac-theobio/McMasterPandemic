devtools::load_all("../McMasterPandemic")
rel_path <- "../macpan_vaccine"

## sim setup (make params and state)
base_params <- read_params(file.path(rel_path,
                                     "params",
                                     "model_params_ON.csv"))

vax_params <- expand_params_vax(
  base_params,
  model_type = "twodose",
  vax_doses_per_day = 1,
  vax_prop_first_dose = 1,
  vax_efficacy_dose1 = 0.6,
  vax_efficacy_dose2 = 0.9,
  vax_avg_response_time = 14,
  vax_avg_response_time_R = 7,
  vax_alpha_dose1 = 0.6,
  vax_alpha_dose2 = 0.9,
  vax_mu_dose1 = 1,
  vax_mu_dose2 = 1
)

params <- expand_params_variant(
  vax_params,
  variant_prop = 0.01,
  variant_advantage = 1.5,
  variant_vax_efficacy_dose1 = 0.3,
  variant_vax_efficacy_dose2 = 0.8
)

## by default, states for vaxified parameters put everyone in the unvax layer first
state <- make_state(params = params)

## run sim
sim <- run_sim(params, state)
