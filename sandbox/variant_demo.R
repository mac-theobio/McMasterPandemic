devtools::load_all("../McMasterPandemic")
rel_path <- "../macpan_vaccine"

## sim setup (make params and state)
base_params <- read_params(file.path(rel_path,
                                     "params",
                                     "model_params_ON.csv"))

base_params <- fix_pars(base_params,target=c(R0=1.1))

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
dd <- data.frame(Date = seq.Date(as.Date('2020-03-01'),as.Date("2020-07-01"),by=1)
				, Symbol = "variant_prop"
				, Relative_value = 1
)
dd2 <- (dd
	%>% mutate(Relative_value = ifelse(Date>=as.Date("2020-04-01"),100,Relative_value))
)


sim <- run_sim(params, state, params_timevar=dd,start_date="2020-03-01",end_date="2020-07-01")
sim2 <- run_sim(params, state, params_timevar=dd2,start_date="2020-03-01",end_date="2020-07-01")

(ggplot()
	+ geom_line(data=sim,aes(x=date,y=report),color="black")
	+ geom_line(data=sim2,aes(x=date,y=report),color="red")
)
