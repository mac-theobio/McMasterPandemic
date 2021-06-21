library(profvis)
library(McMasterPandemic)
library(anytime)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
library(tidyverse)
library(zoo)

for(file in list.files(path='../R')){
  source(paste('../R/', file, sep=''))
}

## get base params
params <- read_params("PHAC.csv")

## vaxify params
params <- expand_params_vax(
  params,
  model_type = "twodose",
  vax_doses_per_day = 1,
  vax_prop_first_dose = 0.9,
  vax_efficacy_dose1 = 0.6,
  vax_efficacy_dose2 = 0.9,
  vax_avg_response_time = 14,
  vax_avg_response_time_R = 7,
  vax_alpha_dose1 = 0.6,
  vax_alpha_dose2 = 0.9,
  vax_mu_dose1 = 1,
  vax_mu_dose2 = 1
)

## make state (will be vaxified because input params are vaxified)
state <- make_state(params = params)

vaxify <- profvis(
  {
    for (i in c(1:10)){
    ## run sim
      run_sim(params, state)
    }
  },
  prof_output = "vaxify.out"
)
vaxify_summary <- summaryRprof("vaxify.out")

# Note: in afuns, expsim is being used by default
write.csv(vaxify_summary$`by.self`, 'vaxify.csv')
htmlwidgets::saveWidget(vaxify, "vaxify_profile.html")

