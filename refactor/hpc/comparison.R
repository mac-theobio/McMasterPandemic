
library(McMasterPandemic)
library(tidyverse)


keep_vars <- c("H", "ICU", "death", "report")
## data since 15 March

ont_all_sub <- (ont_all
%>% mutate_at("var", trans_state_vars)
  %>% filter(var %in% keep_vars)
)

## adjust mean GI
params <- fix_pars(read_params("ICU1.csv"),
  target = c(Gbar = 6),
  pars_adj = list(c("sigma", "gamma_s", "gamma_m", "gamma_a"))
)
params[["N"]] <- 14.57e6 ## reset pop to Ontario



bd <- c("2020-03-17", "2020-03-23", "2020-03-28")


## This is copied over from the old script
opt_pars <- list(
  ## these params are part of the main parameter vector: go to run_sim()
  params = c(
    log_E0 = 4 ## initial exposed
    , log_beta0 = -1 ## initial baseline transmission
    ## fraction of mild (non-hosp) cases
    , log_mu = log(params[["mu"]])
    ## fraction of incidence reported
    ## logit_c_prop=qlogis(params[["c_prop"]]),
    ## fraction of hosp to acute (non-ICU)
    , logit_phi1 = qlogis(params[["phi1"]])
    ## fraction of ICU cases dying
    ## logit_phi2=qlogis(params[["phi2"]])
  ),
  ## changes in beta at breakpoints
  logit_rel_beta0 = rep(-1, length(bd)),
  ## NB dispersion
  log_nb_disp = NULL
)


timevar_df <- data.frame(Date = bd,
                         Symbol = c("beta0", "beta0", "beta0"),
                         Relative_value = c(NA, NA, NA))

do_hazard <- calibrate(
  data = ont_all_sub,
  base_params = params,
  opt_pars = opt_pars,
  time_args = list(params_timevar = timevar_df),
  sim_args = list(step_args = list(do_hazard = TRUE)),
  use_DEoptim = TRUE,
  DE_cores = 7 #future::availableCores(methods='mc.cores')

)

do_hazard$de_time
do_hazard$mle2_time

save("do_hazard",  file="ONcalib_do_hazard.rda")

