##
## possible params to set/command-line arguments:
##  Gbar = 6
##  which vars to fit?
library(McMasterPandemic)
library(tidyverse)
library(anytime)

## 
keep_vars <- c("H","ICU","d","report")
ont_recent_sub <- (ont_recent
    %>% mutate_at("var",trans_state_vars)
    %>% filter(var %in% keep_vars)
)

ont_all_sub <- (ont_all
    %>% mutate_at("var",trans_state_vars)
    %>% filter(var %in% keep_vars)
)

## unique(ont_recent_sub$var)

## adjust mean GI
params <- fix_pars(read_params("ICU1.csv")
    , target=c(Gbar=6)
    , pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))
)
params[["N"]] <- 19.5e6  ## reset pop to Ontario

## breakpoints
schoolClose <- "2020-Mar-17"
countryClose <- "2020-Mar-23"
socialClose <- "2020-Mar-28"

bd <- anydate(c(schoolClose,countryClose,socialClose))
## print(bd)

opt_pars <- list(
    ## these params go to run_sim
    params=c(log_E0=4
             , log_beta0=-1
             ## fraction of mild (non-hosp) cases
             , log_mu=log(params[["mu"]])
             ## fraction of incidence reported
             ## logit_c_prop=qlogis(params[["c_prop"]]),
             ## fraction of hosp to acute (non-ICU)
             , logit_phi1=qlogis(params[["phi1"]])
             ## fraction of ICU cases dying
             ## logit_phi2=qlogis(params[["phi2"]])
             ),
    log_rel_beta0 = rep(-1, length(bd)),
    log_nb_disp=0)

t_ont_cal1 <- system.time(ont_cal1 <- calibrate(data=ont_all_sub
    , base_params=params
    , optim_args=list(control=list(maxit=10000),hessian=TRUE)
    , opt_pars = opt_pars,
    , break_dates = bd
    ## , debug=TRUE
      )
      ) ## system.time

# rdsave("t_ont_cal1","opt_pars","ont_cal1", "bd","ont_recent_sub","params","keep_vars", "ont_all_sub")



