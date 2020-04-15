##
## possible params to set/command-line arguments:
##  Gbar = 6
##  which vars to fit?

library(McMasterPandemic)
library(tidyverse)
library(anytime)

## 
load("notes/ontario_clean.RData")
keep_vars <- c("H","ICU","d","report", "incidence","newTests")
ont_recent_sub <- (ont_recent
    %>% mutate_at("var",trans_state_vars)
    %>% filter(var %in% keep_vars)
)

unique(ont_recent_sub$var)

## adjust mean GI
params <- fix_pars(read_params("ICU1.csv")
    , target=c(Gbar=6)
    , pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))
)
params[["N"]] <- 19.5e6  ## reset pop to Ontario

## breakpoints
schoolClose <- "17-Mar-2020"
countryClose <- "23-Mar-2020"
socialClose <- "28-Mar-2020"

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

t1 <- system.time(g1 <- calibrate(data=ont_recent_sub
    , base_params=params
    , optim_args=list(control=list(maxit=10000),hessian=TRUE)
    , opt_pars = opt_pars,
    , break_dates = bd
    ## , debug=TRUE
      )
      ) ## system.time

ont_recent_hosp <- filter(ont_recent_sub, var=="H")

## FIXME: break this out into a separate file? (risk of atomization/confusion?)
g2 <- update(g1,  data=ont_recent_hosp)
# rdsave("t1","opt_pars","g1", "g2", "bd","ont_recent_sub","params","keep_vars")



