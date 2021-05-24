##
## possible params to set/command-line arguments:
##  Gbar = 6
##  which vars to fit?
library(McMasterPandemic)
library(tidyverse)

## clean Ontario data are automatically loaded
## load("ontario_clean.RData")
## translate variable names to internally used values
## drop unused variables
keep_vars <- c("H","ICU","death","report")
## data since 15 March
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
params[["N"]] <- 14.57e6  ## reset pop to Ontario

## breakpoints
schoolClose <- "2020-03-17"
countryClose <- "2020-03-23"
socialClose <- "2020-03-28"

bd <- c(schoolClose,countryClose,socialClose)
## print(bd)

## choose parameters to fit: starting values
opt_pars <- list(
    ## these params are part of the main parameter vector: go to run_sim()
    params=c(log_E0=4      ## initial exposed
           , log_beta0=-1  ## initial baseline transmission
             ## fraction of mild (non-hosp) cases
           , log_mu=log(params[["mu"]])
             ## fraction of incidence reported
             ## logit_c_prop=qlogis(params[["c_prop"]]),
             ## fraction of hosp to acute (non-ICU)
             , logit_phi1=qlogis(params[["phi1"]])
             ## fraction of ICU cases dying
             ## logit_phi2=qlogis(params[["phi2"]])
             ),
    ## changes in beta at breakpoints
    logit_rel_beta0 = rep(-1, length(bd)),
    ## NB dispersion
    log_nb_disp=NULL)

print(comb_sub)

## do the calibration
t_ont_cal1 <- system.time(ont_cal1 <- calibrate(data=ont_all_sub
    , base_params=params
    , opt_pars = opt_pars
    , time_args=list(break_dates = bd)
      )
      ) ## system.time

save("ont_cal1", file=sprintf("ONcalib_%s.rda",
                              format(Sys.time(),"%Y%b%d")))
# rdsave("t_ont_cal1","opt_pars","ont_cal1", "bd","ont_recent_sub","params","keep_vars", "ont_all_sub")
