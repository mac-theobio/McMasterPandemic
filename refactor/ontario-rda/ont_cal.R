
library(McMasterPandemic)
library(tidyverse)

load("../../inst/testdata/ONcalib_2020Jun01.rda")
old_ont_cal1 <- ont_cal1

#-----------------------------------------------------------------------------------------
# None of this actually matters, we copy over the $mle2 object from the old ont_cal_1
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
t_ont_cal1 <- system.time(ont_cal1 <- calibrate(data=old_ont_cal1$mle2@data$data#ont_all_sub
    , base_params= old_ont_cal1$mle2@data$base_params
    , opt_pars =  old_ont_cal1$mle2@data$opt_pars
    , time_args=list(break_dates = bd)
      )
      ) ## system.time

print(ont_cal1)
#---------------------------------------------------------------------------------------
ont_cal1$mle2@fullcoef <- old_ont_cal1$mle2@fullcoef
ont_cal1$mle2@coef <- old_ont_cal1$mle2@coef
ont_cal1$mle2@vcov <- old_ont_cal1$mle2@vcov
ont_cal1$mle2@min <- old_ont_cal1$mle2@min
ont_cal1$mle2@details <- old_ont_cal1$mle2@details
ont_cal1$mle2@minuslogl <- old_ont_cal1$mle2@minuslogl


save("ont_cal1", "bd", "ont_all_sub", file=sprintf("ONcalib_%s.rda",
                              format(Sys.time(),"%Y%b%d")))
print(ont_cal1)
print(old_ont_cal1)
