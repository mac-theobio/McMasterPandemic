library(McMasterPandemic)

## load("ontario_calibration.RData")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)
print(unique(ont_all_sub$var))
print(opt_pars)  ## original parameter settings
ont_noICU <- dplyr::filter(ont_all_sub, var != "ICU")
opt_pars_2brks <- opt_pars
opt_pars_2brks$logit_rel_beta0 <- rep(-1,2)  ## only two breakpoints (hosp data doesn't even start until after brk 1)
bd2 <- bd[-1]  ## drop first breakpoint
priors <- list(~dlnorm(rel_beta0[1],meanlog=-1,sd=0.5))
ont_cal_noICU_2brks <- update(ont_cal1
                            , opt_pars=opt_pars_2brks
                            , time_args=list(break_dates=bd2)
                            , data=ont_noICU
                            , priors=priors
                              )

# rdsave("ont_cal_noICU_2brks")
