library(McMasterPandemic)

## load("ontario_calibration.RData")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)
## use hospitalization data only
print(unique(ont_all_sub$var))
ont_recent_hosp <- na.omit(dplyr::filter(ont_all_sub, var=="H"))  ## var %in% c("H","ICU","death")  OR var != "report"
print(opt_pars)  ## original parameter settings
opt_pars_2brk <- opt_pars
opt_pars_2brk$log_rel_beta0 <- rep(-1,2)  ## only two breakpoints (hosp data doesn't even start until after brk 1)
bd2 <- bd[-1]  ## drop first breakpoint
## don't try to estimate fraction mild (= !hospital, mu) or fraction in acute care (= !ICU, phi1)
opt_pars_2brk$params <- opt_pars_2brk$params[c("log_E0","log_beta0")]
ont_cal2 <- update(ont_cal1,  data=ont_recent_hosp, opt_pars=opt_pars_2brk, break_dates=bd2)

# rdsave("ont_cal2", "opt_pars_2brk", "bd2")
