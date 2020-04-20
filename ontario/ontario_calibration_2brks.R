library(McMasterPandemic)

## load("ontario_calibration.RData")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)
print(unique(ont_all_sub$var))
print(opt_pars)  ## original parameter settings
opt_pars_2brks <- opt_pars
opt_pars_2brks$log_rel_beta0 <- rep(-1,2)  ## only two breakpoints (hosp data doesn't even start until after brk 1)
bd2 <- bd[-1]  ## drop first breakpoint
ont_cal_2brks <- update(ont_cal1,  opt_pars=opt_pars_2brks, break_dates=bd2)

# rdsave("ont_cal_2brks", "opt_pars_2brks", "bd2")
