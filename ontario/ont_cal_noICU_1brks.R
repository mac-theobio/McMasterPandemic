library(McMasterPandemic)

## load("ontario_calibration.RData")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)
print(unique(ont_all_sub$var))
print(opt_pars)  ## original parameter settings
ont_noICU <- dplyr::filter(ont_all_sub, var != "ICU")
opt_pars_1brks <- opt_pars
opt_pars_1brks$logit_rel_beta0 <- -1
bd1 <- as.Date("2020-03-20")
ont_cal_noICU_1brks <- update(ont_cal1,  opt_pars=opt_pars_1brks,
                              time_args=list(break_dates=bd1), data=ont_noICU)

# rdsave("ont_cal_noICU_1brks")
