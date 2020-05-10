library(McMasterPandemic)

## load("ontario_calibration.RData")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)
print(unique(ont_all_sub$var))
ont_noICU <- dplyr::filter(ont_all_sub, var != "ICU")
print(opt_pars)  ## original parameter settings
ont_cal_noICU <- update(ont_cal1,  data=ont_noICU)

# rdsave("ont_cal_noICU")
