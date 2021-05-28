library(McMasterPandemic)

load("data/ONcalib.rda")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)

opt_pars_2brks <- opt_pars
opt_pars_2brks$logit_rel_beta0 <- rep(-1,2)  ## only two breakpoints (hosp data doesn't even start until after brk 1)
bd2 <- bd[-1]  ## drop first breakpoint
ont_cal_2brks <- update(ont_cal1
                     ,  opt_pars=opt_pars_2brks
                     , time_args=list(break_dates=bd2)
                       )

save("ont_cal_2brks", file=sprintf("data/ONcalib_2brks.rda",
                              format(Sys.time(),"%Y%b%d")))

