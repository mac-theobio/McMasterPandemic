library(McMasterPandemic)

load("ONcalib_2brks_2021May11.rda")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)

load("../../inst/testdata/ONcalib_2020Jun01.rda")
old_ont_cal_2brks <- ont_cal_2brks
#-----------------------------------------------------------------------------------------
# None of this actually matters, we copy over the $mle2 object from the old ont_cal_1
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
opt_pars_2brks <- opt_pars
opt_pars_2brks$logit_rel_beta0 <- rep(-1,2)  ## only two breakpoints (hosp data doesn't even start until after brk 1)
bd2 <- bd[-1]  ## drop first breakpoint
ont_cal_2brks <- update(ont_cal1
                     ,  opt_pars=old_ont_cal_2brks$mle2@data$opt_pars
                     , time_args=list(break_dates=bd2)
                       )
#---------------------------------------------------------
ont_cal_2brks$mle2@fullcoef <- old_ont_cal_2brks$mle2@fullcoef
ont_cal_2brks$mle2@coef <- old_ont_cal_2brks$mle2@coef
ont_cal_2brks$mle2@vcov <- old_ont_cal_2brks$mle2@vcov
ont_cal_2brks$mle2@min <- old_ont_cal_2brks$mle2@min
ont_cal_2brks$mle2@details <- old_ont_cal_2brks$mle2@details
ont_cal_2brks$mle2@minuslogl <- old_ont_cal_2brks$mle2@minuslogl

save("ont_cal_2brks", file=sprintf("ONcalib_2brks_%s.rda",
                                   format(Sys.time(),"%Y%b%d")))

print(ont_cal_2brks)
print(old_ont_cal_2brks)


