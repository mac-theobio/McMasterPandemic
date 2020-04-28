library(McMasterPandemic)
library(bbmle)
library(DEoptim)

L <- load("../ontario/ontario_calibration.RData")
L <- load("../ontario/ontario_calibration_noICU_2brks_prior.RData")
L <- load("alt_cal.RData")
fit1 <- ont_cal_noICU_2brks_prior   ## abbreviation

## refit **from same params**
opt_pars_new  <- opt_pars_2brks
opt_pars_2brks$log_nb_disp <- setNames(rep(0,3),c("H","report","death"))
pars2 <- relist(coef(fit1$mle2), opt_pars_2brks)
fit2 <- update(ont_cal_noICU_2brks_prior, opt_pars=pars2)


## fit from original params, lower tolerance
fit3 <- update(ont_cal_noICU_2brks_prior, mle2_control=list(maxit=1e4,control=list(reltol=1e-12)))
fit4 <- update(ont_cal_noICU_2brks_prior, mle2_control=list(maxit=1e5,control=list(reltol=1e-16)))

