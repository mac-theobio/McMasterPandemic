## reinstall latest version first
# library(McMasterPandemic) ## not sure why but when i reinstall and reload, update.fit_pansim isn't found, but load_all does the trick...
devtools::load_all()

## load sample calibration that used only the vaxified model (not ageified)
load("sample_vax_calib.Rdata")

## set up ageified beta to slot into vaxified calibration
log_beta0 <- c(-0.9990100, -0.6810834)
newcoef <- list(params.log_beta0 = log_beta0)

## attempt to update some coefs
calib$fit <- update(calib$fit, newcoef = newcoef)

## change format of opt_pars so restore works in coef.fit_pansim
opt_pars <- calib$fit$forecast_args$opt_pars
opt_pars$params <- unlist(list(log_beta0 = log_beta0))
calib$fit$forecast_args$opt_pars <- opt_pars

## the update is saved in the attributes of the mle2 object
print(calib$fit$mle2@coef)
print(attr(calib$fit$mle2, "coef"))

## but not when accessed via coef---the bbmle::mle2 method :(
print(coef(calib$fit$mle2))

## so when we try to access the updated coefs using the coef.fit_pansim method, it doesn't work because internally it's calling the bbmle::mle2 coef method
print(coef(calib$fit, "fitted"))

## i could change coef(object$mle2) to object$mle2@coef within coef.fit_pansim, but i'm worried that by not actually editing the bbmle::mle2 object's coefs, i haven't actually fully updated the calib$fit object's coefs, which may cause invisible trouble

