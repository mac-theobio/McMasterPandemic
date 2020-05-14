library(testthat)
library(McMasterPandemic)
do_slow <- FALSE
library(dplyr)
params <- fix_pars(read_params("ICU1.csv"))
opt_pars <- list(params=c(log_E0=4, log_beta0=-1,
         log_mu=log(params[["mu"]]), logit_phi1=qlogis(params[["phi1"]])),
                                   logit_rel_beta0=c(-1,-1),
                                    log_nb_disp=NULL)
dd <- (ont_all %>% trans_state_vars() %>% filter(var %in% c("report", "death", "H")))
if (do_slow) {
    cal1 <- calibrate(data=dd, base_params=params, opt_pars=opt_pars)
}
