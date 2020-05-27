library(McMasterPandemic)
library(testthat)
library(parallel)

context("DEoptim")

library(dplyr)
params <- fix_pars(read_params("ICU1.csv"))
opt_pars <- list(params=c(log_E0=4, log_beta0=-1,
         log_mu=log(params[["mu"]]), logit_phi1=qlogis(params[["phi1"]])),
                                   logit_rel_beta0=c(-1,-1),
                                    log_nb_disp=NULL)
dd <- (ont_all %>% trans_state_vars() %>% filter(var %in% c("report", "death", "H")))

suppressWarnings(cal1_DE <- calibrate(data=dd, base_params=params, opt_pars=opt_pars,
                     use_DEoptim=TRUE,
                     DE_cores=1,
                     DE_args=list(control=list(itermax=5,trace=FALSE)))
                 )

de <- attr(cal1_DE,"de")

test_that("DE components are present", {
    expect_is(attr(cal1_DE,"de_time"),"proc_time")
    expect_is(de,"DEoptim")
    expect(!is.null(de$member$Sigma),failure_message="Sigma component missing")
})
## test parallel?
