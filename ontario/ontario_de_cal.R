library(McMasterPandemic)
library(bbmle)
library(DEoptim)
library(parallel)

set.seed(101)
fit1 <- ont_cal_noICU_2brks_prior   ## abbreviation
ff <- fit1$forecast_args
ff$fixed_pars <- NULL

lwr <- c(params.log_E0=1,params.log_beta0=-1,params.log_mu=-1,
         params.logit_phi1=-1,logit_rel_beta01=-1,logit_rel_beta02=-1,
         log_nb_disp.H=-1,log.nb_disp.report=-1,log.nb_disp_death=-1)

upr <- c(params.log_E0=5,params.log_beta0=1,params.log_mu=1,
         params.logit_phi1=1,logit_rel_beta01=4,logit_rel_beta02=4,
         log_nb_disp.H=5,log.nb_disp.report=5,log.nb_disp_death=5)


## FIXME: inherit from somewhere sensible
cl <- makeCluster(5)
de_arglist <- c(list(fn=mle_fun, lower=lwr, upper=upr,
               debug_plot=TRUE,
               data=fit1$mle2@data$data,
               control=DEoptim.control(storepopfrom=1, parallelType=1,
                                       cluster=cl,
                                       parVar=c("fit1","lwr","upr"),
                                       packages=list("McMasterPandemic","bbmle"))),
               ff)

de_time <- system.time(de_cal1 <- do.call(DEoptim,de_arglist))
## stopCluster(cl)

# rdsave("de_time", "de_cal1")
