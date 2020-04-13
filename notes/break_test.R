library(McMasterPandemic)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
source("notes/calib_funs.R")
params <- fix_pars(read_params("ICU1.csv"), target=c(Gbar=6),u_interval=c(-1,1),
                   pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a")))
params[["N"]] <- 19.5e6  ## reset pop to Ontario

opt_pars <- list(log_E0=4, log_beta0=-1, log_rel_beta0=c(-1,-1), log_nb_disp=0)
p <- c(log_E0=2,log_beta0=-0.24,log_relbeta01=-0.5,log_rel_beta02=-8, log_nb_disp=4)
bd <- ldmy(c("23-Mar-2020","30-Mar-2020"))
f1 <- forecast_sim(p, opt_pars, base_params=params,
             start_date="29-02-2020", end_date="01-06-2020",
             break_dates=bd)
ggplot(f1,aes(date,value,colour=var)) + geom_line() + geom_vline(xintercept=bd,lty=2)
