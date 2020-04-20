library(McMasterPandemic)
library(ggplot2)
library(tidyverse)
library(anytime)

params <- fix_pars(read_params("ICU1.csv"))


start_date <- anydate("2020-01-01")
end_date <- anydate("2020-04-01")

change1 <- "2020-01-22"
bd <- anydate(c(change1))
params[["N"]] <- 1e7

opt_pars <- list(
   params=c(log_E0=4, log_beta0=-0.5, log_mu=log(params[["mu"]]))
   , log_rel_beta0 = rep(-4, length(bd))
   # , log_nb_disp=0
)

pp <- invlink_trans(restore(opt_pars, opt_pars))

simdat <- forecast_sim(opt_pars
   , opt_pars = opt_pars
   , base_params = params
   , start_date=start_date
   ,end_date=end_date
   ,break_dates = bd)


simreport <- (simdat 
   # %>% filter(var == "S")
   %>% filter(var %in% c("I","report"))
)

print(ggplot(simreport,aes(x=date,y=value,color=var))
   + geom_line()
   + geom_vline(xintercept = bd)
)

