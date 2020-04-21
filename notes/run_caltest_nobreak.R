library(McMasterPandemic)
library(ggplot2)
library(tidyverse)
library(anytime)
library(parallel)

use_true_start <- TRUE
nsim <- 100
options(mc.cores=2)

## setup 

params <- fix_pars(read_params("ICU1.csv"))
params[["obs_disp"]] <- 100 ## BMB: less noise

start_date <- anydate("2020-01-01")
end_date <- anydate("2020-03-31") ## BMB: don't run as long

break1 <- NULL
bd <- anydate(c(break1))

if (use_true_start) {
   opt_pars <- list(
      params=c(log_beta0=log(params[["beta0"]]))
      , log_nb_disp = log(params[["obs_disp"]])
   )
}  else {
   opt_pars <- list(
      ## these params go to run_sim
      params=c(log_E0=4, log_beta0=-1)
      , log_rel_beta0 = rep(-1, length(bd))
      , log_nb_disp=0
   )
}

params[["N"]] <- 1e7


sim_cali <- function(x){
   set.seed(x)
   cat("seed ",x,"\n")
   suppressWarnings(simdat <- run_sim(params
      , start_date=start_date
      , end_date=end_date
      , stoch = c(obs = TRUE, proc=FALSE)
   )
   )
   
   ## plot(sim1break,log=TRUE)
   simdat <- (simdat
              %>% condense()
              %>% pivot()
              %>% filter(var %in% c("report"))
   )
   
   
   ## Need to round value because of negative binomial fit
   dd <- (simdat 
          %>% filter(!is.na(value)) 
          %>% mutate(value = round(value))
          # %>% filter(between(date,cutoff_start, cutoff_end))
   )
   
   g1 <- calibrate(data=dd, base_params=params
                   , start_date = start_date
                   , opt_pars = opt_pars
                   , break_dates = NULL
                   , debug_plot=TRUE
   )
   
   pp <- predict(g1)
   
   if (FALSE) {
      plot(sim1break, keep_state=c("incidence","report"),log=TRUE) +
         geom_line(data=filter(pp,var=="report"),lty=2)
      
      (ggplot(simdat, aes(date,value))
         + scale_y_log10()
         + geom_line()
         + geom_line(data=filter(pp,var=="report"),lty=2)
      )
   }
   
   res_dat <- data.frame(bbmle::confint(g1$mle2, method="quad", level=0.95)
                         , estimate = bbmle::coef(g1$mle2)
                         , seed = x
                         , pars = names(g1$mle2@coef)
   )
   return(list(simdat=simdat,fit=g1,pars=res_dat,pred=pp))
}

## mclapply()
res <- mclapply(1:nsim, sim_cali)
save("res", file="run_caltest_nobreak.RData")
