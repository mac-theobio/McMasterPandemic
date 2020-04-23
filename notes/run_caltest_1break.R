library(McMasterPandemic)
library(ggplot2)
library(tidyverse)
library(anytime)
library(parallel)

use_true_start <- TRUE
nsim <- 1
## FIXME option setting for rel_break on log vs logit scale
options(mc.cores=2)

## setup 

params <- fix_pars(read_params("ICU1.csv"))
params[["beta0"]] <- 2
params[["obs_disp"]] <- 100 ## BMB: less noise

start_date <- anydate("2020-01-01")
end_date <- anydate("2020-03-08") ## BMB: don't run as long

break1 <- "2020-02-01"
bd <- anydate(break1)
rel_break1 <- 0.2

opt_pars <- list(
   params=c(log_beta0=log(params[["beta0"]]))
   , logit_rel_beta0 = qlogis(rel_break1)
   , log_nb_disp = log(params[["obs_disp"]])
)

params[["N"]] <- 1e7

sim1break <- forecast_sim(params
   , opt_pars = opt_pars
	, base_params = params
	, start_date=start_date
   , end_date=end_date
   , break_dates = bd
   # , rel_beta0= rel_break1
   , stoch = c(obs = TRUE, proc=FALSE)
)

class(sim1break) <- c(class(sim1break),"pansim")

## plot(sim1break,log=TRUE)
simdat <- (sim1break
	# %>% filter(var=="report")
   # %>% condense()
 	# %>% pivot()
 	%>% filter(var %in% c("report"))
)

    
 ## Need to round value because of negative binomial fit
 dd <- (simdat 
 	%>% filter(!is.na(value)) 
   %>% mutate(value = round(value))
 )

## WORKS FINE if we don't limit the dates?
## OVERSHOOTS TERRIBLY if we do
if (cut_dates) {
	dd <- filter(dd,between(date,cutoff_start, cutoff_end))
 }

 g1 <- calibrate(data=dd
 	, base_params=params
 	, start_date = start_date
	, opt_pars = opt_pars
 	, break_dates = bd
 	, debug_plot=TRUE
 )

    
 pp <- predict(g1)

    
res_dat <- data.frame(bbmle::confint(g1$mle2, method="quad", level=0.95)
	, estimate = bbmle::coef(g1$mle2)
   , seed = x
   , pars = names(g1$mle2@coef)
)

# return(list(simdat=simdat,fit=g1,pars=res_dat,pred=pp, fullsim=sim1break))


## mclapply()
# res <- mclapply(seq(nsim), sim_cali)
# rdsave("bd","params","res", "rel_break1")
