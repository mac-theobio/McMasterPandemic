library(McMasterPandemic)
library(ggplot2)
library(tidyverse)
library(anytime)
library(bbmle)
library(parallel)

use_true_start <- TRUE
cut_dates <- TRUE
nsim <- 24
start_date_offset <- 7
priors <- list(~dnorm(qlogis(rel_beta0[1]),mean=1,sd=1.5))
## +/- 2 sigma = 2% to 88% drop

## true value below is -1.38 (near the lower bound ... pnorm() = 0.056
## FIXME option setting for rel_break on log vs logit scale
## mclapply will use getOption("mc.cores",2) by default
##   set options(mc.cores=...) upstream (e.g. in Rprofile)
## but wrapR seems to ignore this?
options(mc.cores=6)

## setup 

params <- fix_pars(read_params("ICU1.csv"), target=c(R0=3, Gbar=6))
## params[["beta0"]] <- 2
params[["beta0"]] <- 0.9   ## slightly rounded for convenience
params[["obs_disp"]] <- 100 ## BMB: less noise

start_date <- anydate("2020-01-01")
end_date <- anydate("2020-03-08") ## BMB: don't run as long

cutoff_start <- anydate("2020-02-01")
cutoff_end <- anydate("2020-03-08")

break1 <- "2020-02-15"
bd <- anydate(break1)
rel_break1 <- 0.2

if (use_true_start) {
    opt_pars <- list(
        params=c(log_E0=log(params[["E0"]]), log_beta0=log(params[["beta0"]]))
      , logit_rel_beta0 = qlogis(rel_break1)
      , log_nb_disp = log(params[["obs_disp"]])
    )
}  else {
    opt_pars <- list(
        ## these params go to run_sim
        params=c(log_E0=4, log_beta0=-1)
      , logit_rel_beta0 = rep(-1, length(bd))
      , log_nb_disp=0
    )
}

params[["N"]] <- 1e7


sim_cali <- function(x){
    set.seed(x)
    cat("seed ",x,"\n")
    suppressWarnings(
        sim1break <- run_sim_break(params
      , start_date=start_date
      , end_date=end_date
      , break_dates = bd
      , rel_beta0= rel_break1
      , stoch = c(obs = TRUE, proc=FALSE)
        )
        )

   ## plot(sim1break,log=TRUE)
    simdat <- (sim1break
        %>% pivot()
        %>% filter(var %in% c("report"))
    )

    if (cut_dates) {
        dd <- filter(simdat,between(date,cutoff_start, cutoff_end))
    }

    g1 <- calibrate(data=dd
                  , base_params=params
                  , opt_pars = opt_pars
                  , break_dates = bd
                  ## , debug_plot=TRUE
                  , start_date_offset = start_date_offset
                  , priors = priors
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
    return(list(simdat=simdat,fit=g1,pars=res_dat,pred=pp, fullsim=sim1break))
}

## mclapply()
res <- mclapply(seq(nsim), sim_cali)
# rdsave("bd","params","res", "rel_break1")
