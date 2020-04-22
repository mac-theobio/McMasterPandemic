library(McMasterPandemic)
library(tidyverse)
library(anytime)
library(bbmle)
library(parallel)

use_true_start <- TRUE
nsim <- 50
options(mc.cores=6)

## setup 

params <- fix_pars(read_params("ICU1.csv"), target=c(R0=3, Gbar=6))
## params[["beta0"]] <- 2
params[["beta0"]] <- 0.9   ## slightly rounded for convenience
params[["obs_disp"]] <- 100 ## BMB: less noise
params[["N"]] <- 1e7
summary(params)  ## v. high R0 (6.7)
start_date <- anydate("2020-01-01")
end_date <- anydate("2020-03-31") ## BMB: don't run as long

## Start with true parameters
##    optimization breaks if factor of true value is > about 1.5 ... ?
opt_pars <- list(
    params=c(log_beta0=log(params[["beta0"]]*1.2))
  , log_nb_disp = log(params[["obs_disp"]])
)

sim_cali <- function(seed) {
    cat(seed,"\n")
    set.seed(seed)
    simdat <- run_sim(params
                    , start_date=start_date
                    , end_date=end_date
                    , stoch = c(obs = TRUE, proc=FALSE)
                      )
    ## plot(simdat, log=TRUE)
    simdat <- (simdat
	%>% condense()
	%>% pivot()
	%>% filter(var %in% c("report"))
	%>% filter(!is.na(value)) 
	%>% mutate(value = round(value))
    )

    ## print(params)
    ## print(opt_pars)
    g1 <- calibrate(data=simdat, base_params=params
                  , start_date = start_date
                  , opt_pars = opt_pars
                  , break_dates = NULL
                    ## , debug_plot=TRUE
                    ## , debug=TRUE
                    ## , mle2_args=list(browse_obj=TRUE)
                    )

    print(bbmle::coef(g1$mle2))

    pp <- predict(g1)

    res_dat <- data.frame(bbmle::confint(g1$mle2, method="quad", level=0.95)
                        , estimate = bbmle::coef(g1$mle2)
                        , seed = seed
                        , pars = names(g1$mle2@coef)
                          )
    return(list(simdat=simdat,fit=g1,pars=res_dat,pred=pp, fullsim=simdat))
}

res <- mclapply(seq(nsim), sim_cali)
# rdsave(res,params)
