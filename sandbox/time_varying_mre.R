library(McMasterPandemic)
library(zoo)
library(tidyverse)
library(shellpipes)

## Copying example from run_sim and modifying it

##1. testing if Relative_value=1 and non-timevar run_sim are the same (done)
## MLi: This means I can always construct my own timevar for any time-varying simulation 
## Move this to testhat?!? 

params <- read_params("ICU1.csv")
## paramsS <- update(params,c(proc_disp=0.1,obs_disp=100))
## paramsSz <- update(paramsS, zeta=5)

state <- make_state(params=params)
startdate <- as.Date("2020-01-01")
enddate <- as.Date("2020-05-01")

## Is the deterministic sim deterministic?
res1 <- run_sim(params,state,start_date=startdate,end_date=enddate)
res2 <- run_sim(params,state,start_date=startdate,end_date=enddate)
stopifnot(identical(res1,res2))

## Is Relative_value=1 correctly not doing anything?
time_pars <- data.frame(Date=as.Date(startdate:enddate),
                        Symbol="beta0",
                        Relative_value=1,
                        stringsAsFactors=FALSE)
res1_t <- update(res1, params_timevar=time_pars)
stopifnot(identical(c(res1),c(res1_t))) ## Use c() to drop attributes

## Does just changing mu do anything? Yes
new_mu_p <- update(params,c(mu=0.8))
new_mu <- run_sim(new_mu_p, state,start_date=startdate,end_date=enddate)
stopifnot(!identical(c(res1),c(new_mu)))

## Does putting a relative value on mu do anything? No.
tv_mu <- data.frame(Date=as.Date(startdate:enddate),
	Symbol="mu",
	Relative_value=rep(c(0.8,.7),each=length(startdate:enddate)/2),
	stringsAsFactors=FALSE
)

print(tv_mu)
res1_mu <- update(new_mu,params_timevar=tv_mu)
stopifnot(identical(c(res1),c(res1_mu)))

## Extra checks
## print(c(res1_mu))
## stopifnot(identical(c(new_mu),c(res1_mu)))

## Double-check on whether it's update; apparently not
## res2_mu <- run_sim(params, params_timevar=tv_mu,state,start_date=startdate,end_date=enddate)
## stopifnot(identical(c(res1),c(res2_mu)))

## Does putting a relative value on beta do anything? Yes.
tv_beta <- data.frame(Date=as.Date(startdate:enddate),
                         Symbol="beta0",
                         Relative_value=rep(c(1,.7),each=length(startdate:enddate)/2),
                         stringsAsFactors=FALSE)
res1_t3 <- update(res1,params_timevar=tv_beta)
print(identical(c(res1),c(res1_t3)))
