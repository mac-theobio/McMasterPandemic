library(McMasterPandemic)
library(zoo)
library(tidyverse)

## Copying example from run_sim and modifying it

##1. testing if Relative_value=1 and non-timevar run_sim are the same (done)
## MLi: This means I can always construct my own timevar for any time-varying simulation 
## Move this to testhat?!? 

params <- read_params("ICU1.csv")
paramsS <- update(params,c(proc_disp=0.1,obs_disp=100))
paramsSz <- update(paramsS, zeta=5)
state <- make_state(params=params)
startdate <- as.Date("2020-01-01")
enddate <- as.Date("2020-05-01")
time_pars <- data.frame(Date=as.Date(startdate:enddate),
                        Symbol="beta0",
                        Relative_value=1,
                        stringsAsFactors=FALSE)

## This is checking if we can get the same thing if we don't add stoch

res1 <- run_sim(params,state,start_date=startdate,end_date=enddate)
res2 <- run_sim(params,state,start_date=startdate,end_date=enddate)
stopifnot(identical(res1,res2))

## This fits a timevar dataframe where beta0 = 1 
res1_t <- update(res1, params_timevar=time_pars)

## keep only numeric values
stopifnot(identical(c(res1),c(res1_t)))

time_pars2 <- data.frame(Date=as.Date(startdate:enddate),
                        Symbol="mu",
                        Relative_value=rep(c(1,.7),each=length(startdate:enddate)/2),
                        stringsAsFactors=FALSE)
res1_t2 <- update(res1,params_timevar=time_pars2)
stopifnot(identical(c(res1),c(res1_t2)))



time_pars3 <- data.frame(Date=as.Date(startdate:enddate),
                         Symbol="beta0",
                         Relative_value=rep(c(1,.7),each=length(startdate:enddate)/2),
                         stringsAsFactors=FALSE)
res1_t3 <- update(res1,params_timevar=time_pars3)
stopifnot(identical(c(res1),c(res1_t3)))
