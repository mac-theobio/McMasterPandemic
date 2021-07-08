library(McMasterPandemic)
library(zoo)
library(tidyverse)
library(shellpipes)

params <- read_params("ICU1.csv")

state <- make_state(params=params)
startdate <- as.Date("2020-01-01")
enddate <- as.Date("2020-05-01")

base <- run_sim(params,state,start_date=startdate,end_date=enddate)

## Does just changing mu do anything? Yes
new_mu_p <- update(params,c(mu=0.8))
new_mu <- run_sim(new_mu_p, state,start_date=startdate,end_date=enddate)
stopifnot(!identical(c(base),c(new_mu)))


## A tv substitution that works
tv_beta <- data.frame(Date=as.Date(startdate:enddate),
	Symbol="beta0",
	Relative_value=rep(c(1,.7),each=length(startdate:enddate)/2),
	stringsAsFactors=FALSE
)
ts_beta <- run_sim(params, params_timevar=tv_beta,state,start_date=startdate,end_date=enddate, verbose=TRUE)
stopifnot(!identical(c(base),c(ts_beta)))

## mu = fraction mild
tv_mu <- data.frame(Date=as.Date(startdate:enddate),
	Symbol="mu",
	Relative_value=rep(c(0.4,.9),each=length(startdate:enddate)/2),
	stringsAsFactors=FALSE
        )
ts_mu <- run_sim(params, params_timevar=tv_mu,state,start_date=startdate,end_date=enddate, verbose=TRUE)
stopifnot(!identical(c(base),c(ts_mu)))

## mu is the Fraction of symptomatic cases that are mild (or moderate)
## if relative value is < 1, then more people to I_s and hospital
## this looks correct!

combdat <- dplyr::bind_rows(base=base,ts_mu=ts_mu,
                            .id="model")
gg <- (ggplot(combdat)
    + geom_point(aes(x=date,y=H,colour=model))
    + scale_colour_manual(values=c("black","red"))
)

print(gg)

## Qs/To do:

## 1. FIX run_sim_breaks so that it accepts a "Symbol" vector in time_args (default is "beta0")
## if Symbol is length-1, replicate to match break_dates
## otherwise check that length(Symbol)==length(break_dates)
## update calibrate_comb to allow this to be done conveniently

## 2. add 'absolute' possibility to timevar  (combine with logit link to allow calibration of mu near 1 without going >1)?
## troubleshoot calibration, ideally by calibrating to the simulation above ...

## 3. if we fix spikes we can change only at 'changepoints' rather than resetting every day ... and we will be happier that
## the machinery is actually working right!

## 4. what's the best way to document ?

