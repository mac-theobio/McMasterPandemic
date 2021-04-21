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

tv_mu <- data.frame(Date=as.Date(startdate:enddate),
	Symbol="mu",
	Relative_value=rep(c(0.8,.7),each=length(startdate:enddate)/2),
	stringsAsFactors=FALSE
)
ts_mu <- run_sim(params, params_timevar=tv_mu,state,start_date=startdate,end_date=enddate, verbose=TRUE)
stopifnot(!identical(c(base),c(ts_mu)))

## mu is the Fraction of symptomatic cases that are mild (or moderate)
## if relative value is < 1, then more people to I_s and hospital
## this looks correct!

gg <- (ggplot()
       + geom_point(data=base,aes(x=date,y=H),color="black")
       + geom_point(data=ts_mu,aes(x=date,y=H),color="red")
)

print(gg)
=======
stopifnot(!identical(c(base),c(ts_mu)))

params_timevar2

ts_mu <- run_sim(params, params_timevar=tv_mu,state,start_date=startdate,end_date=enddate, verbose=TRUE)
## what if e.g. mu ends up >1 during calibration?
## mu=min(mu) ?
## extend time_pars to specify the scale


## calibration?
>>>>>>> Stashed changes
