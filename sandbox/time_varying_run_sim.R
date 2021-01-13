## Time varying bt using run_sim

library(McMasterPandemic)
library(zoo)
library(tidyverse)

## Copying example from run_sim and modifying it

##1. testing if Relative_value=1 and non-timevar run_sim are the same

params <- read_params("ICU1.csv")
pp <- fix_pars(params, target = c(R0 = 1.3, Gbar=6))
state <- make_state(params=pp)

startdate <- as.Date("2020-01-01")
enddate <- as.Date("2020-10-01")

## This is checking if we can get the same thing if we don't add stoch

sim0 <- run_sim(pp,state,start_date=startdate,end_date=enddate)

gg0 <- (ggplot(sim0,aes(x=date))
	+ geom_point(aes(y=incidence))
)

print(gg0)


## We want a dataframe with time varying relative beta through time
## If relative beta is constant though time, it should give back the same trajectory
time_pars <- data.frame(Date=as.Date(startdate:enddate)
	, Symbol="beta0"
	, Relative_value=1
)
	# , stringsAsFactors=FALSE)


## This fits a timevar dataframe where beta0 = 1 
sim0_t <- update(sim0, params_timevar=time_pars)

print(gg0
	+ geom_point(data=sim0_t, aes(x=date,y=incidence), color="red")
)

## Now we want relative beta to drop by a factor of 2 linearly starting from July 1st to Oct 1st

lockdown <- as.Date("2020-07-01")

time_pars2 <- data.frame(Date=as.Date(startdate:enddate)
	, Symbol="beta0"
	, Relative_value= c(rep(1, length(startdate:lockdown)-1)
		, seq(1,0.5,length.out = length(lockdown:enddate))
	)
)
print(time_pars2)

sim0_t_reduce <- update(sim0, params_timevar=time_pars2)

gg_rel_beta <- (ggplot(time_pars, aes(x=Date))
	+ geom_point(aes(y=Relative_value))
	+ geom_point(data=time_pars2, aes(x=Date, y=Relative_value), color="red")
)

print(gg_rel_beta)

print(gg0
		+ geom_point(data=sim0_t_reduce, aes(x=date,y=incidence), color="red")
)
