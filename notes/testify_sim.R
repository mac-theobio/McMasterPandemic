library(McMasterPandemic)
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()

pp <- read_params(matchFile(".csv$"))

set.seed(0802)

## I don't think we get ensembles unless we have obs error??
nsims <- 20

simdf <- function(sim){
	sim0 <- (run_sim(params = pp, end_date = "2020-06-15")
		%>% mutate(sim = sim
			, type="non_testify"
			, negtest = NA  ## Adding extra columns to match testify frame
			, N = NA
			, prostest = NA
			, P = NA 
		)
	)
	sim_testify <- (
		run_sim(
			params = pp, ratemat_args = list(testify=TRUE)
			, end_date = "2020-06-15"
		) 
		%>% mutate(sim = sim, type="testify")
	)
	return(bind_rows(sim0,sim_testify))
}

sim_summary <- function(nsims){
	simframe <- bind_rows(lapply(seq(1,nsims),simdf))
	redframe <- (simframe
		%>% select(sim,type,date,incidence,report)
		%>% gather(key="var",value="value",-sim,-date, -type)
		%>% group_by(date,type,var)
		%>% summarise(lwr = quantile(value,probs = 0.025,na.rm = TRUE)
				, med = quantile(value, probs = 0.5,na.rm = TRUE)
				, upr = quantile(value, probs = 0.975,na.rm = TRUE)
		)
	)
	return(redframe)
}

simdat <- sim_summary(nsims)

saveVars(simdat)
