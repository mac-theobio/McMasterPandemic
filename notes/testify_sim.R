## library(McMasterPandemic)
library(tidyverse)
library(devtools)

load_all("../")

source("makestuff/makeRfuns.R")
commandEnvironments()

pp <- read_params(matchFile(".csv$"))

set.seed(0802)

## I don't think we get ensembles unless we have obs error??
nsims <- 1

simdf <- function(sim,pars){
	sim0 <- (run_sim(params = pars, end_date = "2020-06-15")
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
			params = pars, ratemat_args = list(testify=TRUE)
			, end_date = "2020-06-15"
		) 
		%>% mutate(sim = sim, type="testify")
	)
	return(bind_rows(sim0,sim_testify))
}

sim_summary <- function(nsims,pars){
	simframe <- bind_rows(lapply(seq(1,nsims),function(x)simdf(x,pars)))
	redframe <- (simframe
		%>% group_by(type,sim)
		%>% transmute(sim,type,date,incidence,report
			, newP = diff(c(0,P))
			)
		%>% ungroup()
		%>% gather(key="var",value="value",-sim,-date, -type)
		%>% group_by(date,type,var)
		%>% summarise(lwr = quantile(value,probs = 0.025,na.rm = TRUE)
				, med = quantile(value, probs = 0.5,na.rm = TRUE)
				, upr = quantile(value, probs = 0.975,na.rm = TRUE)
		)
	)
	return(redframe)
}

simdat <- sim_summary(nsims,pars=pp) %>% mutate(params="default")


## Different testing weights
ppwts <- pp
ppwts[grepl("W",names(pp))] <- 1:length(ppwts[grepl("W",names(pp))])

simdatwts <- sim_summary(nsims,pars=ppwts) %>% mutate(params="wts")


## Different positivity

ppPos <- pp
ppPos[grepl("P",names(pp))] <- rep(c(0.8,0.95),c(2,9))

simdatPos <- sim_summary(nsims,pars=ppPos) %>% mutate(params="pos")

saveVars(simdat, simdatwts, simdatPos)
