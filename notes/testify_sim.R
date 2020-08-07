## library(McMasterPandemic)
library(tidyverse)
library(devtools)

load_all("../")

source("makestuff/makeRfuns.R")
commandEnvironments()

pp <- read_params(matchFile(".csv$"))
pp[["iso_p"]] <- 0
pp[["testing_intensity"]] <- 0.002
pp[["N"]] <- 1.5e7

start <- as.Date("2020-01-01")
end <- as.Date("2020-06-01")

print(pp)

set.seed(0802)

## I don't think we get ensembles unless we have obs error??
nsims <- 1

simdf <- function(sim,pars,untestify){
	sim0 <- (run_sim(params = pars
		, start_date = start
		, end_date = end)
		%>% mutate(sim = sim
			, type="non_testify"
			, negtest = NA  ## Adding extra columns to match testify frame
			, N = NA
			, postest = report ## they should mean the same thing in non-testify world
			, P = NA 
		)
	)
	if(!untestify){
	sim0 <- data.frame()
	}
	sim_testify <- (
		run_sim(
			params = pars, ratemat_args = list(testify=TRUE)
			, start_date = start
			, end_date = end
		) 
		%>% mutate(sim = sim, type="testify")
	)
	return(bind_rows(sim0,sim_testify))
}

sim_summary <- function(nsims,pars,untestify=FALSE){
	simframe <- bind_rows(lapply(seq(1,nsims),function(x)simdf(x,pars,untestify)))
	redframe <- (simframe
		%>% group_by(type,sim)
		%>% transmute(sim,type,date,incidence,postest
			, total_test = postest + negtest
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
ppPos[grepl("P",names(pp))] <- rep(c(0,0.5),c(1,10))


simdatPos <- sim_summary(nsims,pars=ppPos) %>% mutate(params="pos")

simcombo <- bind_rows(simdat, simdatwts, simdatPos)

saveVars(simcombo)
