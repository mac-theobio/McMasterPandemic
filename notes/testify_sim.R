## callArgs only works interactively and is target-dependent
callArgs <- "testwt_N.sims.Rout testify_sim.R testwt_N.rda testify_funs.rda sims.csv"

library(tidyverse)
library(parallel)
source("makestuff/makeRfuns.R")
print(commandEnvironments())

## Double-sourcing will be necessary sometimes until we 
## make makeR a real package
source("makestuff/makeRfuns.R")
commandEnvironments()
## WHICH of these do you want today?
## library(devtools); load_all("../")
library("McMasterPandemic")

if (!exists("keep_all")) keep_all <- FALSE

fn <- matchFile(".csv$")

params <- (read_params(fn)
    %>% fix_pars(target=c(R0=R0, Gbar=Gbar))
    %>% update(
            N=pop
        )
)

paramsw0 <- params[!grepl("^W",names(params))] ## removing all of the regular W-parameters
class(paramsw0) <- "params_pansim"

print(paramsw0)
summary(paramsw0)

print(iso_t)
print(omega)
print(testing_intensity)


## create factorial combination of parameter vectors
pf <- expand.grid(iso_t=iso_t,omega=omega,testing_intensity=testing_intensity)
print(pf)

## run a simulation based on parameters in the factorial frame
update_and_simulate <- function(x){
	cat(pf[x,"iso_t"],pf[x,"omega"],pf[x,"testing_intensity"],"\n")
	paramsw0 <- update(params
		, iso_t=pf[x,"iso_t"]
		, omega = pf[x,"omega"]
		, testing_intensity=pf[x,"testing_intensity"]
	)

	sims <- (simulate_testify_sim(p=paramsw0)
		%>% mutate(iso_t = pf[x,"iso_t"]
      	 , omega = pf[x,"omega"]
          , testing_intensity = pf[x,"testing_intensity"]
          )
	)
	return(sims)
}


simlist <- mclapply(1:nrow(pf),function(x)update_and_simulate(x),mc.cores = 3)
simframe <- bind_rows(simlist)

if (!keep_all) {
    simdat <- (simframe
        %>% transmute(date
                    , incidence
                    , postest
                    , total_test = postest + negtest
                    , pos_per_million = 1e6*postest/total_test
                    , report
                    , iso_t
						  , omega
                    , testing_intensity
                      )
        %>% gather(key="var",value="value",-c(date, iso_t, omega,testing_intensity))
    )
} else {
    
    simdat <- (simframe
        %>% select(-c(D,X,foi,N,P))
        %>% gather(key="var",value="value",-c(date, iso_t, omega, testing_intensity))
        %>% separate(var,c("pref","testcat"),sep="_")
        %>% mutate_at("pref", ~ case_when(grepl("^([Hh]|IC)",.) ~ "hosp",
                                          grepl("^I[ap]",.) ~ "asymp_I",
                                          TRUE ~ .))
        %>% group_by(date, iso_t, omega, testing_intensity, pref, testcat)
        %>% summarise(value=mean(value),.groups="drop")
        %>% mutate_at("pref", factor, levels=c("S","E","asymp_I","Im","Is","hosp","R"))
        %>% mutate_at("testcat", factor, levels=c("u","n","p","t"))
    )
    
}

warnings()

saveVars(simdat, params, simtestify, pf)
