## callArgs only works interactively and is target-dependent
callArgs <- "basic.sims.Rout testify_sim.R basic.rda sims.csv"
callArgs <- "absolute.sims.Rout testify_sim.R absolute.rda sims.csv"

library(tidyverse)
library(parallel)
source("makestuff/makeRfuns.R")
print(commandEnvironments())

## Double-sourcing will be necessary sometimes until we 
## make makeR a real package
source("makestuff/makeRfuns.R")

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

print(W_asymp)
print(iso_t)
print(testing_intensity)

## create factorial combination of parameter vectors
pf <- expand.grid(W_asymp=W_asymp,iso_t=iso_t,testing_intensity=testing_intensity)
print(pf)

## run a simulation based on parameters in the factorial frame
simtestify <- function(x){
	cat(pf[x,1],pf[x,2],pf[x,3],"\n")
	paramsw0 <- update(paramsw0
		, W_asymp=pf[x,1]
		, iso_t = pf[x,2]
		, testing_intensity=pf[x,3]
	)
        sims <- (run_sim(params = paramsw0, ratemat_args = list(testify=TRUE)
                       , start_date = start
                       , end_date = end
                       , use_ode = use_ode
                       , step_args = list(testwt_scale=testwt_scale)
                       , condense_args=list(keep_all=keep_all, add_reports=!keep_all) ## checkout the expanded version
                         )
            %>% mutate(W_asymp = pf[x,1]
                     , iso_t = pf[x,2]
                     , testing_intensity=pf[x,3]
                       )
	)
	return(sims)
}

simlist <- mclapply(1:nrow(pf),function(x)simtestify(x),mc.cores = 3)
simframe <- bind_rows(simlist)

print(simframe)

if (!keep_all) {
    simdat <- (simframe
        %>% transmute(date
                    , incidence
                    , postest
                    , total_test = postest + negtest
                    , pos_per_million = 1e6*postest/total_test
                    , report
                    , W_asymp
                    , iso_t
                    , testing_intensity
                      )
        %>% gather(key="var",value="value",-c(date, W_asymp, iso_t, testing_intensity))
    )
} else {
    
    simdat <- (simframe
        %>% select(-c(D,X,foi,N,P))
        %>% gather(key="var",value="value",-c(date, W_asymp, iso_t, testing_intensity))
        %>% separate(var,c("pref","testcat"),sep="_")
        %>% mutate_at("pref", ~ case_when(grepl("^([Hh]|IC)",.) ~ "hosp",
                                          grepl("^I[ap]",.) ~ "asymp_I",
                                          TRUE ~ .))
        %>% group_by(date,W_asymp, iso_t, testing_intensity, pref, testcat)
        %>% summarise(value=mean(value),.groups="drop")
        %>% mutate_at("pref", factor, levels=c("S","E","asymp_I","Im","Is","hosp","R"))
        %>% mutate_at("testcat", factor, levels=c("u","n","p","t"))
    )
    
}

warnings()

saveVars(simdat, params, simtestify, pf)
