## callArgs only works interactively and is target-dependent
callArgs <- "testwt_N.sims.Rout testify_sim.R testwt_N.rda testify_funs.rda sims.csv"

library(tidyverse)
library(parallel)
library(zoo)
library(McMasterPandemic)
source("makestuff/makeRfuns.R")

## Double-sourcing will be necessary sometimes until we 
## make makeR a real package
source("makestuff/makeRfuns.R")
commandEnvironments()


if (!exists("keep_all")) keep_all <- FALSE

params <- (read_params("PHAC_testify.csv")
	%>% update(W_asymp = W_asymp
	)
)

print(params)

summary(params)


## create factorial combination of parameter vectors
pf <- expand.grid(iso_t=iso_t
	, omega=omega
	, testing_type=testing_type
	, Gbar=Gbar
	, W_asymp = W_asymp
	, testing_intensity = testing_intensity
)

print(pf)

# print(pf <- pf[c(1,2,7,8),])

datevec <- as.Date(start):as.Date(end)
testdat <- data.frame(Date = as.Date(datevec)
        , intensity =  NA
)



## run a simulation based on parameters in the factorial frame
update_and_simulate <- function(x, testdat){
	print(x)
	## Update and fix_pars need to be together right before simulating
	paramsw0 <- update(params
		, iso_t=pf[x,"iso_t"]
		, testing_intensity = pf[x,"testing_intensity"]
		, omega = pf[x,"omega"]
	)
	paramsw0 <- fix_pars(paramsw0, target=c(R0=R0,Gbar=pf[x,"Gbar"]))

	## This can be a bit cleaner
	if(pf[x,"testing_type"] == "constant"){
 		testdat$intensity <- pf[x,"testing_intensity"]
	}
	if(pf[x,"testing_type"] == "linear"){
		testdat$intensity <- seq(min_testing,max_testing,length.out =nrow(testdat))
	}
	if(pf[x,"testing_type"] == "logistic"){
		testdat$intensity <- plogis(seq(qlogis(min_testing/max_testing),qlogis(0.99),length.out = nrow(testdat)))*max_testing
	}

	
	sims <- (simtestify(p=paramsw0,testing_data=testdat)
	%>% mutate(iso_t = pf[x,"iso_t"]
     	 , Gbar = pf[x,"Gbar"]
         , testing_type = pf[x,"testing_type"]
			, testing_intensity = pf[x,"testing_intensity"]
			, omega = pf[x,"omega"]
         )
	)
	return(sims)
}

simlist <- mclapply(1:nrow(pf),function(x)update_and_simulate(x,testdat=testdat),mc.cores = 3)

print(simlist)

simframe <- bind_rows(simlist)


print(simframe)

## What is this keep_all stuff?

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
                    , testing_type
						  , Gbar
						  , testing_intensity
                      )
        %>% gather(key="var",value="value",-c(date, iso_t, omega,testing_type, Gbar, testing_intensity))
    )
} else {
    
    simdat <- (simframe
        %>% select(-c(D,X,foi,N,P))
        %>% gather(key="var",value="value",-c(date, iso_t, omega, testing_type, Gbar, testing_intensity))
        %>% separate(var,c("pref","testcat"),sep="_")
        %>% mutate_at("pref", ~ case_when(grepl("^([Hh]|IC)",.) ~ "hosp",
                                          grepl("^I[ap]",.) ~ "asymp_I",
                                          TRUE ~ .))
        %>% group_by(date, iso_t, omega, testing_type, pref, testcat, Gbar)
        %>% summarise(value=mean(value),.groups="drop")
        %>% mutate_at("pref", factor, levels=c("S","E","asymp_I","Im","Is","hosp","R"))
        %>% mutate_at("testcat", factor, levels=c("u","n","p","t"))
    )
    
}

warnings()

saveVars(simdat, params, update_and_simulate, pf)
