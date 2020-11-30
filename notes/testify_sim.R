## callArgs only works interactively and is target-dependent
callArgs <- "testwt_N.random_sims.Rout testify_sim.R testwt_N.rda testify_funs.rda random.csv"

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


params <- read_params(makeArgs()[5])
print(params)

summary(params)


## create factorial combination of parameter vectors
pf <- expand.grid(R0=R0
	, testing_intensity = testing_intensity
	, testing_type = c("constant","reduceR","increaseT","mix")
	# , testing_type = c("mix")
	, iso_t = iso_t
)

print(pf)

dateVec <- seq.Date(start,end,by=1)

timevars_constant <- data.frame(Date= rep(dateVec,2)
	, Symbol = rep(c("beta0","testing_intensity")
		,each=length(dateVec)
		)
	, Relative_value = c(seq(1,1,length.out=length(dateVec))
		, seq(1,1,length.out=length(dateVec))
	)
)

timevars_reduceR <- data.frame(Date= rep(dateVec,2)
	, Symbol = rep(c("beta0","testing_intensity")
		, each=length(dateVec)
		)
	, Relative_value = c(seq(1,0.5,length.out=length(dateVec))
		, seq(1,1,length.out=length(dateVec))
	)
)

timevars_increaseT <- data.frame(Date= rep(dateVec,2)
	, Symbol = rep(c("beta0","testing_intensity")
		, each=length(dateVec)
		)
	, Relative_value = c(seq(1,1,length.out=length(dateVec))
		, seq(1,5,length.out=length(dateVec))
	)
)

timevars_mix <- data.frame(Date = rep(dateVec,2)
	, Symbol = rep(c("beta0","testing_intensity")
		, each = length(dateVec)
		)
	, Relative_value = c(seq(1,0.5,length.out=length(dateVec))
		, seq(1,5, length.out=length(dateVec))
	)
)



## run a simulation based on parameters in the factorial frame
update_and_simulate <- function(x){
	print(x)
	## Update and fix_pars need to be together right before simulating
	paramsw0 <- update(params
		, iso_t= pf[x,"iso_t"]
		, testing_intensity = pf[x,"testing_intensity"]
		, omega = omega
	)

	if(pf[x,"testing_type"] == "constant"){
		testdat <- timevars_constant
	}
	if(pf[x,"testing_type"] == "reduceR"){
		testdat <- timevars_reduceR
	}
	if(pf[x,"testing_type"] == "increaseT"){
		testdat <- timevars_increaseT
	}
	if(pf[x,"testing_type"] == "mix"){
		testdat <- timevars_mix
	}


	paramsw0 <- fix_pars(paramsw0, target=c(R0=pf[x,"R0"],Gbar=Gbar))

	sims <- (simtestify(p=paramsw0,timevars=testdat)
	%>% mutate(R0 = pf[x,"R0"]
			, testing_intensity = pf[x,"testing_intensity"]
			, testing_type = pf[x,"testing_type"]
			, iso_t = pf[x,"iso_t"]
         )
	)
	return(sims)
}

simlist <- mclapply(1:nrow(pf),function(x)update_and_simulate(x),mc.cores = 3)

print(simlist)

simframe <- bind_rows(simlist)


print(simframe)

## What is this keep_all stuff?

if (!keep_all) {
    simdat <- (simframe
	 	  %>% group_by(R0,testing_intensity,testing_type,iso_t)
        %>% transmute(date
                    , incidence
						  , CumIncidence = cumsum(incidence)
                    , postest
                    , total_test = postest + negtest
#                    , pos_per_million = 1e6*postest/total_test
                    , positivity = postest/total_test
                    , iso_t
						  , omega
                    , testing_type
						  , Gbar
						  , testing_intensity
                      )
		  %>% ungroup()
        %>% gather(key="var",value="value",-c(date, iso_t,testing_type, testing_intensity,R0))
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
