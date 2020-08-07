## library(McMasterPandemic)
library(tidyverse)
library(devtools)

load_all("../")

source("makestuff/makeRfuns.R")
commandEnvironments()

pp <- read_params(matchFile(".csv$"))
pp[["iso_p"]] <- 0
pp[["N"]] <- 1.5e7 ## population of Ontario

ppw0 <- pp[!grepl("^W",names(pp))] ## Copying BMB, removing all of the regular W-parameters

start <- as.Date("2020-01-01")
end <- as.Date("2020-06-01")

print(pp)

set.seed(0807)

Wasymp <- c(0.01, 0.1,1)
isop <- c(0,0.5,0.9,1)

simframe <- data.frame()
for(i in Wasymp){
	for(j in isop){
		ppw0[["W_asymp"]] <- i
		ppw0[["iso_p"]] <- j
			sims <- (run_sim(params = ppw0, ratemat_args = list(testify=TRUE)
			, start_date = start
			, end_date = end
##			, condense_args=list(keep_all=TRUE) checkout the expanded version
		) 
	%>% mutate(W_asymp = paste0("W_asymp = ",i)
			, iso_p = paste0("iso_p = ",j)
			)
	)
	simframe <- bind_rows(simframe,sims)
	}
}

print(simframe)

simdat <- (simframe
	%>% transmute(date,incidence,postest
			, total_test = postest + negtest
			, report
			, W_asymp
			, iso_p
		)
	%>% gather(key="var",value="value",-date, -W_asymp, -iso_p)
)


saveVars(simdat)
