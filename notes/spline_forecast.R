library(McMasterPandemic)
library(dplyr)


# callArgs <- "spline_forecast.Rout spline_forecast.R collect_recalib.rda cachestuff/spline_calib.rda"

source("makestuff/makeRfuns.R")
commandEnvironments()


flist <- list.files(path="cachestuff/",pattern="recalib.no")

forecasting <- function(x){
	seed <- gsub("spline_recalib.noE0","",x)
	seed <- gsub(".RDS","",seed)
	set.seed(as.numeric(seed))

true_forecast <- (predict(ff,ensemble=TRUE,nsim=1,end_date = as.Date("2020-06-01"))
	%>% filter(var == "report")
	%>% mutate(seed = seed
		, type = "true"
	)
)

print(head(true_forecast))

sim_forecast <- (predict(readRDS(paste0("cachestuff/",x)),ensemble=TRUE,end_date=as.Date("2020-06-01"))
	%>% filter(var == "report")
	%>% mutate(seed = seed
		, type = "sim"
		)
)

combodat <- bind_rows(true_forecast,sim_forecast)

return(combodat)

}

forecastdat <- bind_rows(lapply(flist,function(x)forecasting(x)))

saveVars(forecastdat)




