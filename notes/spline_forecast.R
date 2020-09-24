library(McMasterPandemic)
library(dplyr)


# callArgs <- "spline_forecast.Rout spline_forecast.R collect_recalib.rda cachestuff/spline_calib.rda"

source("makestuff/makeRfuns.R")
commandEnvironments()


flist <- list.files(path="cachestuff/",pattern="spline_recalib")

forecasting <- function(x){
	
	modlist <- readRDS(paste0("cachestuff/",x))
	truedat <- (modlist$full_dat
		%>% filter(var == "report")
		%>% mutate(seed = x
			, value = ifelse(is.na(value),NA,value)
			, type = "true"
			, lwr = value
			, upr = value
			)
	)
	
	sim_forecast <- (predict(readRDS(paste0("cachestuff/",x))$fit,ensemble=TRUE,end_date=max(truedat$date)
		, stoch=c(proc=FALSE,obs=TRUE))
	%>% filter(var == "report")
	%>% mutate(seed = x
		, type = "sim"
		)
	%>% select(-vtype)
)
combodat <- bind_rows(truedat,sim_forecast)

return(combodat)

}

forecastdat <- bind_rows(lapply(flist,function(x)forecasting(x)))

saveVars(forecastdat)




