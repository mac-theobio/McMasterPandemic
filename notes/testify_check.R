library(McMasterPandemic)
library(tidyverse);theme_set(theme_bw())
library(zoo)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()


p <- read_params(makeArgs()[3])

print(p)

start <- "2020-01-01"
end <- "2020-12-25"

dateVec <- seq.Date(from=as.Date(start)
	, to = as.Date(end)
	, by = 1
)


intensity <- rep(0.002,length(dateVec))

testing_data <- data.frame(Date = dateVec
	, Symbol = "testing_intensity"
	, Relative_value = c(1, intensity[-1]/intensity[1])
)


sim_args <- list(ratemat_args = list(testing_time = "report")
	, start_date = start
	, end_date = end
	, step_args = list(testwt_scale = "sum_smooth")
	, condense_args = list(keep_all = FALSE
		, add_reports = TRUE)
	, params_timevar = testing_data)

sims <- do.call(run_sim, c(list(params=p), sim_args))

simdat <- (sims
	%>% select(date,postest,hosp,H,report)
	%>% pivot_longer(names_to = "var", -date)
	%>% filter(value>1)
)

print(simdat)

gg <- (ggplot(simdat,aes(x=date,y=value))
	+ geom_point()
	+ geom_line()
	+ facet_wrap(~var,scale="free",nrow=2)
)

print(gg)
print(gg + scale_y_log10())



