library(McMasterPandemic)
library(tidyverse);theme_set(theme_bw())
library(zoo)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()


p <- read_params(makeArgs()[3])
p <- fix_pars(p, target = c(R0=2, Gbar=6))
p <- update(p,c(rho=1/5,c_prop=1,testing_intensity=0.002))

print(p)

start <- "2020-01-01"
end <- "2021-12-25"

dateVec <- seq.Date(from=as.Date(start)
	, to = as.Date(end)
	, by = 1
)


intensity <- rep(1,length(dateVec))

timevars <- data.frame(Date= rep(dateVec,2)
	, Symbol = rep(c("beta0","testing_intensity"),each=length(dateVec))
	, Relative_value = c(seq(1,0.6,length.out=length(dateVec))
		, rep(1,1,length.out=length(dateVec))
	)
)

sim_args <- list(ratemat_args = list(testing_time = "report")
	, start_date = start
	, end_date = end
	, step_args = list(testwt_scale = "N")
	, condense_args = list(keep_all = FALSE
		, add_reports = TRUE)
	, params_timevar = timevars)

sims <- do.call(run_sim, c(list(params=p), sim_args))

simdat <- (sims
	%>% select(date,postest,hosp,H,report,incidence,S)
	%>% pivot_longer(names_to = "var", -date)
	%>% filter(value>1)
)

print(simdat)

gg <- (ggplot(simdat,aes(x=date,y=value))
	+ geom_point()
	+ geom_line()
	+ facet_wrap(~var,scale="free_y",nrow=2)
)

print(gg)
print(gg + scale_y_log10())



