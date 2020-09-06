library(McMasterPandemic)
library(tidyverse); theme_set(theme_bw())
callArgs <- "basic.rda basic.sims.rda"

source("makestuff/makeRfuns.R")
print(commandEnvironments())
if (!interactive()) makeGraphics()

print(params)


paramsw0 <- params
summary(paramsw0)

## run combinations from pf
simdat <- simtestify(1)

dat <- (simdat
	%>% transmute(date
	 	, postest
	, death
	, hosp
	, report
	, incidence
	, H
	)
	%>% gather(key="var", value="value", -date)
	%>% mutate(value = round(value))
	# %>% filter(date > as.Date("2020-02-25"))
)

ggdat <- (ggplot(dat,aes(date,value,color=var)) 
   + geom_point()
   + scale_y_log10()
)
ggdat

ggdat2 <- (ggdat
   + facet_wrap(~var,ncol=1)   
)

ggdat2

pp <- attr(simdat,"params")

opt_pars <- with(as.list(pp),
                 list(params=c(log_beta0=log(beta0)
                             # , log_E0=log(E0)
                             ## , logit_mu = plogis(mu)
                             ## , logit_phi1 = plogis(phi1)
                            # , logit_W_asymp = plogis(W_asymp)
                             # , log_testing_intensity=log(testing_intensity)
                 )
                 )
)
                               

sim_args <- list(ratemat_args = list(testify=TRUE, testing_time="report")) 

testdat <- data.frame(Date=dat$date, intensity=testing_intensity[1])

testify_calib <- do.call(calibrate_comb,
                         c(nlist(params=pp
                               , debug_plot=TRUE
                               , debug=FALSE
                               ## , use_testing = TRUE
                               , use_DEoptim = FALSE
                               ## , testing_data = testdat
                               , data = dat
                               , use_spline = FALSE
                               , opt_pars = opt_pars
                               , sim_args = sim_args
                               # , skip.hessian=TRUE
                               , maxit = 1000
                               # , start_date = "2020-01-27"
                                 ))
                         )

print(testify_calib)

dd <- predict(testify_calib
   , ensemble=FALSE
   , keep_vars=c("postest","death","H")
)

print(dd)

gg2 <- (ggplot(dd,aes(date))
   + geom_line(aes(y=value))
   + geom_point(data=dat,aes(date,value))
   + facet_wrap(~var,scale="free",ncol=2)
   + ggtitle("postest")
)

print(gg2)
