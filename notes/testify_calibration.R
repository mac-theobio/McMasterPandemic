library(McMasterPandemic)
library(tidyverse); theme_set(theme_bw())

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(params)


paramsw0 <- params
summary(paramsw0)

## run the *first* set of parameters from pf
print(pf[9,])
simdat <- simtestify(9)

dat <- (simdat
	%>% transmute(date
	 	, postest
		, death
		, H
	)
	%>% gather(key="var", value="value", -date)
	%>% mutate(value = round(value))
)

pp <- attr(simdat,"params")

opt_pars <- with(as.list(pp),
                 list(params=c(log_E0=log(E0)
                             , log_beta0=log(beta0)
                             ## , logit_mu = plogis(mu)
                             ## , logit_phi1 = plogis(phi1)
                             , logit_W_asymp = plogis(W_asymp)
                             , log_testing_intensity=log(0.002)
                 )
                 )
)
                               
# opt_pars <- list(params=c(log_E0=2, log_beta0=-1, logit_mu = -1
# 	, logit_phi1 = -1
# 	, logit_W_asymp = -1
# 	)
#       , log_nb_disp = c(postest=1)
# )

sim_args <- list(ratemat_args = list(testify=TRUE, testing_time="report")) 

testdat <- data.frame(Date=dat$date, intensity=testing_intensity[1])

testify_calib <- do.call(calibrate_comb,
                         c(nlist(params=paramsw0
                               , debug_plot=FALSE
                               , debug=FALSE
                               ## , use_testing = TRUE
                               , use_DEoptim = FALSE
                               ## , testing_data = testdat
                               , data = dat
                               , use_spline = FALSE
                               , opt_pars = opt_pars
                               , sim_args = sim_args
                                 ))
                         )

print(testify_calib)

dd <- predict(testify_calib
   , ensemble=FALSE
   , keep_vars=c("postest")
)

print(dd)

gg2 <- (ggplot(dd,aes(date))
   + geom_line(aes(y=value))
   + geom_point(data=dat,aes(date,value))
   + ggtitle("postest")
)

print(gg2)


