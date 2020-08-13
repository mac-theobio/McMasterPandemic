library(McMasterPandemic)
library(tidyverse); theme_set(theme_bw())

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

print(params)


paramsw0 <- params
keep_all <- FALSE
print(pf[1,])
simdat <- simtestify(1)

gg <- (ggplot(simdat,aes(date,postest))
	+ geom_point()
)

print(gg)

dat <- simdat2 %>% transmute(date, var="postest", value=postest)
opt_pars <- list(params=c(log_E0=2, log_beta0=-1, logit_mu = -1
	, logit_phi1 = -1
	, logit_W_asymp = -1
	)
	, log_nb_disp = c(postest=1)
)

sim_args <- list(ratemat_args = list(testify=TRUE, testing_time="report")) 

testdat <- data.frame(Date=dat$date, intensity=testing_intensity[1])

testify_calib <- do.call(calibrate_comb, c(nlist(params=paramsw0
	, debug_plot=FALSE
	, use_testing = TRUE
	, use_DEoptim = FALSE
	, testing_data = testdat
	, data = dat
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
   + geom_point(data=res_list[[1]]$data,aes(date,value))
   + ggtitle("postest")
)

print(gg2)


