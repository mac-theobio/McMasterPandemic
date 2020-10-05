library(McMasterPandemic)

callArgs <- "spline_sim.Rout spline_sim.R spline_fit.rda spline.csv"
source("makestuff/makeRfuns.R")

commandEnvironments()
makeGraphics()

R0 <- 2

params <- read_params(matchFile(".csv$"))

X <- cbind(1,mod_ns$model[,-1])

first_date <- as.Date("2020-01-01")

dd <- first_date -1 + 1:fitmax

## Reconstruct the spline fit
bb <- coef(mod_bs)
Rt <-  exp(X %*% matrix(bb, ncol=1))
print(Rt)
plot(Rt)
print(params)

## Time-vary betas are actually Rs, so we set base R to 1
params <- fix_pars(params, target=c(R0=1))
print(params)
summary(params)

nullsim <- run_sim_loglin(params=params
	, time_args=list(X_date=dd, X=X)
	, sim_args=list(start_date=min(dd),end_date=max(dd))
)

sim <- run_sim_loglin(params=params
	, extra_pars=list(time_beta=bb)
	, time_args=list(X_date=dd, X=X)
	, sim_args=list(start_date=min(dd),end_date=max(dd))
)

print(sim)

plot(nullsim$date,nullsim$report, log="y")
lines(sim$date,sim$report)

saveEnvironment()
