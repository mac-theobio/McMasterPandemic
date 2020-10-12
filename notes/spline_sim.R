library(McMasterPandemic)

callArgs <- "spline_sim.Rout spline_sim.R spline_fit.rda spline.csv"
source("makestuff/makeRfuns.R")

commandEnvironments()
makeGraphics()

params <- read_params(matchFile(".csv$"))

X <- cbind(1,mod_bs$model[,-1])

first_date <- as.Date("2020-01-01")

dd <- first_date -1 + 1:fitmax

## Reconstruct the spline fit
bb <- coef(mod_bs)
Rt <-  exp(X %*% matrix(bb, ncol=1))
print(Rt)
plot(Rt)

## Time-varying betas are actually Rs, so we set base R to 1
## nullsim tests that this scaling seems to work
adj_params <- fix_pars(params, target=c(R0=Rt[[1]]))
scaled_params <- fix_pars(params, target=c(R0=1))

nullsim <- run_sim_loglin(params=adj_params
	, time_args=list(X_date=dd, X=X)
	, sim_args=list(start_date=min(dd),end_date=max(dd))
)

sim <- run_sim_loglin(params=scaled_params
	, extra_pars=list(time_beta=bb)
	, time_args=list(X_date=dd, X=X)
	, sim_args=list(start_date=min(dd),end_date=max(dd))
)

print(sim)

plot(nullsim$date,nullsim$report, log="y")
lines(sim$date,sim$report)

saveEnvironment()
