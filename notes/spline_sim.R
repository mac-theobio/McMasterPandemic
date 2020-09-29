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

## Construct a spline
bb <- coef(mod_ns)
bt <-  exp(X %*% matrix(bb, ncol=1))
print(bt)
plot(bt)

print(params)
params <- update(params, beta0=bt[[1]])
params <- fix_pars(params, target=c(R0=R0))
summary(params)

sims <- run_sim_loglin(params=params
	, extra_pars=list(time_beta=bb)
	, time_args=list(X_date=dd, X=X)
	, sim_args=list(start_date=min(dd),end_date=max(dd))
)

print(sims)

print(plot(sims$date,sims$report, log="y"))

saveEnvironment()
