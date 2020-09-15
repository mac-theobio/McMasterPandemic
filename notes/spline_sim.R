library(McMasterPandemic)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()


fitmax <- 126 ## Copied from dependency 

params <- read_params("ICU1.csv")

X <- cbind(1,mod_ns$model[,-1])

first_date <- as.Date("2020-01-01")

dd <- first_date -1 + 1:fitmax

bb <- coef(mod_ns)

sims <- run_sim_loglin(params=params
	, extra_pars=list(time_beta=bb)
	, time_args=list(X_date=dd, X=X)
	, sim_args=list(start_date=min(dd),end_date=max(dd))
)


print(sims)

print(plot(sims$date,sims$report))

saveVars(sims, mod_ns)
