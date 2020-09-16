library(McMasterPandemic)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

R0 <- 2

fitmax <- 126 ## Copied from dependency 

params <- read_params(matchFile(".csv$"))

X <- cbind(1,mod_ns$model[,-1])

first_date <- as.Date("2020-01-01")

dd <- first_date -1 + 1:fitmax

bb <- coef(mod_ns)

bt <-  exp(X %*% matrix(bb, ncol=1))

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


## Check if MLi's X is the same as BMB's way to get X
X2 <- calibrate_comb(data=sims, params=params
	, use_spline=TRUE
	, spline_type="ns"
	, spline_df = 7
	, spline_extrap="constant"
	, return="X"
)

X2 <- cbind(1,X2)

print(X2)

print(dim(X2))

print(all.equal(X,X2))

saveVars(sims, mod_ns, params, X, bt,bb)
