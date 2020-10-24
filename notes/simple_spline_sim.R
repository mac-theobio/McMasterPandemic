library(McMasterPandemic)
library(splines)

callArgs <- "simple_spline_sim.Rout simple_spline_sim.R  spline.csv"
source("makestuff/makeRfuns.R")

commandEnvironments()
makeGraphics()

params <- read_params(matchFile(".csv$"))

first_date <- as.Date("2020-01-01")
fitmax <- 100
dd <- first_date -1 + 1:fitmax
ndf <- 6

R0 <- 2

scaled_params <- fix_pars(params, target=c(R0=R0))

X <- bs(dd, df=ndf)


## Reconstruct the spline fit
bb <- c(-1,-0.8,-0.5,-1,-0.3,-.5)
Rt <-  R0*exp(X %*% matrix(bb, ncol=1))
plot(Rt)

scaled_params["obs_disp"] <- 50

print(scaled_params)
print(summary(scaled_params))

saveVars(bb, X,Rt, scaled_params, dd, fitmax)
