library(McMasterPandemic)

callArgs <- "spline_sim.Rout spline_sim.R spline_fit.rda spline.csv"
source("makestuff/makeRfuns.R")

commandEnvironments()
makeGraphics()

params <- read_params(matchFile(".csv$"))

print(head(mod_bs$model))

## The first column in mod_bs$model is logss

X <- mod_bs$model[,-1]

first_date <- as.Date("2020-01-01")

dd <- first_date -1 + 1:fitmax

## Reconstruct the spline fit
bb <- coef(mod_bs)

### This is not Rt, B(0) is 1, it is relative to R0, Rt = R0*spline_shape, if there is an intercept parameter, then Rt = R0*spline_shape/B0

Rt <-  R0*exp(X %*% matrix(bb, ncol=1))
plot(Rt)


## Everything below is experimenting, I think all we need is X, bb, and Rt from above.

saveVars(bb,X,Rt,dd, first_date, fitmax)
