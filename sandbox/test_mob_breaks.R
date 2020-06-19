## run from head directory
library(McMasterPandemic)
load("ontario/ontario_clean.RData")
head(comb_sub) ## from package
pars <- read_params("ICU1.csv")


## beta0 acts as the 'true' intercept term in each of these cases.
## the model matrix shown here sweeps up what's left.

## mobility: 1 break, single intercept, piecewise
X1 <- calibrate_comb(ont_recent_nt,
               params=pars,
               use_mobility=TRUE,
               mob_data=comb_sub,
               use_spline=FALSE,
               mob_breaks="2020-04-15",
               return_X=TRUE)
matplot(X1,type="l",lty=1)

## 1 break, 2 intercepts (= 1 in model matrix), piecewise
## param 1 = period 2 change in intercept
## param 2 = period 1 mobility power
## param 3 = period 2 mobility power
X2 <- calibrate_comb(ont_recent_nt,
               params=pars,
               use_mobility=TRUE,
               mob_data=comb_sub,
               use_spline=FALSE,
               mob_breaks="2020-04-15",
               mob_breaks_int=TRUE,
               return_X=TRUE)
matplot(X2,type="l",lty=1)

## 1 break, single intercept, logistic
## param 1: period-1 mobility (drops down as mobility decreases,
##  then increases smoothly back to 1 as it 'shuts off'
## param 2: period-2 mobility (at zero initially, then drops down
##  as it takes effect
X3 <- calibrate_comb(ont_recent_nt,
               params=pars,
               use_mobility=TRUE,
               mob_data=comb_sub,
               use_spline=FALSE,
               mob_breaks="2020-04-15",
               mob_logist_scale=3,
               return_X=TRUE)
matplot(X3,type="l",lty=1)


## param 1: indicator variable for period 2, but takes effect smoothly
## param 2, 3: equiv to params 1/2 in previous example
X4 <- calibrate_comb(ont_recent_nt,
               params=pars,
               use_mobility=TRUE,
               mob_data=comb_sub,
               use_spline=FALSE,
               mob_breaks="2020-04-15",
               mob_logist_scale=3,
               mob_breaks_int=TRUE,
               return_X=TRUE)

matplot(X4,type="l",lty=1)


