## run from head of repo
source("sandbox/tmb-sandbox/mre_tmb_obj_fun_nans.R", echo=TRUE)
obj_fun$fn()  ## NaN  (bad pars are cached from last time)
obj_fun$par   ## original default/'good' pars
obj_fun$fn(obj_fun$par) 
obj_fun$fn()
bad_pars

## try out 
test_par <- function(lastpar) {
  p <- obj_fun$par
  p[length(p)] <- lastpar
  obj_fun$fn(p)
}
bad_vec <- seq(1,5, by = 0.01)
bad_vals <- sapply(bad_vec, test_par)
plot(bad_vec, bad_vals)

## there are good ways to catch/store warnings and
##  errors in R but I do it too rarely/am too lazy to
##  work it out here.

## original values:
s1 <- obj_fun$simulate(obj_fun$par)
## no warnings

m <- function(x) {
  matplot(x$simulation_history, type = "l", lty = 1, log = "y")
}
m(s1)  ## already have zero values 
##  (e.g. starting values before epidemic takes off?)


s2 <- obj_fun$simulate(c(obj_fun$par[1:7],0.8))
## negative rate marix element warnings start between
##  last-element == 0.5 and 0.8

s3 <- try(obj_fun$simulate(c(obj_fun$par[1:7],2)))
## 'negative simulation value' error