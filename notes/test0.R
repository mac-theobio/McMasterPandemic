library(McMasterPandemic)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
## load("forecast.RData")

co <- structure(c(beta0 = 0.337130113192467, Ca = 0.75, Cp = 1, Cs = 0.06, 
Cm = 0.92, alpha = 0.75, sigma = 0.333333333333333, gamma_p = 0.5, 
gamma_s = 0.5, gamma_m = 0.14, gamma_a = 0.2, rho = 0.08333, 
mu = 0.8, N = 57780000, E0 = 836.992371104812, iso_m = 0, iso_s = 0, 
phi1 = 0.74, phi2 = 0.28, psi1 = 0.25, psi2 = 0.2, psi3 = 0.25, 
c_prop = 0.1, c_delay_mean = 5, c_delay_cv = 0.25), class = "params_pansim")

set.seed(101)
plot(run_sim(
    params=update(co, c(obs_disp=20, proc_disp=1))
  , stoch=c(proc=TRUE,obs=FALSE)
  , stoch_start = c(proc="2020-04-12", obs="2020-04-02")
), log=TRUE)

sessionInfo()
