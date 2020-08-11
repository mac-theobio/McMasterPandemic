# This is to setup a simple SIR model and run the simulations

require(deSolve)
require(ggplot2)
require(tidyr)

N0 <- 1e+06 #total population

sir.model <- function(t,x,params){
  ## extract the state variables
  S <- x[1]
  I <- x[2]
  R <- x[3]
  ## extract the parameters
  beta <- params["beta"] #transition rate
  gamma <- params["gamma"] #recovery rate
  N <- 1e+06 #total population
  ## the model equations
  dS.dt <- -beta*S*I/N
  dI.dt <- beta*S*I/N-gamma*I
  dR.dt <- gamma*I
  ## combine results into a single vector
  dxdt <- c(dS.dt,dI.dt,dR.dt)
  ## return result as a list!
  list(dxdt)
}

dd <- (data.frame(date=seq(as.Date("2020-02-14"),
                           as.Date("2020-03-14"),
                           by="1 day"))
)

params0 <- c(beta=1,gamma=1/3)

xstart <- c(S=N0,I=1,R=0) ## initial conditions


out <- as.data.frame(
  ode(
    func=sir.model,
    y=xstart,
    times=seq(1:length(dd$date)),
    parms=params0
  )
)

out2 <- out %>%
  pivot_longer(c(S, I, R), names_to = "compartment", values_to = "value")

ggplot(data=out2, aes(x=time,y=value,col=compartment))+geom_line()  

R0 <- params0["beta"]/params0["gamma"]

###
# step 2: using run_sim_loglin machinery with changing the params to the corresponding deterministic SIR model
# goal: Can we get the same result as of the deterministic SIR model?
library(McMasterPandemic)
library(dplyr)
devtools::load_all() ## update code if necessary


my_inf <- 100 #Inf
params <- read_params("ICU1.csv")
attributes(params)
# set all params to be 0 except for the corresponding params in SIR
params["Ca"] <-1
params[3:5] <- 0 # all relative transmissions, i.e., C_ 
  
params["beta0"] <- 1
params["alpha"] <-1 #assuming that all exposed people go to I_a

params["sigma"] <- my_inf #Inf

params["gamma_a"] <- 1/3
params["gamma_m"] <- my_inf

params["gamma_s"] <- my_inf
params["gamma_p"] <- my_inf
params["rho"] <- my_inf

params["delta"] <-0
params["mu"] <-0

params["N"] <- N0
params["E0"] <- 1





  
## run calibrate_comb() with interesting spline settings to
##  return the time-varying beta log-linear model matrix:
X <- calibrate_comb(data=dd, params=params,
                    use_spline=TRUE,
                    spline_type="ns",
                    spline_setback=14,
                    spline_extrap="constant",
                    return="X")
par(las=1)
matplot(X, type="l",lty=1, lwd=2,
        col=palette(),
        xlab="day",ylab="basis function value")

# debugonce(run_sim_loglin)
sims <- run_sim_loglin(params=params,
                       # extra_pars= NULL,
                       # time_args=list(X_date=dd$date, X=X),
                       sim_args=list(start_date=dd$date[1],end_date=tail(dd$date,1))
)

out_temp <- na.omit(as.data.frame(sims[c("date","S","I","R")]))
out3 <- out_temp %>%
  pivot_longer(c(S, I, R), names_to = "compartment", values_to = "value")

ggplot(data=na.omit(out3[1:25,]), aes(x=date,y=value,col=compartment))+geom_line()
