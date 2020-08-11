# This is to setup a simple SIR model and run the simulations

require(deSolve)
require(ggplot2)
require(tidyr)

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

params0 <- c(beta=1,gamma=1/3)

xstart <- c(S=1e+06,I=1,R=0) ## initial conditions

dd <- (data.frame(date=seq(as.Date("2020-02-14"),
                                    as.Date("2020-05-31"),
                                    by="1 day"))
)


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

R0 <- params0["beta"]/params["gamma"]

###





