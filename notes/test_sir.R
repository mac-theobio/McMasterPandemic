## This is to setup a simple SIR model and run the simulations

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
  ## BMB: with(as.list(c(x,params)), { ... }) OR McMasterPandemic::unpack(c(x,params)) can simplify this slightly
  ## BMB: specify as parameter (don't hard-code it!); also, already specified outside the function?
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

xstart <- c(S=N0,I=1,R=0) ## initial conditions

dd <- (data.frame(date=seq(as.Date("2020-02-14"),
                                    as.Date("2020-05-31"),
                                    by="1 day"))
)


out <- as.data.frame(
  ode(
    func=sir.model,
    y=xstart,
    times=seq(1:length(dd$date)), ## BMB: use seq() *or* colon
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

params <- read_params("ICU1.csv")
# set all params to be 0 except for the corresponding params in SIR
## BMB:
params[1:length(params)] <- 0  
params["beta"] <- 1  ## BMB: this isn't used by MacPan. beta0?
params["gamma_a"] <- 1/3
params["N"] <- N0
params["E0"] <- 1

## BMB: how do you want to simplify the flow chart?
## setting alpha=0 will mean there are no asymptomatic cases,
##  so everyone flows S -> E -> Ip -> {Im, Is} -> (H/ICU/etc.)
## setting mu=0 will mean there are only severe cases, so everyone flows
##  S -> E -> Ip -> Is
## setting all the other gamma parameters to zero will mean no one ever flows
## out of the exposed class ... (or out of other compartments if they
## could ever get there)
##
## maximal simplification is probably from
params <- update(params,
                 alpha=0,  ## no asymptomatic
                 mu=1,     ## all mild cases (thus no hospital/ICU/etc.)
                 Cp=0)     ## no presymptomatic transmission
## this means we have two parallel 'exposed' (but noninf.) classes
## i.e. the distribution of exposed time is Gamma(2) (or Erlang(2)) rather
##  than exponential
##
## alternately you could make *everyone* asymptomatic (alpha=1)
## and set Ca=1 (fully infectious), in which case we would have a
## straight SEIR model (but it could get weird if we want to introduce testing)
##
## the only way to get to a true SIR model is to make sigma arbitrarily large,
##  (i.e. time in exposed class -> 0)
##  which is numerically problematic

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

sims <- run_sim_loglin(params=params,
                       ## BMB: you do need to supply a coef vector
                       ## if not (if (length(extra_pars$time_beta)) ==0)
                       ## then the spline component will be skipped
                       ## (we'll run with constant Rt)
                       extra_pars=list(time_beta=NULL),
                       time_args=list(X_date=dd$date, X=X),
                       sim_args=list(start_date=dd$date[1],
                                     end_date=tail(dd$date,1),
                                     condense_args=list(add_reports=FALSE))
)
## BMB: get object of type 'closure' is not subsettable
## replace date with dd$date
## add condense_args() to skip computing reports (on a currently empty
##  vector ...)

## right now beta0=0, time_beta=NULL, and sigma=0; so there is no
## transmission, *and* no-one ever leaves the E compartment ...

## BMB: Depending on what you are aiming to test here (time-varying
##  beta?  testing? both?, you may be able to simplify. In particular,
## if you want to run the model with no exogenous/non-autonomous
## components, you can use run_sim_range() for a greatly simplified
## interface (see ?run_sim_range)

## BMB: while we are somewhat interested in the comparison of an SIR
## with a simplified version of the model, our main goal here is to
## compare a simple version of the model with _testing machinery_
## see MacPan/notes/ratemat_vis.Rmd and notes/testing_flow_graph.R
##
## in particular you will want to allow for different ways to
## scale the weights, e.g. see R/sim_funs.R line 228, and make_test_wtsvec
## in R/testify.R ..
