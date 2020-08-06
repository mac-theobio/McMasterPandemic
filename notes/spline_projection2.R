# Function to extrapolate given Rt and tplus to go beyond range of time scale of Rt

library(splines)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

source("makestuff/makeRfuns.R")
makeGraphics()

t <- 1:10
logRt <- t^2 
tplus <- 7 # time steps for extrapolation (in days)
dof <- 7 #degree of freedom of ns()

spline_extrap <- function(logRt, tplus, dof){
  fitmax <- length(logRt)
  tmax <- fitmax + tplus
  t <- 1:tmax
  
  Rf <- data.frame(t=t, logRt=ifelse(t>fitmax, NA, logRt))
  
  mod_ns <- lm(logRt~ns(t,df=dof),data=na.omit(Rf))
  co_ns <- coef(mod_ns) 
  mod_ns <- model.frame(mod_ns)
  
  pv_ns <- attr(terms(mod_ns),"predvar")
  
  bns <- eval(pv_ns[[3]])
  bns <- cbind(1, bns)
  # projection for ns
  proj_ns <- bns %*% co_ns
  
  Rf <- (Rf %>% mutate(ns_projection = proj_ns))
  return(Rf)  
}

# testing the function
out <- spline_extrap(logRt,tplus,dof)

print(gg <- ggplot(out, aes(x=t,y=ns_projection)) 
            +geom_point()
            +geom_point(aes(x=t,y=t^2),color="red")
)

