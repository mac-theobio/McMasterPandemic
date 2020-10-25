# Goal: given time series of Rt, we want to reconstruct Rt by using spline fit and predict
# Spline methods: bs() and ns() with specified degree of freedom df
library(splines)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

source("makestuff/makeRfuns.R")

makeGraphics()

set.seed(101)

# parameters 
fitmax <- 200
t <- 1:fitmax
r <- 0.04
h <- 2
R0 <- 2
c <- 2

ndf <- 7

## The susceptible effect makes R smaller, and increases through time
## The infections makes R bigger, and grows and saturates
Q <- exp(r*t)
S <- 1/(1+Q)
I <- Q/(1+Q^2)
Rt = R0*S^h*I^(-c)
logRt <- log(Rt)

Rf <- data.frame(t, Rt)
Sf <- data.frame(t
	, logss = log(Rt/R0)
)

print(plot(Sf))

# model fits
mod_bs <- lm(logss~bs(t,df=ndf)-1,data=Sf)
mod_ns <- lm(logss~ns(t,df=ndf)-1,data=Sf)

Sfpredict <- (Sf
	%>% mutate(bs_fit = predict(mod_bs)
		, ns_fit = predict(mod_ns))
	%>% gather(key = "spline_type", value="pred", -t, -logss)
)

print(Sfpredict)

gg <- (ggplot(Sfpredict, aes(t))
	+ geom_point(aes(y=logss), color="black")
	+ geom_line(aes(y=pred,color=spline_type))
	+ scale_color_manual(values=c("red","blue"))
)

print(gg)

gg2 <- (ggplot(Sfpredict, aes(t))
	+ geom_point(aes(y=exp(logss)),color="black")
	+ geom_line(aes(y=exp(pred),color=spline_type))
	+ scale_color_manual(values=c("red","blue"))
	+ ylab("Rt")
)

print(gg2)
# bs is doing better than ns
saveEnvironment()
