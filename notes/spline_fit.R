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
fitmax <- 90
t <- 1:fitmax
r <- 0.04
h <- 2
R0 <- 2
c <- 2

ndf <- 5

## The susceptible effect makes R smaller, and increases through time
## The infections makes R bigger, and grows and saturates
Q <- exp(r*t)
S <- 1/(1+Q)
I <- Q/(1+Q^2)
logRt <- log(R0*S^h*I^(-c))

Rf <- data.frame(t, logRt)

# model fits
mod_bs <- lm(logRt~bs(t,df=ndf),data=Rf)
mod_ns <- lm(logRt~ns(t,df=ndf),data=Rf)

Rfpredict <- (Rf
	%>% mutate(bs_fit = predict(mod_bs)
		, ns_fit = predict(mod_ns))
	%>% gather(key = "spline_type", value="pred", -t, -logRt)
)

print(Rfpredict)

gg <- (ggplot(Rfpredict, aes(t))
	+ geom_point(aes(y=logRt), color="black")
	+ geom_line(aes(y=pred,color=spline_type))
	+ scale_color_manual(values=c("red","blue"))
)

print(gg)
# bs is doing better than ns
saveVars(Rf, Rfpredict, mod_bs, mod_ns)
