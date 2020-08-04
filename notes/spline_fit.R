library(splines)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)

source("makestuff/makeRfuns.R")

makeGraphics()

## We want to generate Rt first 

fitmax <- 126 # (note fitmax<=tmax)

t <- 1:fitmax
r <- 0.04
h <- 2
R0 <- 3
c <- 2

Q <- exp(r*t)
S <- 1/(1+Q)
I <- Q/(1+Q^2)
logRt <- log(R0*S^h*I^(-c))

Rf <- data.frame(t, logRt)

mod_bs <- lm(logRt~bs(t,df=7),data=Rf)
mod_ns <- lm(logRt~ns(t,df=7),data=Rf)

Rfpredict <- (Rf
	%>% mutate(bs_fit = predict(mod_bs)
		, ns_fit = predict(mod_ns)
	)
	%>% gather(key = "spline_type", value="pred", -t, -logRt)
)

print(Rfpredict)

gg <- (ggplot(Rfpredict, aes(t))
	+ geom_point(aes(y=logRt), color="black")
	+ geom_line(aes(y=pred,color=spline_type))
	+ scale_color_manual(values=c("red","blue"))
)

print(gg)

saveVars(Rf, Rfpredict, mod_bs, mod_ns)
